#!/usr/bin/env python
import pyranges as pr
import pandas as pd
import numpy as np
from easyterm import command_line_options

try:
    import pyranges as pr
    import pyfaidx

except:
    # printerr
    raise notracebackException(
        "ERROR pyranges & pyfaidx must be installed to use selenoprofiles assess! Try this command:\npip install pyranges\n\nor follow instructions at https://pyranges.readthedocs.io/en/latest/installation.html"
    )


help_msg = """This program obtains two tables containing the annotations for all selenoprotein predicted genes from Selenoprofiles of the input genome.

### Input/Output:
-s Specify the input Selenoprofiles file in GTF or GFF format.
-e Specify the input genome file in GFF3 or GTF format. 
-f Specify the input genome Fasta file.
-o Specify the name of the output csv table containing annotation for each Ensembl transcript.
-agg Specify the name of the output csv aggregate table.
-cs Specify the name of selenoprofiles column which will be taken as ID to work with. Default is transcript_id.
-cg Specify the name of input genome's column which will be taken as ID to work with. Default is ID.
-stop Specify if stop codons are removed and how. Three options available: 'auto', 'all' and 'no'. Default is 'auto'.
      'auto': In case some genome transcripts have stop codons and some not. Searches for transcripts with stop codons and removes them.
      'all': The last three positions from each genome transcript are removed, without performing any search.
      'no': Assumes there are no stop codons in genome transcripts, so they aren't removed.

Note that if any the input or output files are not specified the script will crash

### Script description:
This script has been initially designed to work with Selenoprofiles GTF file format and genome GFF3 file format. If another type of formats are given to these input options, 
-cs and -cg must be provided to run correctly. It works only with coding sequences (CDS), so 3' and 5' UTRs are discarded. 
This script assumes each Selenoprofiles ID contains a described Selenocysteine position (Genomic intervals, strand...). This script compares genomic intervals between Selenoprofiles and genomes. 
Thus, annotations are based on overlaps between Selenoprofiles predicted genes and genome transcripts.

### Requirements:

Pyranges version: 0.0.120
Pandas version: 1.1.5
Numpy version: 1.19.2
Easyterm version: 0.7.2

See https://github.com/maxtico/assess_annotation for more information.
"""

def_opt = {
    "s": "",
    "e": "",
    "f": "",
    "o": "run_multiple.csv",
    "agg": "run_aggregate.csv",
    "cs": "transcript_id",
    "cg": "ID",
    "stop": "auto",
    "cmd": "assess",
}

def calculate_frame(self, by):
    """Assess the frame of each genomic interval, assuming all are coding sequences.

    The input Pyranges contains an added "Frame" column, which determines the base of the CDS that is the first base of a codon.
    Resulting values are in range between "0" and "2" included. "0" indicates that the first base of the CDS is the first base of a codon,
    "1" indicates the second base and "2" indicates the third base of the CDS.
    "by" argument allows to calculate the frame for each transcript.

       Parameters
       ----------
       by : str or list of str

          Column(s) to group by to calculate the frame for each transcript.

       Returns
       -------
       PyRanges
           The function does not return anything because it adds a "Frame" column inplace.


       Examples
       --------
       >>> p= pr.from_dict({"Chromosome": [1,1,1,2,2],
       ...                  "Strand": ["+","+","+","-","-"],
       ...                  "Start": [1,31,52,101,201],
       ...                  "End": [10,45,90,130,218],
       ...                  "transcript_id": ["t1","t1","t1","t2","t2"] })
       >>> p
       +--------------+--------------+-----------+-----------+-----------------+
       |   Chromosome | Strand       |     Start |       End | transcript_id   |
       |   (category) | (category)   |   (int32) |   (int32) | (object)        |
       |--------------+--------------+-----------+-----------+-----------------|
       |            1 | +            |         1 |        10 | t1              |
       |            1 | +            |        31 |        45 | t1              |
       |            1 | +            |        52 |        90 | t1              |
       |            2 | -            |       101 |       130 | t2              |
       |            2 | -            |       201 |       218 | t2              |
       +--------------+--------------+-----------+-----------+-----------------+
       Stranded PyRanges object has 5 rows and 5 columns from 2 chromosomes.
       For printing, the PyRanges was sorted on Chromosome and Strand.

       >>> p.calculate_frame(by='transcript_id')
       +--------------+--------------+-----------+-----------+-----------------+-----------+
       |   Chromosome | Strand       |     Start |       End | transcript_id   |     Frame |
       |   (category) | (category)   |   (int32) |   (int32) | (object)        |   (int32) |
       |--------------+--------------+-----------+-----------+-----------------+-----------|
       |            1 | +            |         1 |        10 | t1              |         0 |
       |            1 | +            |        31 |        45 | t1              |         0 |
       |            1 | +            |        52 |        90 | t1              |         2 |
       |            2 | -            |       101 |       130 | t2              |         2 |
       |            2 | -            |       201 |       218 | t2              |         0 |
       +--------------+--------------+-----------+-----------+-----------------+-----------+
       Stranded PyRanges object has 5 rows and 6 columns from 2 chromosomes.
       For printing, the PyRanges was sorted on Chromosome and Strand.

    """
    # Column to save the initial index
    self.__index__ = np.arange(len(self))

    # Filtering for desired columns
    if type(by) == list:
        l = by
    else:
        l = by.split()
    sorted_p = self[["Strand", "__index__"] + l]

    # Sorting by 5' (Intervals on + are sorted by ascending order and - are sorted by descending order)
    sorted_p = sorted_p.sort(by="5")

    # Creating a column saving the length for the intervals (for selenoprofiles and ensembl)
    sorted_p.__length__ = sorted_p.End - sorted_p.Start

    # Creating a column saving the cummulative length for the intervals
    for k, df in sorted_p:
        sorted_p.dfs[k]["__cumsum__"] = df.groupby(by=by).__length__.cumsum()

    # Creating a frame column
    sorted_p.Frame = (sorted_p.__cumsum__ - sorted_p.__length__) % 3

    # Appending the Frame of sorted_p by the index of p
    sorted_p = sorted_p.apply(lambda df: df.sort_values(by="__index__"))

    self.Frame = sorted_p.Frame

    # Drop __index__ column
    self.apply(lambda df: df.drop(["__index__"], axis=1, inplace=True))


def ass_ann(df,sec):
   """Assesses the genome annotation from selenoprofiles predictions.
   It takes as input a df which has normal and _ens columns (Normal for selenoprotein predictions and _ens for genome transcripts.
   "sec" option is needed always, as it contains the selenocysteines of the df for each predicted transcript.
   """
   # Assert control condition
   assert len(df["Strand"].unique()) == 1, "Transcript has more than one strand!"

   # Get the selenocysteines which match with our transcripts
   Sec = sec[sec['transcript_id'].isin(df['transcript_id'])]  # Add ensembl id

   # Get the selenocysteines for our specific transcript_id_ens
   Sp_sec = Sec[Sec["transcript_id_ens"].isin(df["transcript_id_ens"])]

   if (df["Strand"] == df["Strand_ens"]).all():
       # Condition for missing transcripts
       if (df["transcript_id_ens"] == "-1").all():
           df["Type_annotation"] = "Missing"

       # Condition for out of frame transcripts
       elif (
           (df["Feature"] == "CDS")
           & (df["Frame_genome"] != df["Frame_genome_ens"])
       ).any():
           df["Type_annotation"] = "Out of frame"

       # Condition for well annotated transcripts
       elif (
           not Sp_sec.empty
           and (
               (Sp_sec["Feature"] == "Selenocysteine")
               & (Sp_sec["transcript_id_ens"] != "-1")
               & (Sp_sec["Overlap"] == 3)
           ).all()):
           df["Type_annotation"] = "Well annotated"

       elif (
           (df["Feature"] == "Selenocysteine")
           & (df["transcript_id_ens"] != "-1")
           & (df["Overlap"] != 3)
       ).any():
           df["Type_annotation"] = "Spliced"

       # Check for different types of misannotations
       else:
           # Creating a temporary dataframe to store the df lines with the selenocysteine line/s
           temp_df = pd.merge(Sec, df, on='transcript_id', suffixes=("_Sec", "_all"))

           # Checking for + strand
           if (temp_df["Strand_Sec"] == "+").all():
                # For stop codon cases
                if (temp_df["Start_Sec"] == temp_df["End_ens_all"]).any() and not (temp_df["End_Sec"]<temp_df["Start_ens_all"]).any():
                   df["Type_annotation"] = "Stop codon"

                # For skipped cases
                elif (temp_df["Start_Sec"] > temp_df["End_ens_all"]).any() and (
                   temp_df["End_Sec"] <= temp_df["Start_ens_all"]
                ).any():
                   df["Type_annotation"] = "Skipped"

                # For upstream cases
                elif (temp_df["Start_Sec"] > temp_df["End_ens_all"]).any():
                   df["Type_annotation"] = "Upstream"

                # For downstream cases
                elif (temp_df["End_Sec"] <= temp_df["Start_ens_all"]).any():
                   df["Type_annotation"] = "Downstream"

                # For other cases
                else:
                   df["Type_annotation"] = "Other"

           # Checking for - strand
           elif (temp_df["Strand_Sec"] == "-").all():
                if (temp_df["End_Sec"] == temp_df["Start_ens_all"]).any() and not (temp_df["Start_Sec"]<temp_df["End_ens_all"]).any():
                   df["Type_annotation"] = "Stop codon"

                # For skipped cases
                elif (temp_df["Start_Sec"] >= temp_df["End_ens_all"]).any() and (
                   temp_df["End_Sec"] < temp_df["Start_ens_all"]
                ).any():
                   df["Type_annotation"] = "Skipped"

                # For upstream cases
                elif (temp_df["End_Sec"] < temp_df["Start_ens_all"]).any():
                   df["Type_annotation"] = "Upstream"

                # For downstream cases
                elif (temp_df["Start_Sec"] >= temp_df["End_ens_all"]).any():
                   df["Type_annotation"] = "Downstream"

                # For other cases
                else:
                   df["Type_annotation"] = "Other"
   else:
      df["Type_annotation"] = "Other"

   return df


def main(args={}):
    if not args:
        opt = command_line_options(def_opt, help_msg)

    else:
        opt = args

    if (opt["s"] == None) | (opt["e"] == None) | (opt["f"] == None):
        raise NoTracebackError("ERROR options -s/-e/-f are compulsory!")


    print("*******************************")
    print("Reading input files...")
    print("*******************************")

    # File inputs for ensembl ans selenoprofiles
    if ".gff3" or ".gff" in opt["e"]:
        genome_pyr = pr.read_gff3(opt["e"])
    else:
        genome_pyr = pr.read_gtf(opt["e"])

    if ".gtf" in opt["s"]:
        seleno_pyr = pr.read_gtf(opt["s"])
    else:
        seleno_pyr = pr.read_gff3(opt["s"])

    # Saving all the exons from ensembl genome
    CDS_df = genome_pyr[genome_pyr.Feature.str.startswith("CDS")]
    CDS_df = CDS_df.as_df()
    CDS_df.rename(columns={opt["cg"]: "transcript_id_ens"}, inplace=True)
    CDS_df["Strand"] = CDS_df["Strand"].astype("string")
    CDS_df = pr.PyRanges(CDS_df, int64=True)

    ##REMOVING STOP CODONS
    if opt["stop"] == "auto":
        print("\n*******************************")
        print("Removing stop codons...")
        print("*******************************")

        # Getting the stop codons from the transcripts
        last_codons = CDS_df.spliced_subsequence(
            -3, by="transcript_id_ens"
        )  # Aqui selecciones els intervals de les ulimes tres posicions del gff3

        # We get the spliced sequence for each transcript
        last_codons_seq = pr.get_transcript_sequence(
            last_codons, path=opt["f"], group_by="transcript_id_ens"
        )

        # Create a column to know which are stop codons and which not
        last_codons_seq.Sequence = last_codons_seq.Sequence.str.upper()
        has_stop = last_codons_seq.Sequence.isin({"TGA", "TAA", "TAG"})

        # Save the IDs of the sequences which contain stop codons
        has_stop_ids = last_codons_seq[has_stop].transcript_id_ens

        # Removing those transcripts which contain stop codons
        nc_df = CDS_df[
            CDS_df.transcript_id_ens.isin(has_stop_ids)
        ]  # Transcripts WITH stop codons
        sc_df = CDS_df[
            ~(CDS_df.transcript_id_ens.isin(has_stop_ids))
        ]  # Transcripts with NO stop codons

        # Removing the stop codons from the transcripts
        no_stop_df = nc_df.spliced_subsequence(0, -3, by="transcript_id_ens")

        no_stop_df = no_stop_df.as_df()
        sc_df = sc_df.as_df()
        # Concatenation of the two dataframes in a single ENSEMBL GFF in order to apply our conditions
        ensembl_nosc = pd.concat([sc_df, no_stop_df], ignore_index=True)
        ensembl_nosc = pr.PyRanges(ensembl_nosc)

    elif opt["stop"] == "all":
        ensembl_nosc = CDS_df.spliced_subsequence(-3, by="transcript_id_ens")

    else:
        ensembl_nosc = CDS_df
    # PREPROCESSING

    print("\n*******************************")
    print("Preprocessing the files...")
    print("*******************************")

    # To work only with desired colums where ID is ensembl transcript ID and transcript_id is from Selenoprofile transcripts
    seleno_pyr=seleno_pyr.as_df()
    seleno_pyr.rename(columns={opt["cs"]: "transcript_id"}, inplace=True)
    seleno_pyr = pr.PyRanges(seleno_pyr, int64=True)
    sel_cor = seleno_pyr[["Source", "Strand", "Feature", "transcript_id"]]
    ens_cor = ensembl_nosc[["Source", "Strand", "Feature", "transcript_id_ens"]]

    # Converting ensembl Feature column into category type
    ens_cor = ens_cor.as_df()
    ens_cor["Strand"] = ens_cor["Strand"].astype("string")
    ens_cor["Feature"] = ens_cor["Feature"].astype("category")

    # Renaming the transcript_id names from both dataframes
    ens_cor["transcript_id_ens"] = ens_cor["transcript_id_ens"].str.replace("CDS:", "")
    ens_cor = pr.PyRanges(ens_cor, int64=True)
    sel_cor = sel_cor.as_df()
    secs = sel_cor[sel_cor['transcript_id'].str.contains(":")]
    if secs.empty:
        sel_cor = pr.PyRanges(sel_cor, int64=True)
    else:
        secs['transcript_id'] = secs['transcript_id'].str.split(":", expand=True)[1]
        sel_cor.loc[(sel_cor.index.isin(secs.index)), 'transcript_id'] = secs['transcript_id']
        sel_cor = pr.PyRanges(sel_cor, int64=True)

    # CALCULATING FRAMES

    print("\n*******************************")
    print("Calculating frames...")
    print("*******************************")

    ####Ensembl
    # Calculate frame for ensembl transcripts
    calculate_frame(ens_cor, by="transcript_id_ens")
    # Creating genome frame
    ens_cor = ens_cor.as_df()
    ens_cor.loc[ens_cor.Strand == "+", "Frame_genome"] = (
        ens_cor[ens_cor.Strand == "+"]["Start"]
        - ens_cor[ens_cor.Strand == "+"]["Frame"]
    ) % 3
    ens_cor.loc[ens_cor.Strand == "-", "Frame_genome"] = (
        ens_cor[ens_cor.Strand == "-"]["End"] + ens_cor[ens_cor.Strand == "-"]["Frame"]
    ) % 3
    ens_cor = pr.PyRanges(ens_cor, int64=True)

    ####Selenoprofiles
    # Dataframes only for CDS of selenoprofiles
    sel_CDS = sel_cor[sel_cor.Feature.str.startswith("CDS")]
    # Dataframe for storing Selenocysteines
    sel_Sec = sel_cor[sel_cor.Feature.str.startswith("Selenocysteine")]
    # Computing frames
    calculate_frame(sel_CDS, by='transcript_id')
    # Join Selenocysteines df with CDS+frame to assess Sec frame
    sel_frame = sel_Sec.join(sel_CDS, suffix="_CDS")
    sel_frame = sel_frame.as_df()
    # Now I have many useless columns, drop everything except the original ones df_sec and Frame (Filtering all the CDS columns)
    sel_frame = sel_frame.drop(sel_frame.filter(regex="_CDS").columns, axis=1)
    # Add output selenocysteines to sel_CDS
    sel_CDS = sel_CDS.as_df()
    ## Now Frame is the frame of CDS intervals for those for which Feature == 'CDS', and for Selenocysteine features, it is the frame of the CDS interval that contains them
    sel_cor = pd.concat([sel_CDS, sel_frame])
    # Creating genome frame column
    sel_cor.loc[sel_cor.Strand == "+", "Frame_genome"] = (
        sel_cor[sel_cor.Strand == "+"]["Start"]
        - sel_cor[sel_cor.Strand == "+"]["Frame"]
    ) % 3
    sel_cor.loc[sel_cor.Strand == "-", "Frame_genome"] = (
        sel_cor[sel_cor.Strand == "-"]["End"] + sel_cor[sel_cor.Strand == "-"]["Frame"]
    ) % 3
    sel_cor = pr.PyRanges(sel_cor, int64=True)

    # Join the ensembl transcript + selenocysteine position in a new pyranges
    annotation = sel_cor.join(ens_cor, how="left", suffix="_ens", report_overlap=True)

    # Useful for out_of_frame cases
    cds = annotation[annotation.Feature == "CDS"]
    sec = annotation[annotation.Feature == "Selenocysteine"]

    # ADD FILTER ON FRAME
    annotation = annotation.as_df()

    # When filtering TAKING INTO ACCOUNT THE CASES OF OUT OF FRAME. We added selenocysteine conditions because imagine: We have an ensembl tid which overlaps with selid but have different frame, this one will be wrongly annotated
    annotation_final = annotation.loc[
        (annotation["Frame_genome"] == annotation["Frame_genome_ens"])
        | (annotation["transcript_id_ens"] == "-1")
        | (annotation["Feature"] == "Selenocysteine")
    ]

    #####OUT OF FRAME CASES

    # Finding CDS which contain selenocysteines
    if sec.empty == False:
        merge = cds.join(sec, suffix="_sel")
        merge = merge.as_df()

        # Dropping columns
        result_clean = merge.drop(merge.filter(regex="_sel").columns, axis=1)

        # Filtering for out of frame CDS
        exons = result_clean[
            result_clean["Frame_genome"] != result_clean["Frame_genome_ens"]
        ]

        # Saving those cases
        exons = exons[exons["Start_ens"] != -1]
        out_of_frame = exons.drop_duplicates()

        # Convert annotation into dataframe to manipulate it correctly
        annotation_ok = pd.concat([annotation_final, out_of_frame])
    else:
        annotation_ok = annotation_final
    annotation_final = annotation_ok.reset_index(drop=True)

    # ASSESSING ANNOTATIONS

    print("\n*******************************")
    print("Assessing annotations...")
    print("*******************************")

    # Creating a dataframe with Selenocysteines
    df_sec = annotation_final[annotation_final["Feature"] == "Selenocysteine"]
    df_sec.reset_index(drop=True)

    # Creating a new column (Type of annotation) to annotate the different types of transcripts
    annotation_final["Type_annotation"] = "0"

    # Apply groupby + if for good annotations
    final_ann_df = annotation_final.groupby(['transcript_id', "transcript_id_ens"]).apply(
        lambda x: ass_ann(x,df_sec)
    )

    print("\n*******************************")
    print("Creating output tables...")
    print("*******************************")

    min_df = final_ann_df[
        ['transcript_id', "transcript_id_ens", "Type_annotation", "Feature"]
    ].drop_duplicates()

    # Check hierarchy

    def fn(df):
        if not (df.Type_annotation == "Missing").all():
            return df[
                (df.Type_annotation != "Missing") & (df.Feature != "Selenocysteine")
            ]
        else:
            if (df.Feature == "Selenocysteine").any():
                return df[df.Feature != "Selenocysteine"]
            else:
                return df

    min_df2 = min_df.groupby('transcript_id', as_index=False).apply(lambda x: fn(x))
    min_df2.reset_index(drop=True, inplace=True)
    min_df2.drop("Feature", axis=1, inplace=True)

    # Type annotation + hierarchy value
    agg_df = min_df.groupby('transcript_id', as_index=False).apply(lambda x: fn(x))
    agg_df["Type_annotation"] = agg_df["Type_annotation"].replace(
        ["Upstream", "Downstream", "Stop codon", "Out of frame", "Skipped", "Spliced"],
        "Missannotation",
    )

    # Sorting by hierarchy
    hierarchy = ["Well annotated", "Missannotation", "Missing"]
    mapping_dict = {value: index for index, value in enumerate(hierarchy)}
    agg_df["hierarchy"] = agg_df.Type_annotation.map(mapping_dict)
    agg_df.sort_values(by="hierarchy", inplace=True)
    agg = agg_df.groupby('transcript_id').first()
    agg = agg[["transcript_id_ens", "Type_annotation"]]

    # Saving to an output file
    min_df2.to_csv(opt["o"], sep="\t", index=False)
    agg.to_csv(opt["agg"], sep="\t")


if __name__ == "__main__":
    main()
