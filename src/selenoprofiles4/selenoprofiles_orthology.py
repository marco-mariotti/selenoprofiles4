#!/usr/bin/env python
import pandas as pd
from pyaln import Alignment
from easyterm import command_line_options
import os, subprocess, shlex
from .MMlib3 import *

help_msg = """selenoprofiles orthology: utility to classify selenoprofiles results according to orthologous groups.

Usage:   selenoprofiles orthology -i  fam1.ali [fam2.ali ... famN.ali]  [other options]

This utility takes the .ali files produced by selenoprofiles join as input.
As output, it produces one tsv file per input file, classifying the subfamily of each selenoprofiles prediction.

Output columns:
 - Candidate   sequence ID in selenoprofiles format (family.numericID.label.species.target_name)
 - Subfamily   subfamily assigned to this sequence
 - Similarity  similarity score between this sequence and the built-in "anchor" subfamily sequences
The output of selenoprofiles orthology can be used as input of selenoprofiles evolution (run: selenoprofiles evolution -h)

##### Compulsory input:
-i   selenoprotein family .ali file(s) produced by selenoprofiles join (run: selenoprofiles join -h)

### Output:
-of  output folder. Default: current directory
-o   suffix of the output file. Default: .orthology.

## Score similarity parameters:
For meaning of values, see https://pyaln.readthedocs.io/en/latest/alignment.html#pyaln.Alignment.score_similarity

-g   how to take into account gaps when comparing sequences. Possible values: {y,n,t,a} Default: n
-m   which similarity score metrics is used. Possible values: {i, w} Default: w
-w   if AWSI is computed (-m w), define weights per alignment column. Possible values: {m, i, q} Default: m

## Anchor input:
-af  override built-in path where anchor alignments are searched and loaded
"""

def_opt = {
    "i": [],
    "o": ".orthology.",
    "of": "./",
    "g": "n",
    "m": "w",
    "w": "m",
    "ap": "",
    "cmd": "orthology",
}


# Function to see if selenoprofiles output
def is_selenoprofiles_output_title(title):
    return (
        "chromosome:" in title
        and "target:" in title
        and "positions:" in title
        and "strand:" in title
    )


def candidate_score_similarity(ali, family, opt):
    anchor_path = opt["ap"]
    # Focusing on multimember families (Checking if sp family has anchor == Means it is multimember)
    anchor_ali_filename = anchor_path + family + ".anchor_ali.fa"
    if os.path.isfile(anchor_ali_filename):
        #### ALIGNMENT GENERATION from candidate.fa + anchor.fa
        ref_ali = alignment(anchor_ali_filename)

        # We need to convert Pyaln alignment to old class alignment() to use transfer_alignment function
        candidates_2 = alignment()
        for name, seq in ali:
            candidates_2.add(name + " " + ali.get_desc(name), seq)

        service(" ... transferring alignment ... ")
        # Run it
        joined_ali = ref_ali.transfer_alignment(candidates_2, dont_shrink=True)
        # Removing gap columns
        joined_ali.shrink()
        # COnverting transfered alignment to Pyaln Alignment
        cand_ali = Alignment(
            [(name, joined_ali.seq_of(name)) for name in joined_ali.titles()]
        )

        service(" ... computing similarity ... ")

        #### COMPUTE SCORE SIMILARITY
        # Getting candidate names
        names = cand_ali.names()
        # Filtering by only candidates
        names_filtered = [
            item
            for item in names
            if is_selenoprofiles_output_title(cand_ali.get_desc(item))
            and "SEED" not in item
        ]
        # We select only the anchor sequences to produce subfamily alignments
        anchor_names = [
            name.split(" ")[0] for name in ref_ali.titles() if "SEED" in name
        ]
        # First we're going to create a dictionary storing subfamily:list of sequences of that subfamily (GPX1:GPX1.SEED...)
        subfam_dict = {}
        # Then we're going to save each of the subfamilies (by iterating through the anchor names)
        for anc in anchor_names:
            # Selecting ANCHOR SEQUENCES only
            if not "SEED" in anc:
                continue
            # We get the subfamily of that specific sequence
            subfam = anc.split(".")[1]
            # Adding subfamily to the dict in case it isn't there
            if not subfam in subfam_dict:
                subfam_dict[subfam] = []
            # In case it is, we add the anchor sequence to the subfamily list
            subfam_dict[subfam].append(anc)

        df_dict = {"Candidate": [], "Subfamily": [], "Similarity": []}

        ##### pseudocode:
        
        # ## precompute candidate subalignment   ...
        # cand_only_ali=cand_ali[ candNames, :]

        # sim_df_list=[]
        # for subfam in subfamilies:
        #     #precompute subfam alig
        #     sbf_ali = cand_ali[ subfam_names, :]

        #     sbf_sim_df=sbf_ali.score_similarity(
        #         targets=cand_only_ali,
        #         gaps... )
        #     # rename  sbf_sim_df first and only col to subfam name
        #     sim_df_list.append( sbf_sim_df )

        # full_sim= pd.concat(sim_df_list, axis=1)
        # # index: cand seq names
        # # cols: one per subfam, each containing the score similarity


        
            
            
        
        for cand in names_filtered:
            # Automatic build and process of subalignments
            for i in list(subfam_dict.keys()):
                # Instead of iterating, we get the sequences from specific subfamily
                names_sbf = subfam_dict.get(i, [])
                names_sbf_copy = names_sbf.copy()
                # Add to the end of the alignment the candidate sequence
                names_sbf_copy.append(cand)

                # Generate the alignment for subfamily
                Sbf_alignment = cand_ali[names_sbf_copy, :].trim_gaps(1)

                # Compute the scoring similarity matrix
                a = Sbf_alignment[:-1, :].score_similarity(
                    targets=Sbf_alignment[[cand], :],
                    gaps=opt["g"],
                    metrics=opt["m"],
                    weights=opt["w"],
                )

                # Rename the index according to subfamily
                a = a.rename(index={(cand): (i)})

                # Save it to a dict
                df_dict["Candidate"].append(cand)
                df_dict["Subfamily"].append(i)
                metrics_name = a.columns[0]
                df_dict["Similarity"].append(a[metrics_name].values[0])

    else:
        profile_ali_filename = anchor_path + family + ".fa"
        # Calling the anchor alignment & candidate alignment
        ref_ali = alignment(profile_ali_filename)

        candidates_2 = alignment()
        for name, seq in ali:
            candidates_2.add(name + " " + ali.get_desc(name), seq)

        # Run it
        joined_ali = ref_ali.transfer_alignment(candidates_2, dont_shrink=True)
        # Removing gap columns
        joined_ali.shrink()
        # Converting alignment to Alignment from Pyaln
        cand_ali = Alignment(
            [(name, joined_ali.seq_of(name)) for name in joined_ali.titles()]
        )

        names = cand_ali.names()
        # We want to get rid of profile sequences
        names_filtered = [
            item
            for item in names
            if is_selenoprofiles_output_title(cand_ali.get_desc(item))
        ]

        # We select only the anchor sequences to produce subfamily alignments
        anchor_names = []
        for name in ref_ali.titles():
            anchor_names.append(name.split(" ")[0])

        df_dict = {"Candidate": [], "Subfamily": [], "Similarity": []}
        # Compute score similarity
        for cand in names_filtered:
            anchor_names.append(cand)

            # Creating alignment of anchor+candidate
            alig_c = cand_ali[anchor_names, :]
            a = alig_c[:-1, :].score_similarity(
                targets=alig_c[[cand], :],
                gaps=[opt["g"]],
                metrics=opt["m"],
                weights=opt["w"],
            )
            # Rename the index according to subfamily
            a = a.rename(index={(cand): (family)})
            # Save it to a dict
            df_dict["Candidate"].append(cand)
            df_dict["Subfamily"].append(family)
            metrics_name = a.columns[0]
            df_dict["Similarity"].append(a[metrics_name].values[0])

    # Save to output file the scoring similarity matrix
    score_subaln = pd.DataFrame.from_dict(df_dict)
    # Get the maximum score for each candidate
    max_idx = score_subaln.loc[score_subaln.groupby("Candidate")["Similarity"].idxmax()]
    max_idx.columns = ["Candidate", "Subfamily", "Similarity"]

    # TEMPORAL selection of candidates
    df_cand = max_idx[max_idx["Candidate"].isin(ali.names())]

    # #Adding species name (Gives a warning -- Talk with Marco)
    # df_cand["Species"] = df_cand["Candidate"].apply(lambda x: x.split(".")[3])

    return df_cand


def main(args={}):
    # if not args:
    #     opt = command_line_options(def_opt, help_msg)
    # else:
    #     opt = args
    opt = args

    if not opt["i"]:
        print(help_msg)
        sys.exit(1)

    for ffile in opt["i"]:
        file_name = os.path.basename(ffile)

        # Get the family : It needs the file name to be like GPx.ali, DI.ali, SelO.ali...
        fam = file_name.split(".")[0]

        # Read the alignment
        align = Alignment(ffile, fileformat="fasta")
        n_candidates = len(
            [None for title in align.titles() if is_selenoprofiles_output_title(title)]
        )
        write(
            f"Loaded alignment {ffile} : {n_candidates} predictions and {align.n_seqs()-n_candidates} profile seqs",
            1,
        )

        # Execute the function
        out = candidate_score_similarity(align, fam, opt)

        # Save output
        if not os.path.isdir(opt["of"]):
            os.mkdir(opt["of"])
        outfile = opt["of"].rstrip("/") + "/" + fam + "." + opt["o"].strip(".") + ".tsv"
        write(f"--> writing output: {outfile}", 1)
        out.to_csv(outfile, sep="\t", index=False)


# if __name__ == "__main__":
#     main()
