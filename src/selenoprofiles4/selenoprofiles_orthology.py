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

        cand_only_ali = cand_ali[names_filtered, :]

        # Initialize a list where we will save the results
        df_list = []

        # Automatic build and process of subalignments
        for i in list(subfam_dict.keys()):
            # Instead of iterating, we get the sequences from specific subfamily
            names_sbf = subfam_dict.get(i, [])
            # Subfamily alignment
            sbf_ali = cand_ali[names_sbf, :]
            # Compute score similarity
            sbf_sim_df = sbf_ali.score_similarity(
                targets=cand_only_ali, gaps=opt["g"], metrics=opt["m"], weights=opt["w"]
            )
            # Saving the score as subfamily name
            sbf_sim_df.columns = [i]
            # Df save to later concat
            df_list.append(sbf_sim_df)

    else:
        names = ali.names()
        # We want to get rid of profile sequences
        names_filtered = [
            item for item in names if is_selenoprofiles_output_title(ali.get_desc(item))
        ]

        # We select only the anchor sequences to produce subfamily alignments
        anchor_names = [
            item
            for item in names
            if not is_selenoprofiles_output_title(ali.get_desc(item))
        ]

        # Saving candidate alignment only
        cand_only_ali = ali[names_filtered, :]
        # Saving profile alignment
        prof_ali = ali[anchor_names, :]

        df_list = []
        # Computing score similarity
        sbf_ali = prof_ali.score_similarity(
            targets=cand_only_ali, gaps=opt["g"], metrics=opt["m"], weights=opt["w"]
        )
        # Renaming as selenoprotein family
        sbf_ali.columns = [family]
        # Save df to concat
        df_list.append(sbf_ali)

    # Save to output file the scoring similarity matrix
    result_df = pd.concat(df_list, axis=1)

    # Get the value and the column name with the highest score
    max_col = result_df.idxmax(axis=1)
    max_val = result_df.max(axis=1)

    # Create an output dataframe containing the candidate+subfamily+similarity
    df_cand = pd.DataFrame(
        {"Candidate": result_df.index, "Similarity": max_val, "Subfamily": max_col}
    )
    df_cand.reset_index(drop=True, inplace=True)

    df_cand["Species"] = df_cand["Candidate"].str.split(".").str[3]
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
