#!/usr/bin/env python
import pandas as pd
from pyaln import Alignment, pyaln_folder
from easyterm import command_line_options
import os, subprocess, shlex
from ncbi_db import ncbi_taxonomy_tree
from MMlib3 import *

help_msg = """selenoprofiles evolution: utility to exclude non-expected genes predicted by Selenoprofiles.

Usage: selenoprofiles evolution -i fam1.orthology.tsv [fam2.orthology.tsv ... famN.orthology.tsv]  [other options]

This utility takes the .tsv files produced by selenoprofiles orthology as input.
As output, it produces one tsv file per input file, filtering non-expected predictions.

Output columns:
 - Candidate            sequence ID in selenoprofiles format (family.numericID.label.species.target_name)
 - Subfamily            subfamily assigned to this sequence
 - Similarity           similarity score between this sequence and the built-in "anchor" subfamily sequences
 - Species              species name of the predicted sequence
 - Pass_filter          boolean column indicating if the sequence is expected or not
 - Discard_description  column that provides an explanation for why the filter criteria were not met

#### Compulsory input
-i   selenoprotein family fam.orthology.tsv file(s) produced by selenoprofiles orthology (run: selenoprofiles orthology -h)

### Output
-of    output folder. Default: current directory.
-o     suffix of the output file. Default: .evolution.
-temp  temporal folder to save intermediate files.
-a     optional output .ali file. Requires the input .ali file used in selenoprofiles orthology. Outputs an alignment of the
       filtered sequences.

### Optional parameters
-switch  decides whether to keep or not selenoprotein homologs
-exp     provide own expectation table
-map     map manually species to lineage. Avoids using NCBI_DB. User needs to provide a species /t lineage table.
-pexp    print the expectation table

"""

def_opt = {
    "i": [],
    "o": ".evolution.",
    "of": "./",
    "temp": "temp",
    "exp": "/Lobster/mtico/expectation_table/Expectation_table_2.csv",
    "switch": 0,
    "a": "",
    "map": 0,
    "pexp": 0,
}


def is_selenoprofiles_output_title(title):
    return (
        "chromosome:" in title
        and "target:" in title
        and "positions:" in title
        and "strand:" in title
    )


# Function used to obtain the lineage for each specie
def assign_lineage(species, taxonomy_lineages, mapping_df, usr):
    if species not in taxonomy_lineages.keys():
        # Avoid repeating a lot of times the same question for the same species
        if species in usr:
            return 0
        else:
            # Check for the species which go from Anas_pl to Anas_pl_pl
            if len(species.split("_")) == 2:
                # Code to search for dict keys which contain that species
                matching_key = [
                    key for key in taxonomy_lineages.keys() if species in key
                ]
                # In case the species is not present
                try:
                    lineage_list = [
                        item.strip(";")
                        for item in taxonomy_lineages[matching_key[0]].split(" ")
                    ]
                    for lineage in lineage_list:
                        if lineage in mapping_df["Lineage"].values:
                            return lineage

                except KeyError:
                    # raise Exception(f"'{matching_key[0]}' not found in the taxonomy lineage")
                    user_choice = input(
                        f"'{matching_key[0]}' not found in the taxonomy lineage. Do you want to continue (Y/N)? "
                    )
                    if user_choice.lower() == "y":
                        # Assign lineage as 0
                        usr.append(species)
                        return 0
                    else:
                        print("Analysis cancelled by the User")
                        exit()

            # Code for those species which are sus_scrofa_hampshire to search only for sus_scrofa
            else:
                sp_words = species.split("_")
                shortened_sp = "_".join(sp_words[:-1])
                # In case shortened_sp does not exist in the lineage dictionary
                try:
                    lineage_list = [
                        item.strip(";")
                        for item in taxonomy_lineages[shortened_sp].split(" ")
                    ]
                    for lineage in lineage_list:
                        if lineage in mapping_df["Lineage"].values:
                            return lineage
                except KeyError:
                    # raise Exception(f"'{shortened_sp}' not found in the taxonomy lineage")
                    user_choice = input(
                        f"'{shortened_sp}' not found in the taxonomy lineage. Do you want to continue (Y/N)? "
                    )
                    if user_choice.lower() == "y":
                        # Assign lineage as 0
                        usr.append(species)
                        return 0
                    else:
                        print("Analysis cancelled by the User")
                        exit()
    else:
        # If species is present
        lineage_list = [
            item.strip(";") for item in taxonomy_lineages[species].split(" ")
        ]
        for lineage in lineage_list:
            if lineage in mapping_df["Lineage"].values:
                return lineage


# This function is used for creating a description for discarded sequences
def create_discard_description(row):
    if not row["Pass_filter"]:
        if row["Subfamily"] == "0":
            # In case the predicted SUbfamilies were GPX5,GPX7,GPX8
            return "Wrong family classification"
        else:
            return "Low similarity score"
    else:
        return ""


def expectations(table, family, opt):
    if not os.path.exists(opt["temp"]):
        os.mkdir(opt["temp"])
    # Loading expectation table
    expectation_table = pd.read_csv(opt["exp"])

    # Print expectation table if user wants
    if opt["pexp"]:
        service(" ... printing expectation table ... ")
        print(expectation_table)
        exit()

    # Option to select or not selenocysteines
    if opt["switch"]:
        candidates = table
    else:
        # Selecting only selenocysteines
        candidates = table[table["Candidate"].str.contains("selenocysteine")]

    if family == "GPx":
        # Selenoprofiles orthology can score sequences to GPX8,GP7,GPX5. We should get rid of them
        #   candidates['Subfamily']=candidates['Subfamily'].replace(['GPX5','GPX7','GPX8'],'0')
        candidates.replace(["GPX1B", "GPX3B"], ["GPX1", "GPX3"], inplace=True)

    # First group by subfamily and species name to index best values
    sorted_df = candidates.sort_values(by=["Species", "Subfamily", "Similarity"])
    # Now assign index best by sort position
    sorted_df["Index_best"] = sorted_df.groupby(["Species", "Subfamily"]).cumcount(
        ascending=False
    )
    # Check if the user wants to manually map species to lineage
    if not opt["map"]:
        ## We need taxonomy lineages to be able to compare with expectation table. So we are going to use ncbi_taxonomy.
        # To NCBI format
        sorted_df["Species"] = sorted_df["Species"].str.capitalize()
        sorted_df["Species"].to_csv(
            opt["temp"] + "/species_list.csv", header=None, index=False
        )
        # Running ncbi_taxonomy
        args = ncbi_taxonomy_tree.parse_opts(
            "-n temp/species_list.csv -u --lineage -temp temp/ "
        )
        lineage_annotated_tree = ncbi_taxonomy_tree.main(args)

        # Getting the lineage for each node
        taxonomy_lineages = {}
        missing_species = {}
        for species in sorted_df["Species"]:
            try:
                unmasked_species = species.replace("_", " ")
                leaf = lineage_annotated_tree & unmasked_species
                taxonomy = "; ".join(leaf.lineage.split("||")[:-1])
                taxonomy_lineages[species] = taxonomy
            except:
                missing_species[species] = 1

        # We add a lineage column to the df in order to be able to compare with the expectations
        user_choices = []
        sorted_df["Lineage"] = sorted_df["Species"].apply(
            lambda x: assign_lineage(
                x, taxonomy_lineages, expectation_table, user_choices
            )
        )

    else:
        inp = sorted_df
        manual_tab = pd.read_csv(opt["map"])
        sorted_df = pd.merge(inp, manual_tab, on="Linage", how="outer")

    # Merging df to compare subfamilies
    joined = pd.merge(sorted_df, expectation_table, on="Lineage", how="outer")
    # Comparing Expected values to Index best similarity scores
    joined_sec = joined.astype({"Index_best": "int"})
    joined_sec["Pass_filter"] = joined_sec.apply(
        lambda row: row.get(row["Subfamily"], 0) > row["Index_best"], axis=1
    )
    joined_sec["Pass_filter"] |= ~joined_sec["Subfamily"].isin(joined_sec.columns)

    # Now we should do a table with the number of discarded sequences for each subfamily
    numbers = (
        joined_sec.groupby("Subfamily")["Pass_filter"]
        .value_counts()
        .unstack()
        .reset_index()
    )

    final_table = joined_sec[
        ["Candidate", "Subfamily", "Similarity", "Species", "Index_best", "Pass_filter"]
    ]

    final_table["Discard_description"] = final_table.apply(
        create_discard_description, axis=1
    )

    return final_table


def main(args={}):

    opt = args

    # Iterate through each family alignment
    for sp in opt["i"]:
        filename = os.path.basename(sp)
        # Get the selenoprotein family name (in case it is GPx.candidate_matrix.csv)
        fam = filename.split("_")[0]
        # Read the selenoprotein family input file
        candidates = pd.read_csv(sp, sep="\t")
        # Execute the filter
        out = expectations(candidates, fam, opt)

        # Save output
        if not os.path.isdir(opt["of"]):
            os.mkdir(opt["of"])
        outfile = opt["of"].rstrip("/") + "/" + fam + "." + opt["o"].strip(".") + ".tsv"
        write(f"--> writing output: {outfile}", 1)
        out.to_csv(outfile, sep="\t", index=False)

    if opt["a"]:
        ali = Alignment(opt["a"], fileformat="fasta")
        # List to save true selenocysteines + profile seqs
        filt_names = [
            item
            for item in ali.names()
            if not is_selenoprofiles_output_title(ali.get_desc(item))
        ]
        prof_cand = filt_names + out[out["Pass_filter"] == True]["Candidate"].tolist()
        # Creating alignment file
        ali_filt = ali[prof_cand, :]
        # Saving output
        outfile = opt["of"].rstrip("/") + "/" + fam + "." + opt["o"].strip(".") + ".ali"
        ali_filt.write(fileformat="fasta", to_file=outfile)


#if __name__ == "__main__":
 #   main()
