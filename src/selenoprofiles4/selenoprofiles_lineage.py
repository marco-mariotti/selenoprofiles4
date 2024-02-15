#!/usr/bin/env python
import pandas as pd
from pyaln import Alignment, pyaln_folder
from easyterm import command_line_options
import os, subprocess, shlex
from ncbi_db import ncbi_taxonomy_tree
import numpy as np
from MMlib3 import *

help_msg = """selenoprofiles lineage: utility to exclude non-expected genes predicted by Selenoprofiles.

Usage: selenoprofiles lineage -i fam1.orthology.tsv [fam2.orthology.tsv ... famN.orthology.tsv]  [other options]

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
-o     suffix of the output file. Default: .lineage.
-temp  temporal folder to save intermediate files.
-a     optional output .ali file. Requires the input .ali file used in selenoprofiles orthology. Outputs an alignment of the
       filtered sequences.

### Optional parameters
-all     decides whether to keep or not selenoprotein homologs
-exp     provide own expectation table
-map     map manually species to lineage. Avoids using NCBI_DB. User needs to provide a species /t lineage table.
-pexp    print the expectation table
-l       include lineage in the output table

"""

def_opt = {
    "i": [],
    "o": ".lineage.",
    "of": "./",
    "temp": "temp",
    "exp": "/Lobster/mtico/expectation_table/Expectation_table_2.csv",
    "switch": 0,
    "a": "",
    "map": 0,
    "pexp": 0,
    "l": 0,
}


def is_selenoprofiles_output_title(title):
    return (
        "chromosome:" in title
        and "target:" in title
        and "positions:" in title
        and "strand:" in title
    )


def species2lineage(df, d, opt):
    unique_species = df["Species"].unique()
    if not all(species in d for species in unique_species):
        ## We need taxonomy lineages to be able to compare with expectation table. So we are going to use ncbi_taxonomy.
        # To NCBI format
        df["Species"] = df["Species"].str.capitalize()
        df["Species"].to_csv(
            opt["temp"] + "/species_list.csv", header=None, index=False
        )
        args = ncbi_taxonomy_tree.parse_opts(
            "-n temp/species_list.csv -u --lineage -temp temp/ "
        )
        lineage_annotated_tree = ncbi_taxonomy_tree.main(args)

        # Getting the lineage for each node
        missing = {}
        for species in df["Species"]:
            try:
                unmasked_species = species.replace("_", " ")
                leaf = lineage_annotated_tree & unmasked_species
                taxonomy = "; ".join(leaf.lineage.split("||")[:-1])
                d[species] = taxonomy
            except:
                missing[species] = 1
        return d
    else:
        return d


# Function used to obtain the lineage for each specie
def assign_lineage(species, taxonomy_lineages, mapping_df):
    if species not in taxonomy_lineages.keys():
        # Code for those species which are sus_scrofa_hampshire to search only for sus_scrofa
        if len(species.split("_")) == 3:
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
                raise Exception(f"'{shortened_sp}' not found in the taxonomy lineage")
                return 0
        else:
            raise Exception(f"'{species}' not found in the taxonomy lineage")
            return 0

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
        if pd.isna(row["Candidate"]):
            # In case the predicted SUbfamilies were GPX5,GPX7,GPX8
            return "Missing prediction"
        else:
            return "Low similarity score"
    else:
        return ""


def expectations(table, family, opt, d):
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
    if opt["all"]:
        candidates = table
    else:
        # Selecting only selenocysteines
        candidates = table[table["Candidate"].str.contains("selenocysteine")]

    if family == "GPx":
        # candidates['Subfamily']=candidates['Subfamily'].replace(['GPX5','GPX7','GPX8'],'0')
        candidates.replace(["GPX1B", "GPX3B"], ["GPX1", "GPX3"], inplace=True)

    elif family == "SelW":
        candidates.replace(["SelW1", "SelW2"], ["SelW","SelW"], inplace=True)

    # First group by subfamily and species name to index best values
    sorted_df = candidates.sort_values(by=["Species", "Subfamily", "Similarity"])
    # Now assign index best by sort position
    sorted_df["Index_best"] = sorted_df.groupby(["Species", "Subfamily"]).cumcount(
        ascending=False
    )
    # Check if the user wants to manually map species to lineage
    if not opt["map"]:
        # We add a lineage column to the df in order to be able to compare with the expectations
        sorted_df["Lineage"] = sorted_df["Species"].apply(
            lambda x: assign_lineage(x, d, expectation_table)
        )

    else:
        inp = sorted_df
        manual_tab = pd.read_csv(opt["map"])
        sorted_df = pd.merge(inp, manual_tab, on="Linage", how="inner")

    # Merging df to compare subfamilies
    joined = pd.merge(sorted_df, expectation_table, on="Lineage", how="inner")
    # Comparing Expected values to Index best similarity scores
    joined_sec = joined.astype({"Index_best": "int"})
    joined_sec["Pass_filter"] = joined_sec.apply(
        lambda row: row.get(row["Subfamily"], 0) > row["Index_best"], axis=1
    )

    ## Counting missing predictions
    # Counting the number of rows per species + subfamily
    grouped_counts = (
        joined_sec.groupby(["Species", "Subfamily"]).size().reset_index(name="count")
    )
    # Merging it to the previous df
    test = pd.merge(joined_sec, grouped_counts, on=["Species", "Subfamily"], how="left")
    # Determining if rows are missing or not
    test["Missings"] = test.apply(
        lambda row: row.get(row["Subfamily"], 0) > row["count"], axis=1
    )
    # Getting the rows which are Missing
    new_rows = test[test["Missings"]].copy()
    # Just maintaining Species & Subfamily : other values NA
    new_rows.loc[:, new_rows.columns.difference(["Species", "Subfamily"])] = np.nan
    # Adding Pass_filter to false
    new_rows["Pass_filter"] = False

    # Adding to the end missing predictions
    joined_sec = pd.concat([joined_sec, new_rows], ignore_index=True, sort=False)

    joined_sec["Pass_filter"] |= ~joined_sec["Subfamily"].isin(joined_sec.columns)

    # Now we should do a table with the number of discarded sequences for each subfamily
    numbers = (
        joined_sec.groupby("Subfamily")["Pass_filter"]
        .value_counts()
        .unstack()
        .reset_index()
    )

    if opt["l"]:
        final_table = joined_sec[
            [
                "Candidate",
                "Subfamily",
                "Similarity",
                "Species",
                "Lineage",
                "Pass_filter",
            ]
        ]
    else:
        final_table = joined_sec[
            ["Candidate", "Subfamily", "Similarity", "Species", "Pass_filter"]
        ]

    # We are going to add description to our table
    final_table["Discard_description"] = final_table.apply(
        create_discard_description, axis=1
    )

    return final_table


def main(args={}):
    opt = args
    sp2lin = {}

    # Iterate through each family alignment
    for sp in opt["i"]:
        filename = os.path.basename(sp)
        # Get the selenoprotein family name (in case it is GPx.candidate_matrix.csv)
        fam = filename.split(".")[0]
        # Read the selenoprotein family input file
        candidates = pd.read_csv(sp, sep="\t")
        # Creating here species to lineage
        species2lineage(candidates, sp2lin, opt)
        # Execute the filter
        out = expectations(candidates, fam, opt, sp2lin)

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


# if __name__ == "__main__":
#    main(	)
