#!/usr/bin/env python
import pandas as pd
from pyaln import Alignment, pyaln_folder
from easyterm import command_line_options
import os, subprocess, shlex
from ncbi_db import ncbi_taxonomy_tree
import numpy as np
from .MMlib3 import *

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
    "exp": "",
    "switch": 0,
    "a": [],
    "map": 0,
    "pexp": 0,
    "cmd": "lineage",
    "ann": "",
    "l":False,
}


def is_selenoprofiles_output_title(title):
    return (
        "chromosome:" in title
        and "target:" in title
        and "positions:" in title
        and "strand:" in title
    )

def rename(self,d):
   self._ord = [d[n] if n in d else n for n in self._ord]
   self._seqs={d.get(n, n):s for n,s in self._seqs.items()}
   self._desc={d.get(n, n):des for n,des in self._desc.items()}

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
        else:
            raise Exception(f"'{species}' not found in the taxonomy lineage")

    else:
        # If species is present
        lineage_list = [
            item.strip(";") for item in taxonomy_lineages[species].split(" ")
        ]
        for lineage in reversed(lineage_list):
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
        joined = pd.merge(sorted_df, expectation_table, on="Lineage", how="inner")

    else:
        inp = sorted_df
        manual_tab = pd.read_csv(opt["map"])
        joined = pd.merge(inp, manual_tab, on="Species", how="inner")

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
    if family in ['GPx','DI','TXNRD']:
        fam_cols = [ col for col in joined_sec.columns.tolist() if family.upper() in col]
    else:
        fam_cols = [ col for col in joined_sec.columns.tolist() if family in col]
    missing = pd.DataFrame(columns=['Species', 'Subfamily', 'Count'])

    # Merging it to the previous df
    test = pd.merge(joined_sec, grouped_counts, on=["Species", "Subfamily"], how="left")
    # We have a problem, we are not taking into account those subfamilies which are not in that species but they should be
    seen_species = {}
    for index, row in test.iterrows():
      species = row['Species']
      subfamily = row['Subfamily']
      if species not in seen_species:
         seen_species[species] = []
         if subfamily not in seen_species[species]:
           seen_species[species].append(subfamily)
           temp_sp = test[(test['Species'] == species)]
           missing_fam_cols = [col for col in fam_cols if col not in temp_sp['Subfamily'][temp_sp['Species']==species].unique()]
           # If any value is missing, add a single row to the missing DataFrame
           for col in missing_fam_cols:
               count_value = row[col]
               for i in range(count_value):
                    missing = missing.append({'Species': species, 'Subfamily': col, 'Count': i}, ignore_index=True)

           # Checking if rows are missing in nonmultimember families
           temp = test[(test['Species'] == species) & (test['Subfamily'] == subfamily)]
           expected_count = temp[subfamily].unique()[0]
           actual = temp['count'].unique()[0]
           if int(expected_count) > actual:
              # Calculate the difference
              difference = int(expected_count) - actual
              for i in range(difference):
                missing = missing.append({'Species': row['Species'], 'Subfamily': row['Subfamily'], 'Count': i}, ignore_index=True)
      else:
         if subfamily not in seen_species[species]:
           seen_species[species].append(subfamily)
           # Checking if rows are missing in nonmultimember families
           temp = test[(test['Species'] == species) & (test['Subfamily'] == subfamily)]
           expected_count = temp[subfamily].unique()[0]
           actual = temp['count'].unique()[0]
           if int(expected_count) > actual:
              # Calculate the difference
              difference = int(expected_count) - actual
              for i in range(difference):
                missing = missing.append({'Species': row['Species'], 'Subfamily': row['Subfamily'], 'Count': i}, ignore_index=True)

    # Now we need to add the missings for those species absent in the current file
    absent = {species : d[species] for species in d if species not in test['Species'].unique()}
    # Creating artificial df to assign lineage
    abs_df = pd.DataFrame(list(absent.keys()), columns=['Species'])
    # Now we can assign lineage
    abs_df['Lineage'] = abs_df['Species'].apply(
            lambda x: assign_lineage(x, absent, expectation_table)
    )
    # Joining absents with expectation
    abs_join = pd.merge(abs_df, expectation_table, on="Lineage", how="inner")
    # Creating missing sequences
    # If any value is missing, add a single row to the missing DataFrame
    for _,row in abs_join.iterrows():
        species = row['Species']
        for col in fam_cols:
            count_value = row[col]
            for i in range(count_value):
                missing = missing.append({'Species': species, 
                'Subfamily': col, 
                'Count': i}, ignore_index=True
                )

    # Adding Pass_filter to false
    missing["Pass_filter"] = False

    # Adding to the end missing predictions
    joined_sec = pd.concat([joined_sec, missing], ignore_index=True, sort=False)

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

    os.makedirs("temp", exist_ok=True)
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

        # Creating sp drawwer input
        if (len(opt["a"]) > 0) and (len(opt["ann"]) > 0):
            annotation = pd.read_csv(opt['ann'],sep="\t",header=None)
            # Reading the alignment
            for aln in opt['a']:
              if fam in aln:
                ali = Alignment(aln, fileformat="fasta")

            # Producing ann + filt
            final_df = pd.DataFrame()
            # Producing tsv file from lineage output + annotations 
            missing_sp = out[out['Candidate'].isna()]
            out['Candidate'] = out['Candidate'].str.split('.').str[:4].str.join('.')
            merged_df = pd.merge(annotation, out, left_on=0, right_on='Candidate', how='inner')
            final_df = pd.concat([merged_df, missing_sp]) # DI_filtered.tsv / annotations

            # Creating names+seqs for missing predicitons
            names = []
            seqs = []
            occurrences = {}
            for index, row in final_df[final_df[0].isna()].iterrows():
                subfamily_letters = ''.join([c for c in row['Subfamily'] if c.isalpha()])
                subfamily_numbers = ''.join([c for c in row['Subfamily'] if c.isdigit()])
                species= row['Species']
                subfamily = row['Subfamily']
                species_lower = species[0].lower() + species[1:]
                hyper_str = "-"*(ali.ali_length()-1) + "A" # Sequence str
                #name = f"{subfamily_letters}.{subfamily_numbers}.Missing.{species}.unfiltered" # Name str
                key = (species, subfamily_letters)
                if key in occurrences:
                    occurrences[key] += 1
                    name = f"{subfamily}.{subfamily_numbers}.Missing.{species_lower}.{species}_target{occurrences[key]}.unfiltered chromosome:Missing strand:+ positions:1-50,60-110,120-170"
                else:
                    occurrences[key] = 1
                    name = f"{subfamily}.{subfamily_numbers}.Missing.{species_lower}.{species}_target{occurrences[key]}.unfiltered chromosome:Missing strand:+ positions:1-50,60-110,120-170"

                names.append(name)
                seqs.append(hyper_str)

            # Getting those names in alignment which are in the annotations  
            #matching = [name for name in ali.names() if any(str(subs) in str(name) for subs in annotations['0'] if pd.notna(subs))]
            f_dict = {}
            matching = []
            for name in ali.names():
               for subs in final_df[0]:
                   if pd.notna(subs) and (str(subs) in str(name)) and ((str(subs)+'_') not in str(name)):
                       matching.append(name)
                       if final_df.loc[final_df[0]==subs,'Pass_filter'].bool():
                           f_dict[name] = name+'.unfiltered'
                       else:
                           f_dict[name] = name+'.filtered'
                       break

            ali_f = ali[matching,:]

            # Adding missing sequences
            for i in range(len(seqs)):
                ali_f.add_seq(names[i],seqs[i])

            rename(ali_f,f_dict)
            # Changing selenocysteine by readthrough, arginine, unaligned
            seq_dict = {}
            replace_dict={
                row[0]: "arginine" if row[2] == "Well" else "readthrough" if row[2] == "Missannotation" else "unaligned"
                for index, row in final_df.iterrows()
            }

            # There are some families with empty series at the end
            last_key, last_value = list(replace_dict.items())[-1]
            if str(last_key) == 'nan':
               replace_dict.popitem()

            #Replacing each selenocysteine by readthrough, arginine, unaligned
            for name in ali_f.names():
                for key in replace_dict.keys():
                    if (str(key) in str(name)) and ((str(key) + "_") not in str(name)):
                        new_name=name.replace("selenocysteine",replace_dict[key])
                        if fam != "TXNRD":
                          new_name_2 = new_name.replace(fam,final_df['Subfamily'][final_df[0] == key].values[0])
                        else:
                          new_name_2 = new_name.replace("TR",final_df['Subfamily'][final_df[0] == key].values[0])
                        seq_dict[name]=new_name_2

            rename(ali_f,seq_dict)

            # Saving output
            outfile = opt["of"].rstrip("/") + "/" + fam + "." + opt["o"].strip(".") + ".drawer.ali"
            ali_f.write(fileformat="fasta", to_file=outfile)

        elif opt["a"]:
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
 #   main(	)
