#!/usr/bin/env python
import pandas as pd
from pyaln import Alignment, pyaln_folder
from easyterm import command_line_options
import os
from ncbi_db import ncbi_taxonomy_tree

help_msg='''selenoprofiles evolution: utility to exclude non-expected genes predicted by Selenoprofiles.

### Input/Output:
-i Selenoprotein family .tsv files produced by Selenoprofiles orthology feature.
-o Specify the suffix of the output file.
-of Specify output folder

### Script description:
This script has been designed to take all .tsv files produced by Selenoprofiles orthology feature as input. The
output produced is a tsv file containing the input columns from Selenoprofiles orthology (Candidate sequence, Subfamily, 
AWSI (score similarity), Species) + three new columns: Index_best, Discarded & Discard_description.

'''

def_opt = {
     "i":[],
     "o":"_expectations",
     "of":"",
}

# This function is used for creating a description for discarded sequences
def create_discard_description(row):
    if not row['Expected']:
        if row['Subfamily'] == '0':
            # In case the predicted SUbfamilies were GPX5,GPX7,GPX8
            return 'Wrong family classification'
        else:
            return 'Low similarity score'
    else:
        return ''

def expectations(table, family):
  # Loading expectation table
  expectation_table=pd.read_csv("/Lobster/mtico/expectation_table/Expectation_table.csv")

  # Selecting only selenocysteines
  candidates = table[table['Candidate'].str.contains('selenocysteine')]
  if family == 'GPx':
      #Selenoprofiles orthology can score sequences to GPX8,GP7,GPX5. We should get rid of them
      candidates['Subfamily']=candidates['Subfamily'].replace(['GPX5','GPX7','GPX8'],'0')
      #Look at Subfamily column + column with the value from this one and compare Index_best with value from column taken before
      candidates.replace(['GPX1B','GPX3B'],['GPX1','GPX3'],inplace=True)

  # First group by subfamily and species name to index best values
  sorted_df=candidates.sort_values(by=['Species','Subfamily','AWSI'])
  # Now assign index best by sort position
  sorted_df['Index_best']=sorted_df.groupby(['Species','Subfamily']).cumcount(ascending=False)
  ## We need taxonomy lineages to be able to compare with expectation table. So we are going to use ncbi_taxonomy.
  #To NCBI format
  sorted_df['Species']=sorted_df['Species'].str.capitalize()
  sorted_df['Species'].to_csv("species_list.csv",header=None, index=False)
  #Running ncbi_taxonomy
  args = ncbi_taxonomy_tree.parse_opts( '-n species_list.csv -u --lineage -temp temp/ ')
  lineage_annotated_tree =  ncbi_taxonomy_tree.main(args)

  #Getting the lineage for each node
  taxonomy_lineages={}
  missing_species={}
  for species in sorted_df['Species']:
   try:
     unmasked_species=species.replace('_',' ')
     leaf = lineage_annotated_tree & unmasked_species
     taxonomy = '; '.join(leaf.lineage.split('||')[:-1])
     taxonomy_lineages[species]=taxonomy
   except:
    missing_species[species]=1

  #We add a lineage column to the df in order to be able to compare with the expectations
  lineages_dict={'Sauropsida':'Birds','Actinopterygii':'Fish','Anura':'Frog','Eutheria':'Placentals','Metatheria':'Marsupials','Prototheria':'Platypus'}
  for species in taxonomy_lineages:
    for key in lineages_dict:
      if key in taxonomy_lineages[species]:
        sorted_df.loc[sorted_df['Species']==species,'Lineage']=lineages_dict[key]

  #Merging df to compare subfamilies
  joined=pd.merge(sorted_df,expectation_table,on="Lineage")

  #Comparing Expected values to Index best similarity scores
  joined_sec = joined.astype({'Index_best':'int'})
  joined_sec['Expected'] = joined_sec.apply(lambda row: row.get(row['Subfamily'], 0) > row['Index_best'], axis=1)

  #Now we should do a table with the number of discarded sequences for each subfamily
  numbers=joined_sec.groupby('Subfamily')['Expected'].value_counts().unstack().reset_index()

  final_table = joined_sec[['Candidate','Subfamily','AWSI','Species','Index_best','Expected']]

  final_table['Discard_description'] = final_table.apply(create_discard_description, axis=1)

  return final_table


def main(args={}):

  if not args:
    opt = command_line_options(def_opt, help_msg)
  else:
    opt = args

  # Iterate through each family alignment
  for sp in opt['i']:
   filename = os.path.basename(sp)
   # Get the selenoprotein family name (in case it is GPx.candidate_matrix.csv)
   fam = filename.split("_")[0]
   # Read the selenoprotein family input file
   candidates = pd.read_csv(sp,sep="\t")
   # Execute the filter
   expected = expectations(candidates,fam)
   # Save it to a file
   expected.to_csv(opt['of']+fam+opt['o']+'.tsv', sep="\t", index=False)

if __name__ == "__main__":
    main()
