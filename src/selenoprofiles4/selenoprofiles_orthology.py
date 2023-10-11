#!/usr/bin/env python
import pandas as pd
from pyaln import Alignment, pyaln_folder
from easyterm import command_line_options
import os
import subprocess
import shlex
from MMlib3 import *

help_msg="""selenoprofiles orthology: utility to classify selenoprofiles results according to orthologous groups.

### Input/Output:
-i Selenoprotein family .ali files produced by Selenoprofiles.
-o Specify the suffix of the output file.
-of Specify output folder
-g Specify how to take into account gaps when comparing sequences with anchor alignments.
   Possible values:
   - ‘y’ : gaps are considered and considered mismatches. Positions that are gaps in both sequences are ignored.
   - ‘n’ : gaps are not considered. Positions that are gaps in either sequences compared are ignored.
   - ‘t’ : terminal gaps are trimmed. Terminal gap positions in either sequences are ignored, others are considered as in ‘y’.
   - ‘a’ : gaps are considered as any other character; even gap-to-gap matches are scored as identities.
   Default is 'n'.

-m Specify which metrics are computed.
   Possible values:
   - ‘i’ : average sequence identity, aka ASI.
   - ‘w’ : weighted sequence identity, aka AWSI.
   Multiple arguments may be concatenated (e.g. ‘iw’) to compute all of the possibilities.
   Default is 'w'.

-w If AWSI is computed, defines how weights are computed for each alignment column.
   Possible values:
   - ‘m’ : maximum frequency of non-gap character in self.
   - ‘i’ : information content, i.e. 2- sum(p*log2(p)) where p is frequency of non-gap characters.
   - ‘q’ : quadratic sum, i.e. sum(p*p) where p is frequency of non-gap characters.
   Multiple arguments may be concatenated (e.g. ‘mi’) to compute all of the possibilities.
   Default is 'm'.

### Script description:
This script has been designed to take the .ali files produced by Selenoprofiles join as input. The output produced is a tsv file containing
four main columns: Candidate sequence, Subfamily, AWSI (score similarity), Species. This output can be used for Selenoprofiles evolution feature.

"""

def_opt = {
    "i":[],
    "o":"_Classified",
    "of":"",
    "g":"n",
    "m":"w",
    "w":"m"}

# Function to see if selenoprofiles output
def is_selenoprofiles_output_title(title):
    return (
        "chromosome:" in title
        and "target:" in title
        and "positions:" in title
        and "strand:" in title
    )

def candidate_score_similarity(ali,family,opt):
  anchor_path = '/Lobster/mtico/selenoprofiles_orthology/anchor_alignments/'
  # Focusing on multimember families (Checking if sp family has anchor == Means it is multimember)
  if os.path.isfile(anchor_path+family+"_anchor.ali"):
    #### ALIGNMENT GENERATION from candidate.fa + anchor.fa
    ref_ali = alignment(anchor_path+family+"_anchor.ali")

    # We need to convert Pyaln alignment to old class alignment() to use transfer_alignment function
    candidates_2 = alignment()
    for name, seq in ali:
       candidates_2.add(name+' '+ali.get_desc(name), seq)

    #Run it
    joined_ali = ref_ali.transfer_alignment(candidates_2,dont_shrink = True)
    # Removing gap columns
    joined_ali.shrink()
    # COnverting transfered alignment to Pyaln Alignment
    cand_ali=Alignment( [ (name, joined_ali.seq_of(name)) for name in joined_ali.titles()] )

    #### COMPUTE SCORE SIMILARITY
    # Getting candidate names
    names=cand_ali.names()
    # Filtering by only candidates
    names_filtered = [item for item in names if is_selenoprofiles_output_title(cand_ali.get_desc(item)) and 'SEED' not in item]
    #We select only the anchor sequences to produce subfamily alignments
    anchor_names = [name.split(' ')[0] for name in ref_ali.titles() if 'SEED' in name]
    # First we're going to create a dictionary storing subfamily:list of sequences of that subfamily (GPX1:GPX1.SEED...)
    subfam_dict = {}
    # Then we're going to save each of the subfamilies (by iterating through the anchor names)
    for anc in anchor_names:
      #Selecting ANCHOR SEQUENCES only
      if not 'SEED' in anc:
        continue
      #We get the subfamily of that specific sequence
      subfam = anc.split('.')[1]
      #Adding subfamily to the dict in case it isn't there
      if not subfam in subfam_dict:
        subfam_dict[subfam]=[]
      #In case it is, we add the anchor sequence to the subfamily list
      subfam_dict[subfam].append(anc)

    df_dict={'Candidate':[],'Subfamily':[],'AWSI':[]}

    for cand in names_filtered:
      #Automatic build and process of subalignments
        for i in list(subfam_dict.keys()):
          #Instead of iterating, we get the sequences from specific subfamily
          names_sbf = subfam_dict.get(i, [])
          names_sbf_copy = names_sbf.copy()
          #Add to the end of the alignment the candidate sequence
          names_sbf_copy.append(cand)

          #Generate the alignment for subfamily
          Sbf_alignment=cand_ali[names_sbf_copy,:].trim_gaps(1)

          #Compute the scoring similarity matrix
          a=Sbf_alignment[:-1,:].score_similarity(targets=Sbf_alignment[[cand],:],gaps=opt["g"],metrics=opt["m"],weights=opt["w"])

          #Rename the index according to subfamily
          a=a.rename(index={(cand):(i)})

          #Save it to a dict
          df_dict['Candidate'].append(cand)
          df_dict['Subfamily'].append(i)
          df_dict['AWSI'].append(a['AWSI'].values[0])

  else:
    #Calling the anchor alignment & candidate alignment
    ref_ali = alignment("/software/selenoprotein_profiles/"+family+".fa")

    candidates_2 = alignment()
    for name, seq in ali:
       candidates_2.add(name+' '+ali.get_desc(name), seq)

    #Run it
    joined_ali = ref_ali.transfer_alignment(candidates_2,dont_shrink = True)
    # Removing gap columns
    joined_ali.shrink()
    # Converting alignment to Alignment from Pyaln
    cand_ali=Alignment( [ (name, joined_ali.seq_of(name)) for name in joined_ali.titles()] )

    names = cand_ali.names()
    #We want to get rid of profile sequences
    names_filtered = [item for item in names if is_selenoprofiles_output_title(cand_ali.get_desc(item)) ]

    # We select only the anchor sequences to produce subfamily alignments
    anchor_names = []
    for name in ref_ali.titles():
      anchor_names.append(name.split(' ')[0])

    df_dict = {'Candidate':[],'Subfamily':[],'AWSI':[]}
    # Compute score similarity
    for cand in names_filtered:
      anchor_names.append(cand)

      # Creating alignment of anchor+candidate
      alig_c = cand_ali[anchor_names,:]
      a = alig_c[:-1,:].score_similarity(targets=alig_c[[cand],:],gaps=[opt['g']],metrics=opt['m'],weights=opt['w'])
      # Rename the index according to subfamily
      a = a.rename(index={(cand):(family)})
      # Save it to a dict
      df_dict['Candidate'].append(cand)
      df_dict['Subfamily'].append(family)
      df_dict['AWSI'].append(a['AWSI'].values[0])

  #Save to output file the scoring similarity matrix
  score_subaln = pd.DataFrame.from_dict(df_dict)
  #Get the maximum score for each candidate
  max_idx=score_subaln.loc[score_subaln.groupby('Candidate')['AWSI'].idxmax()]
  max_idx.columns=['Candidate','Subfamily','AWSI']

  #TEMPORAL selection of candidates
  df_cand=max_idx[max_idx['Candidate'].isin(ali.names())]

  #Adding species name (Gives a warning -- Talk with Marco)
  df_cand['Species']=df_cand['Candidate'].apply(lambda x : x.split('.')[3])

  return df_cand


def main(args={}):

  if not args:
    opt = command_line_options(def_opt, help_msg)
  else:
    opt = args

  for file in opt['i']:
    file_name = os.path.basename(file)
    # Get the family : It needs the file name to be like GPx.ali, DI.ali, SelO.ali...
    fam = file_name.split(".")[0]
    # Read the alignment
    align = Alignment(file,fileformat="fasta")
    # Execute the function
    out = candidate_score_similarity(align,fam,opt)
    # Save output
    out.to_csv(opt['of']+fam+opt['o']+'.tsv',sep="\t",index=False)

if __name__ == "__main__":
    main()

