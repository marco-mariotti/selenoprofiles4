Getting started
---------------


Command line structure
======================

In a single run, selenoprofiles can search for multiple profiles in a single target.
Here's a minimal command line::

  selenoprofiles  -o output_folder  -t target_file  -s "species name"  -p profile  [options]

These are the compulsory arguments:

 * **\-o**  output folder, will be created if non-existing. 
 * **\-t**  target_file = a (multi-)fasta file containing nucleotide sequences
 * **\-s**  a species descriptor with no restrictions. For multiple words, use quotes 
 * **\-p**  the profile(s) to be searched. Multiple comma-separated arguments are accepted. Each argument can be:

   * a profile name: invokes a built-in profile (located in the profiles_folder defined in the config file)
   * a profile set: a keyword expanded to a list of profiles. Profile sets are defined in the config file
   * a path to a profile alignment: to create one, see: selenoprofiles build -h

By default, selenoprofiles assumes the target is an eukaryotic genome sequence, and it will attempt to predict
introns. If you're searching (intronless) prokaryotes or eukaryotic mRNA sequences, use option *-no_splice*.

Blast searches can use multiple CPUs. To control how many, use option *-ncpus*. This is just one of the many non-compulsory options and parameters.
The config file *~/.selenoprofiles_config.txt* defines default values for all of them. The full list of options can be inspected with::

  selenoprofiles -h full

Searching for selenoproteins and selenium markers
=================================================

The common usage of selenoprofiles is the prediction of selenoprotein families and selenium metabolism genes.
For this task, you must specify appropriate profile sets: the selenoprotein families expected in a genome depend
on its taxonomy.

To search for metazoan selenoprotein families and selenium usage markers, use::
  
  -p metazoa,machinery

To search for all eukaryotic families, use::
  
  -p eukarya

To search for all prokaryotic families, use::
  
  -p prokarya


Selenoprofiles output
=====================

Upon completion, you will find output files in a subfolder of the output folder::

  output_folder/species_name.target_file_name/output/

Selenoprofiles produces one file per gene prediction, per requested format.
The output files are named after **prediction identifier**, which have the syntax::
  
  profile_name.index.label

Where:

 - *profile_name* identifies the source profile for the prediction
 - *index* is an arbitrary numeric id
 - *label* identifies the class of predicted gene.

For selenoprotein families, the label can be:
   
    - "selenocysteine" for selenoproteins: UGA is aligned to the Sec position of the profile
    - "cysteine" for cysteine-homologs of selenoproteins (i.e. Cys aligned to Sec position)
    - if any other amino acid is aligned to Sec, the label takes its name
    - for predictions that do not include the Sec position of the profile:

      - "unaligned" if the prediction ends upstream or starts downstream of the Sec position
      - "gapped" if Sec is not aligned but there are homologous regions on both sides
      
    - "pseudo" for predictions containing inframe stops or frameshifts (likely pseudogenes)
    - "uga_containing" for predictions whose only pseudogene features are in-frame UGAs

For profiles that are not selenoprotein families (i.e. do not contain Sec), the label can be:

    - "homologue" for standard predictions
    - "pseudo" (see above)

Output files are named after the prediction identifier and the format::

  profile_name.index.label.format   # e.g. GPx.1.selenocysteine.p2g
      
These output **formats** are supported:

 - *p2g*:         native format with query-target alignment, coordinates and other info; :download:`see an example here<files/example.p2g>`.
 - *fasta*:       fasta file with the predicted protein sequence
 - *gff*:         gff file with the genomic coordinates of the prediction
 - *gtf*:         analog to gff, but with this last field: gene_id "prediction_id"; transcript_id "prediction_id";
 - *cds*:         fasta file with coding sequence, including the stop codon if applicable
 - *three_prime*: fasta with the sequence immediately at 3' of the prediction. Length is defined by -three_prime_length
 - *five_prime*:  fasta with the sequence immediately at 5' of the prediction. Length is defined by -five_prime_length
 - *dna*:         fasta file with the full nucleotide sequence, including introns (and frameshift-causing insertions if any)
 - *introns*:     fasta file with the sequence of the introns, split into different fasta headers      
      
Additionally, a fasta alignment called *profile_name.ali* is created.
Only one such *ali* file is produced per profile, containing the sequences of all predictions plus the profile sequences.

On the command line, option -output_FORMAT activates the corresponding output for each prediction, e.g. -output_gff will produce gff files.
By default, only the ali and p2g formats are active, as visible in the config file::

  ### active output format
  output_ali=1
  output_p2g=1

To create a single output file for all predictions, use -output_FORMAT_file providing as argument the file that will be created,
e.g. ``-output_fasta_file all_predicted_proteins.fa``


