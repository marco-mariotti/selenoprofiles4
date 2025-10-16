
Welcome to selenoprofiles's documentation!
==========================================

Selenoprofiles is a profile-based gene finding pipeline specialized in selenoprotein and selenocysteine machinery.

The program takes two inputs per run:

 - one or more *profile* alignments = the protein families to search for,
 - a genome (or any other nucleotide database) = the *target* you want to scan.

By combining the homology-based gene finding tools blast, exonerate, and genewise, selenoprofiles identifies all genes
in the target that are homologous to the input profile(s).

Selenoprofiles comes with built-in profiles for selenoproteins and other genes related to selenium metabolism, allowing
out-of-the-box characterization of the *selenoproteome* of any organism, prokaryotic or eukaryotic.

However, selenoprofiles may be used to search for any protein family (i.e. unrelated to selenoproteins):
users may easily build their own profile(s).

To install selenoprofiles, get started or learn advanced usage, check the documentation pages below.

.. toctree::
   :maxdepth: 1
   :caption: Contents:
      
   installation
   get_started
   citation
   utilities
   workflow
   user_profiles
   pipeline_detailed
   advanced
   troubleshooting
   docker
