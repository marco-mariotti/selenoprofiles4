Searching multiple species
==========================

Selenoprofiles was designed for comparative studies across multiples species.

If you search multiple targets (e.g. different species), you may use the same output folder (as long as target and species names are different).
As seen above, a subfolder within the output folder is created for each.

Collect results: selenoprofiles join
++++++++++++++++++++++++++++++++++++

To collect results from runs on multiple targets, you can use the *selenoprofiles join* utility. See its usage with::

  selenoprofiles join -h

As mentioned, each run of selenoprofiles produces one *ali* file for each profile with at least one result.
The *selenoprofiles join* utility collects *ali* files resulting from searching different targets, and
outputs a single *joined*  *ali* file per profile.

**Prediction identifiers** are augmented in *joined*  *ali* files to discriminate among targets.
They take this form:   *profile_name.index.label.species_name.target_file_name*

Display an overview of results: selenoprofiles drawer
+++++++++++++++++++++++++++++++++++++++++++++++++++++
  
The *selenoprofiles drawer* utility takes these *joined*  *ali* files and displays a graphical summary of results. ::

  selenoprofiles drawer -h

*Selenoprofiles drawer* also takes as input the tree of the species searched in newick format. Therefore, it shows the presence of
homologs of the protein families that were searched across the phylogenetic tree of those species.

To obtain (rough) species trees including virtually all known organisms, you may use NCBI taxonomy.
This can be done directly at its portal at http://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi,
or with more automated tools such as https://github.com/didacs/ncbi_taxonomy.

Here is an example of the output of *selenoprofiles drawer*:

.. figure:: images/drawer.example1.png 
   :width: 450

Each result is shown as a colored rectangle. A numeric tag at its left indicates its selenoprofiles numeric index.
Each column corresponds to a different profile.

The rectangle width and position indicate the prediction coverage and horizontal span when mapped in the profile alignment.
Inside each rectangle you have the chromosome (or contig) name, and the genomic coordinate boundaries,
separated with "+" for results on the plus strand, and "-" for results on the minus strand.

Finally, the intron positions as relative to the protein alignment are shown as vertical white lines.
When frameshifts are present, they are shown as vertical red lines

**Color** is used to depict the **label** of results, e.g. "selenocysteine" results are green, "cysteine" are red,
"homologue" are yellow, "pseudo" are grey. To see a list of colors, run::

  selenoprofiles drawer -colors
	   
*Selenoprofiles drawer* becomes more concise with option *-a*, as it shows only the count of results for each label:
	   
.. figure:: images/drawer.example2.png 
   :width: 450

