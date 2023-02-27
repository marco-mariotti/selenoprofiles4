Workflow
========

In summary, the pipeline is structured as follows (see also figure below).

The program **psitblastn** is used with a PSSM derived from the profile alignment to identify matches in the target genome.
These matches are then used, through the two splice alignment programs **exonerate** and **genewise**,
to deduce the exonic structure of the candidate genes.

The predictions of these three programs
are analyzed and the best prediction is **chosen**. Then, it is then **labelled** through a dedicated procedure.
The predictions are also sometimes **improved** by a few procedures, e.g. they are completed at the 3' looking for a stop codon
and the 5' looking for a methionine codon.

In the workflow, there are three **filtering** steps. First, the blast filtering, which controls how many gene candidates will be processed.
Then the p2g filtering and p2g refiltering, both of which are at the end of the pipeline.
All filtering steps are user definable for each profile.
We provide a sensible default filtering for user input families: each alignment is examined and, based on its sequence conservation,
a similarity threshold is chosen (AWSI filter). This means that a very conserved profile outputs only very similar sequences.

A third filter is applied when multiple profiles are searched, and there are overlapping matches across profiles.
These are assigned to the highest scoring family, and the others are dismissed.

A graphical summary:

.. figure:: images/workflow.pdf

Selenoprofiles normally performs the full pipeline, taking care of
skipping the steps executed previously. The steps of selenoprofiles are:
blast, exonerate, genewise, prediction choice, prediction filtering,
database storage, output. These are denoted respectively by the
step-options *-B -E -G -C -F -D -O* (see figure).

After the filtering step, results are stored in a SQLite **database**.
When selenoprofiles is run, it checks first if the results database contain already the 
results requested, and in that case it passes directly to the output step.

If the
user specify any step-option, the execution of the corresponding step
and of all next ones is forced. This is necessary if you changed
parameters or profile specific procedures. If for example you changed
some parameter relative to the filtering phase, you must force filtering
and output with *-F*.

**Important:** when output is forced,
selenoprofiles  overwrites previous output files, but **it will never delete any**.
This may lead to overlapping predictions in the output.
If you re-run selenoprofiles adding profiles to the search, we recommend to delete all output files first. 


.. warning::

   In the selenoprofiles paper, the program SECISearch was included in the workflow,
   but it is not anymore. To run its successor SECISearch3, visit the webserver at https://seblastian.crg.es/.
   We recommend to produce *-out_three_prime* output with selenoprofiles to get
   the region downstream of predictions. Collect the *.three_prime* files
   corresponding to "selenocysteine" labels, upload them to the webserver, and select
   "SECIS prediction" only. 
   

	    
Full description of the pipeline
--------------------------------

The following detailed description is only meant for advanced users.

TBD
