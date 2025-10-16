.. _docker_usage:

Using selenoprofiles in Docker
==============================

Pull the image
--------------
The Docker image of Selenoprofiles provides a complete, ready-to-use environment with all optional dependencies already installed.
This is the recommended way to run Selenoprofiles if you want to avoid managing dependencies:

  The official Docker image is available on Docker Hub:

  .. code-block:: bash

    docker pull maxtico/selenoprofiles_container:latest

Run Selenoprofiles in Docker
----------------------------
You can run Selenoprofiles directly within the container.
In most cases, you'll want to mount your working directory so results and input files can be accessed outside of the container. 
Here's a complete example command:

  .. code-block:: bash

    docker run --rm -v $(pwd):/mnt maxtico/selenoprofiles_container:latest \
      selenoprofiles \
      -o /mnt/{selenoprofiles_output_folder} \
      -t /mnt/{species_genome} \
      -s {species} \
      -p {profile} \
      -temp /mnt/{temp_folder}

Let's go step by step with a simple example.

**Searching for selenoproteins in the human genome**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Let’s imagine we are working with the human genome and we want to check for the presence of *SELENOP* 
and *SELENOS* using Selenoprofiles. First, make sure your genome file (e.g. genome.fa) is in your current working directory.
Then, run the following command:

  .. code-block:: bash
  
    docker run --rm -v $(pwd):/mnt maxtico/selenoprofiles_container:latest \
        selenoprofiles \
        -o /mnt/selenoprofiles_out \
        -t /mnt/genome.fa \
        -s Homo_sapiens \
        -p GPX,DIO \
        -temp /mnt/temp
        -output_gtf_file /mnt/homo_sapiens.gtf

The Docker image also includes all Selenoprofiles utilities, which can be run in the same way.

**Collecting results from runs on multiple targets**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When running Selenoprofiles on several target genomes (for example, across different species), each run produces 
one or more *.ali* files — one for every profile where at least one hit was found. To combine the results from different 
runs into a single set of alignments, you can use the built-in **join** utility.

For example, suppose you have already run Selenoprofiles on several genomes and stored all outputs under `selenoprofiles_out`.
You can merge the alignments for specific profiles (e.g. *GPX* and *DIO*) as follows:

  .. code-block:: bash
    
      docker run --rm -v $(pwd):/mnt maxtico/selenoprofiles_container:latest \
          selenoprofiles \
          join \
          -d /mnt/selenoprofiles_out \
          -p GPX,DIO \
          -o /mnt/join_results

**Classifying predictions based on orthologous groups**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
After merging the results from multiple runs with the `join` utility, you can classify the candidate sequences 
into **orthologous subfamilies** using the **`selenoprofiles orthology`** utility.

For example, continuing from the previous step where we merged profiles *GPX* and *DIO* into 
`join_results/`, you can run:

  .. code-block:: bash

    docker run --rm -v $(pwd):/mnt maxtico/selenoprofiles_container:latest \
        selenoprofiles \
        orthology \
        -i /mnt/join_results/GPX.ali /mnt/join_results/DIO.ali \
        -of /mnt/orthology_results

**Refining selenoprotein predictions**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
After classifying candidate sequences with **`selenoprofiles orthology`**, you may want to 
**filter out predictions that are not expected** in a given lineage. This is especially useful to remove false 
positives or sequences that correspond to pseudogenes or retrotransposons.

The **`selenoprofiles lineage`** utility performs this filtering. It takes the `.tsv` files produced 
by `orthology` as input and outputs one filtered `.tsv` file per input, excluding non-expected genes 
based on lineage expectations.

Continuing our tutorial example:

  .. code-block:: bash

    docker run --rm -v $(pwd):/mnt maxtico/selenoprofiles_container:latest \
        selenoprofiles \
        lineage \
        -i /mnt/orthology_results/GPX.orthology.tsv /mnt/orthology_results/DIO.orthology.tsv \
        -o /mnt/lineage_results

.. warning::

  When `lineage` utility is executed without `-map` option, the program uses a pre-installed module called `*ncbi_db*`.
  This happens because the tool needs to get the lineage for each species in the input file. 
  The `*ncbi_db*` program requires a database file to infer the lineage. Thus, before running the `lineage` utility,
  you need to follow the following workflow:

  * 1. Create a db file (*tax.db*) using **ncbi_taxonomy_tree**:

    .. code-block:: bash

      docker run --rm -v $(pwd):/mnt maxtico/selenoprofiles_container:latest \
          ncbi_taxonomy_tree \
          --makedb \
          -o /mnt/
  
  * 2. Configure the ncbi_config file located in your home directory (e.g. ~/.ncbi_config) to point to the db file created in step 1. You can do this by changing the following line to the file:

    .. code-block:: text

      ncbi_taxdb = /path/to/your/tax.db
  
  * 3. Now you can run the `lineage` utility as shown above with a slight modification to the command:

    .. code-block:: bash

      docker run --rm -v $(pwd):/mnt -v ~/.ncbi_config:/root/.ncbi_config maxtico/selenoprofiles_container:latest \
          selenoprofiles \
          lineage \
          -i /mnt/orthology_results/GPX.orthology.tsv /mnt/orthology_results/DIO.orthology.tsv \
          -o /mnt/lineage_results \
          -map ncbi_db      

**Assessing selenoprotein annotations**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Additionally, we can include in our workflow a final step to evaluate genome annotations, using **`selenoprofiles assess`**.  
This utility compares existing genome annotations with the **Selenoprofiles predictions**,
producing a comprehensive overview of how well the genome is annotated for selenoproteins.

Countinuing our tutorial example, let's assess how well the human genome is annotated for selenoproteins (note we require
the genome annotation file `genome.gff3` in addition to the genome sequence `genome.fa`):

  .. code-block:: bash

    docker run --rm -v $(pwd):/mnt maxtico/selenoprofiles_container:latest \
        selenoprofiles \
        assess \
        -s /mnt/homo_sapiens.gtf \
        -e /mnt/genome.gff3 \
        -f /mnt/genome.fa \
        -o /mnt/homo_sapiens.complete.tsv \
        -agg /mnt/homo_sapiens.aggregate.tsv


For more information regarding the selenoprofiles4 Docker, refer to the documentation: 
https://hub.docker.com/r/maxtico/selenoprofiles_container/

