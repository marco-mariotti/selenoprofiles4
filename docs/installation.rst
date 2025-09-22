Installation
============

Docker installation
-------------------
Docker installation contains all information available from selenoprofiles, including optional dependencies (see below).
To use Docker for running Selenoprofiles, you should follow the next steps:

1. **Pull the Docker image**:

   You can find the Docker image for Selenoprofiles on Docker Hub. To pull the image, run the following command:

   .. code-block:: bash

      docker pull maxtico/selenoprofiles_container:latest

2. **Running Selenoprofiles on the Docker container**

  Once you've pulled the Docker image, you can run the Selenoprofiles package inside the Docker container.
  In most cases, you’ll want to mount your working directory so results and input files can be accessed outside of the container. Here’s a complete example command:

  .. code-block:: bash

    docker run --rm -v $(pwd):/mnt maxtico/selenoprofiles_container:latest \
      selenoprofiles \
      -o /mnt/{selenoprofiles_output_folder} \
      -t /mnt/{species_genome} \
      -s {species} \
      -p {profile} \
      -output_gtf_file {gtf_file} \
      -temp /mnt/{temp_folder}


For more information on using selenoprofiles Docker, refer to the documentation: 
https://hub.docker.com/repository/docker/maxtico/selenoprofiles_container/

Conda installation
------------------

We recommend to use the conda package manager to install selenoprofiles4
(check `this page to install conda <https://docs.conda.io/en/latest/miniconda.html>`_).

We recommend to create a new dedicated enviroment, e.g. called sp4, then activate it. In a terminal, run::

  conda create -yn sp4
  conda activate sp4

Then **install selenoprofiles4** and its dependencies in the sp4 environment::

    conda install -c mmariotti -c anaconda  -c bioconda -c biobuilds selenoprofiles4

If everything worked correctly, the selenoprofiles command is now available, but it is still not setup.
Run this and follow instructions::
  
  selenoprofiles -setup

Finally, run this to automically download the latest built-in profiles 
to search for known selenoproteins and selenium markers::

  selenoprofiles -download


A couple of notes:

 - the ``-setup`` command creates the config file: *~/.selenoprofiles_config.txt* -- you may edit this to change default options
 - the ``-download`` command downloads data in the folder *~/selenoprofiles_data/*. If necessary, you may move it to another location and change the value of *selenoprofiles_data_dir* in the config file.
  
You're ready to go! Remember to activate the sp4 environment everytime you run selenoprofiles.

At any time, run this to access the help page::

  selenoprofiles -h

Check the :doc:`get_started` page to start using selenoprofiles.


Optional dependencies for selenoprofiles utilities
--------------------------------------------------

Selenoprofiles comes with some utilities: join, build, assess, orthology, lineage, drawer, database.
Some of them have additional dependencies which are not strictly required for selenoprofiles and are not automatically installed:

   .. code-block::

     pip install selenoprofiles4[addons]

