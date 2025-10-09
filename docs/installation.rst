Installation
============

.. note::

  **Platform compatibility:**

  Selenoprofiles4 is officially supported on **Linux** and **macOS** systems.
  It is not compatible with Windows.
  Windows users are advised to run selenoprofiles4 using Docker or Windows Subsystem for Linux (WSL2).

Docker installation
-------------------
The docker image of selenoprofiles contains a complete installation, including all optional dependencies (see below).
To use Docker for running Selenoprofiles, you should follow the next steps:

**Pull the Docker image**:

   You can find the Docker image for Selenoprofiles on Docker Hub. To pull the image, run the following command:

   .. code-block:: bash

      docker pull maxtico/selenoprofiles_container:latest

For more information on using selenoprofiles Docker, refer to the documentation: 
https://hub.docker.com/r/maxtico/selenoprofiles_container/

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

.. warning::

  **Note for Mac users (Apple Silicon, M1/M2/M3/M4):**

  Some bioconda packages required by *selenoprofiles4* (such as *mafft*) are not yet available for the ARM64 (Apple Silicon) architecture.
  If you see an error such as::

      LibMambaUnsatisfiableError: Encountered problems while solving:
        - nothing provides mafft needed by selenoprofiles4-4.x.x-py_0

  This happens because Conda is searching only for ARM64 builds.
  To solve this, use **Rosetta**, which allows Conda to run Intel (x86_64) packages on your Mac:

  1. Make sure Rosetta is installed (only once per system). See: https://support.apple.com/en-us/102527


  2. Create and activate the environment as Intel (osx-64)

    .. code-block:: bash

        CONDA_SUBDIR=osx-64 conda create -n sp4
        conda activate sp4

  3. Install selenoprofiles4 using Intel packages

    .. code-block:: bash

        CONDA_SUBDIR=osx-64 conda install -c mmariotti -c anaconda -c bioconda -c biobuilds selenoprofiles4

  After this, installation should succeed and *selenoprofiles4* will work normally on Apple Silicon.

