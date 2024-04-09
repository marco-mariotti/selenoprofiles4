Installation
------------

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
++++++++++++++++++++++++++++++++++++++++++++++++++

Selenoprofiles comes with some utilities: join, build, drawer, database.
Some of them have additional dependencies which are not strictly required for selenoprofiles and are not automatically installed:

 - selenoprofiles drawer requires ete3. Install it with::

     conda install -c etetoolkit ete3

 - selenoprofiles assess requires pyranges and pyfaidx. Install them with::

     pip install pyranges

     conda install -c bioconda pyfaidx

 - selenoprofiles build offers a graphical interface which requires pylab. Install it with::

     conda install matplotlib

     
