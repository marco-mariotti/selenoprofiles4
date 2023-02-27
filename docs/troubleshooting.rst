Troubleshooting
===============

If you have problems running selenoprofiles, make sure to inspect the full list of its options to see if there's anything relevant::

  selenoprofiles -h full

A common problem is a database lockdown. Selenoprofiles uses a sqlite database to store results.
If killed, sometimes such database gets locked. If this happens, run::

  selenoprofiles database -i output_folder/species_name.target_file_name/results.sqlite -clean

Bug reports
-----------

If you cannot find solution to a problem by yourself, you are encouraged to file an issue at
`its github page <https://github.com/marco-mariotti/selenoprofiles4/issues>`_.
You can also use the same system to file feature requests.


