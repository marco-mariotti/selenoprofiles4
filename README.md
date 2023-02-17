# selenoprofiles4
Gene finding pipeline for selenoprotein and selenocysteine machinery.

## Installation

If you don't have conda, [install it following these instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Create a new dedicated enviroment called sp4, then activate it:

```
conda create -yn sp4
conda activate sp4
```

Install selenoprofiles4 and its dependencies in the sp4 environment:

```conda install -c mmariotti -c anaconda  -c bioconda -c biobuilds selenoprofiles4```

If everything worked correctly, the ```selenoprofiles``` command is now available, but it is still not setup. Run this and follow instructions:

```selenoprofiles -setup```

Finally, run this to automically download the latest built-in profiles ([github repo](https://github.com/marco-mariotti/selenoprotein_profiles)) to search for known selenoproteins and selenocysteine machinery:

```selenoprofiles -download```


You're ready to go! Remember to activate the sp4 environment everytime you run selenoprofiles.

## Usage
For up-to-date usage information, run:

```selenoprofiles -h```

## Bug report
Please report any bugs at the [github issues page](https://github.com/marco-mariotti/selenoprofiles4/issues).

## Citation
Selenoprofiles: profile-based scanning of eukaryotic genome sequences for selenoprotein genes.
Mariotti M, Guigó R.
Bioinformatics. 2010 Nov 1;26(21):2656-63. doi: 10.1093/bioinformatics/btq516

