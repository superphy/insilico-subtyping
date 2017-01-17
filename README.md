PHYLOTYPER
==========

Phylotyper predicts biological subtypes from gene sequence data. It
build a reference phylogenetic tree from genes with known subtype and
then uses ancestral reconstruction to assign likelihoods of each subtype
to the branch points in the tree. A new unknown sequence can be assigned
a subtype based on the extrapolated value from its ancestors in the
tree.

https://github.com/superphy/insilico-subtyping

INSTALLATION
============

Phylotyper is a python package that uses several external programs.
Phylotyper requires the following external programs:

1.  FastTree (tested with version 2.1) -
    http://www.microbesonline.org/fasttree/
2.  MAFFT (tested with version 7.1) -
    http://mafft.cbrc.jp/alignment/software/
3.  trimAI (tested with version 1.2) - http://trimal.cgenomics.org/
4.  BLAST+ (tested with version 2.2.28) -
    ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
5.  R (tested with version 3.3.2) - https://cran.r-project.org/

Phylotyper was built and tested with python 2.7. To install Phylotyper
using python virtual environments:

-   `virtualenv -p python2.7 phylotyper` replacing `phylotyper` with the
    directory name you want to use for phylotyper. Use this directory in
    subsequent steps.
-   `phylotyper/bin/activate`. Activate the virtual enviroment.
-   `git clone https://github.com/superphy/insilico-subtyping.git .`.
    Clone the phylotyper git repository into the directory.
-   `cd phylotyper/`. Move into the phylotyper directory.
-   `python setup.py install`. Install dependencies and create
    `phylotyper` executable.

Running tests:

-   After installing external packages, create config file. See section
    [Setting up Phylotyper config
    file](#setting-up-phylotyper-config-file)
-   `PHYLOTYPER_CONFIG=/path/to/config/file python setup.py test`. Run
    tests.

RUNNING PHYLOTYPER
==================

Setting up Phylotyper config file
---------------------------------

Phylotyper has user-defined settings and also requires several external
programs. These are provided to Phylotyper through a INI-config file.
The config file must contain the following settings:

~~~~ {include="phylotyper_example.ini"}
; Sample Phylotyper config file
;

# External programs
[external]
fasttree=/path/to/FastTree
mafft=/path/to/mafft
trimal=/path/to/trimal
makeblastdb=/path/to/makeblastdb ; part of the BLAST+ suite
blastn=/path/to/blastn ; part of the BLAST+ suite
blastx=/path/to/blastx ; part of the BLAST+ suite
blastdbcmd=/path/to/blastdbcmd ; part of the BLAST+ suite

# R settings
[R]
lib=/path/to/local/R/libs/ ; Location where downloaed R packages can be saved
repo=http://cran.stat.sfu.ca/ ; HTTP address of R repository to use for downloading R packages
rscript=/usr/bin/Rscript ; Part of the R suite

# Phylotyper settings
[phylotyper]
prediction_threshold=0.9 ; Cutoff for calling subtype assignments
~~~~

Set up environment variable
`export PHYLOTYPER_CONFIG=/path/to/your/phylotyper/config/file` to
indicate location of config file to Phylotyper. To make this variable
persistent, the line
`export PHYLOTYPER_CONFIG=/path/to/your/phylotyper/config/file` can be
added to your \~/.bashrc.

Alternatively, the location of the config file can be set using the
argument `--config` when running `phylotyper`.

Running Phylotyper using built-in subtype scheme
------------------------------------------------

Adding new subtype scheme to Phylotyper
---------------------------------------

ABOUT PHYLOTYPER
================

CONTACT
=======

Matt Whiteside matthew.whiteside@phac-aspc.gc.ca

TODO
====

1.  Get tests working
2.  Install stx1
3.  Install eae
4.  Install flic
5.  Install wyz

