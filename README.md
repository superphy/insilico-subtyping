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

See section [Setting up Phylotyper config
file](#setting-up-phylotyper-config-file).

Phylotyper can first identify putative subtype loci in input genomes and
then use the detected loci to predict subtypes. To run phylotyper with
genome DNA sequences as input:

`phylotyper genome subtype_scheme output_directory genome_input1.fasta [genome_input2.fasta genome_input3.fasta...][--noplots][--config /path/to/config/file][--index /path/to/index/file]`

### Arguments:

1.  `genome` - Sub-command for running phylotyper on genome input. Other
    sub-commands are `new` for building new subtype schemes.
2.  `subtype_scheme` - The name for the subtype scheme. See [Built-in
    subtype schemes](#built-in-subtype-schemes) for a list of the
    subtypes currently packaged with phylotyper. E.g. value `stx2`.
3.  `output_directory` - Results and graphics will be output to this
    directory. Files will be overwritten.
4.  `genome_input1.fasta` - Genome multi-fasta DNA sequence input.
    Multiple genomes should be provided as separate files.

### Options:

1.  `--noplots` - Turn off figures.
2.  `--config /path/to/file` - Phylotyper config file. Can also be set
    in environment variable `PHYLOTYPER_CONFIG`. See [Setting up
    Phylotyper config file](#setting-up-phylotyper-config-file).
3.  `--index /path/to/file` - YAML-format file that lists the locations
    of subtype reference files. You can specify a file that is not the
    default file packaged with phylotyper.

Built-in subtype schemes
------------------------

Listed below are the subtypes, their sequence type, number of loci used
in the prediction and a description of subtype. Inputs for genomes is
always nucleotide multi-fasta sequences. Detected loci in genomes are
translated if the subtype sequence type is amino acid. If the input are
gene sequences, the correct sequence type and number of loci are
required.

~~~~ {include="available_subtypes.md"}
- stx2
..+ sequence type: aa
..+ number of loci: 2
..+ description: Escherichia coli Stx2 subtypes
~~~~

Adding new subtype scheme to Phylotyper
---------------------------------------

See section [Setting up Phylotyper config
file](#setting-up-phylotyper-config-file).

New subtype schemes can be added to phylotyper. The reference inputs
will be processed and saved in the locations specified in the
YAML-format file: `subtypes_index.yaml`. Optionally, you can define a
non-default subtype data directory by providing your own YAML file
using `--index /path/to/index/file`.

To create a new subtype for use in phylotyper:

`phylotyper new subtype_scheme subtype_assignment_file output_directory reference_loci1.fasta [reference_loci2.fasta reference_loci3.fasta...][--aa][--config /path/to/config/file][--index /path/to/index/file][--description "Help description"]`

### Arguments:

1.  `new` - Sub-command for adding new subtype in phylotyper. Other
    sub-commands are `genome` for running phylotyper on genome input.
2.  `subtype_scheme` - The name for the subtype scheme.
3.  `subtype_assignment_file` - Subtypes for input loci sequences.
4.  `output_directory` - Results and graphics will be output to this
    directory. Files will be overwritten.
5.  `reference_loci1.fasta` - Fasta DNA or amino-acid sequence input.
    Multiple loci should be provided as separate files. Each entry in
    each fasta file should have a subtype assignment in the
    `subtype_assignment_file`.

### Options:

1.  `--aa` - Set flag when input is amino-acid sequences. Default is
    nucleotide.
2.  `--config /path/to/file` - Phylotyper config file. Can also be set
    in environment variable `PHYLOTYPER_CONFIG`. See [Setting up
    Phylotyper config file](#setting-up-phylotyper-config-file).
3.  `--index /path/to/file` - YAML-format file that lists the locations
    of subtype reference files. You can specify a file that is not the
    default file packaged with phylotyper.
4.  `--description` - A help description for the subtype scheme.

If you would like to contribute a subtype scheme to the main repository, please contact us.

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

