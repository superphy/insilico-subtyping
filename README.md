PHYLOTYPER
==========

Phylotyper predicts biological subtypes from gene sequence data. It
build a reference phylogenetic tree from genes with known subtype and
then uses ancestral reconstruction to assign likelihoods of each subtype
to the branch points in the tree. A new unknown sequence can be assigned
a subtype based on the extrapolated value from its ancestors in the
tree.

https://github.com/superphy/insilico-subtyping

CONTENTS
========

1. [Examples](#examples)
1. [Installation](#installation)
1. [Running Phylotyper](#running-phylotyper)
    1. [Setting up Phylotyper config file](#setting-up-phylotyper-config-file)
    1. [Running Phylotyper using built-in subtype scheme](#running-phylotyper-using-built-in-subtype-scheme)
    1. [Built-in subtype schemes](#built-in-subtype-schemes)
    1. [Adding new subtype scheme to Phylotyper](#adding-new-subtype-scheme-to-phylotyper)
1. [About Phylotyper](#about-phylotyper)
1. [Contact](#contact)

EXAMPLES
========

After [installing](#installation) phylotyper and its dependencies, here is a quick example to demonstrate its usage:

1. Inform phylotyper of the locations for the dependencies by setting up a [INI config file](https://raw.githubusercontent.com/superphy/insilico-subtyping/master/phylotyper_example.ini)

Details [here](#setting-up-phylotyper-config-file).

2. Predict a subtype for one of the amino acid schemes packaged in phylotyper:

    `phylotyper genome stx2 example_data/output/ example_data/genome.fasta \`
    `--config yourconfigfile.ini`

See `example_data/output/` for the results. Find the full list of available schemes and their sequence types [here](#built-in-subtype-schemes). Details on the results can be found [here](#running-phylotyper-using-built-in-subtype-scheme). Note: the config file can also be specified using the enviromentment variable `PHYLTYPER_CONFIG` instead of the `--config` option.

3. Add your own DNA subtype scheme called `myexample`:

    `phylotyper new myexample phylotyper/example_data/example_subtypes.tsv phylotyper/example_data/output/ \`
    `phylotyper/example_data/dna_example_genes.fasta --config yourconfigfile.ini`

4. Use this new DNA subtype scheme to predict subtypes:

    `PHYLOTYPER_CONFIG=yourconfigfile.ini phylotyper genome myexample example_data/output/ example_data/genome.fasta`

See `example_data/output/` for the results. Details on the results can be found [here](#running-phylotyper-using-built-in-subtype-scheme). Note: the config file can also be specified using the enviromentment variable `PHYLTYPER_CONFIG` instead of the `--config` option.

5. Add your own multi-loci amino acid subtype scheme called `myexample2`:

    `phylotyper new myexample2 phylotyper/example_data/example_subtypes.tsv phylotyper/example_data/output/ \`
    `phylotyper/example_data/aa_example_genes_loci1.fasta phylotyper/example_data/aa_example_genes_loci2.fasta --aa \`
    `--config yourconfigfile.ini`

6. Use this new amino-acid subtype scheme to predict subtypes:

    `phylotyper genome myexample2 example_data/output/ example_data/genome.fasta --config yourconfigfile.ini`

See `example_data/output/` for the results.

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

-   `mkdir phylotyper; git clone https://github.com/superphy/insilico-subtyping.git phylotyper`.
    Clone the phylotyper git repository into the directory replacing
    `phylotyper` with the directory name you want to use for phylotyper.
    Use this directory in subsequent steps.
-   `virtualenv -p python2.7 phylotyper`. Install virtual environment.
-   `source phylotyper/bin/activate`. Activate the virtual enviroment.
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

### Outputs:

The output directory is defined in the command-line arguments. The
results of the phylotyper run is given in the file:

`<output_directory>/subtype_predictions.tsv`

This tab-delimited result file contains the following columns:

1.  `genome`: Genome name
2.  `tree_label`: Unique label for a given gene copy in a genome
3.  `subtype probability`: The marginal likelihood from the Phylotyper
    analysis (or `identical` if an identical match was found in the
    reference set)
4.  `phylotyper_assignment`: The subtype assignment provided that the
    probability is above the pre-defined cutoff
5.  `loci`: A list of subtype genes found in the genome. Each list item
    contains a two-part tuple: 1. lcl|genome|unique\_name followed by 2.
    lcl|genome|contig:start-stop.

Additional analysis files that are output include:

1.  `<genome>_loci<nloci>_step2_alignment_input.fasta` Fasta-file
    containing subtype gene found in input genome by BLAST step. There
    can be multiple copies/alleles for a single gene loci (If there is
    only one genome and one loci, the file will be called
    '<input_filename>.locus1\`)
2.  `<genome>_step3_alignment_trimming_summary.html` HTML output from
    trimAI indicating trimmed columns in the alignment
3.  `<genome>_step4_profile_alignment_output.fasta` Fasta-file
    containing aligned input genes and reference genes. MAFFT is used
    for the alignment.
4.  `<genome>_step5_subtype_tree.newick` Newick-file containg
    phylogenetic tree. FastTree is used to build the tree.
5.  `<genome>_step5_posterior_probability_tree.png` Image file showing
    the phylogenetic tree and marginal likelihoods for the unknown gene.

Built-in subtype schemes
------------------------

Listed below are the subtypes, their sequence type, number of loci used
in the prediction and a description of subtype. Inputs for genomes is
always nucleotide multi-fasta sequences. Detected loci in genomes are
translated if the subtype sequence type is amino acid. If the input are
gene sequences, the correct sequence type and number of loci are
required.

~~~~ {include="available_subtypes.md"}
- stx1
..+ sequence type: nt
..+ number of loci: 2
..+ description: Escherichia coli Shiga-toxin 1 (Stx2) subtype
- stx2
..+ sequence type: aa
..+ number of loci: 2
..+ description: Escherichia coli Shiga-toxin 2 (Stx2) subtype
- eae
..+ sequence type: nt
..+ number of loci: 1
..+ description: Escherichia coli Initimin (eae) subtype
- wz
..+ sequence type: nt
..+ number of loci: 2
..+ description: Escherichia coli O-serotype based on wzy and wzx genes
- flic
..+ sequence type: nt
..+ number of loci: 1
..+ description: Escherichia coli H-serotype based on fliC gene
~~~~

Adding new subtype scheme to Phylotyper
---------------------------------------

See section [Setting up Phylotyper config
file](#setting-up-phylotyper-config-file).

New subtype schemes can be added to phylotyper. The reference inputs
will be processed and saved in the locations specified in the
YAML-format file: `subtypes_index.yaml`. Optionally, you can define a
non-default subtype data directory by providing your own YAML file using
`--index /path/to/index/file`.

To create a new subtype for use in phylotyper:

`phylotyper new subtype_scheme subtype_assignment_file output_directory reference_loci1.fasta [reference_loci2.fasta reference_loci3.fasta...][--aa][--config /path/to/config/file][--index /path/to/index/file][--description "Help description"]`

### Arguments:

1.  `new` - Sub-command for adding new subtype in phylotyper. Other
    sub-commands are `genome` for running phylotyper on genome input.
2.  `subtype_scheme` - The name for the subtype scheme.
3.  `subtype_assignment_file` - Subtypes for input loci sequences. This
    is a tab-delimited file with genomes in column 1 and subtypes in
    column 2.
4.  `output_directory` - Results and graphics will be output to this
    directory. Files will be overwritten.
5.  `reference_loci1.fasta` - Fasta DNA or amino-acid sequence input.
    Multiple loci should be provided as separate files. Each genome in
    the fasta files should have a subtype assignment in the
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

If you would like to contribute a subtype scheme to the main repository,
please contact us.

### Input Formats for New Pipeline:

Input sequence files are required to be in fasta format. If the subtype
scheme involves multiple gene/loci, there should be one file per gene
with all copies of the gene in the same file. The order of the gene
files determines their order in the superalignment. A genome can have
multiple copies of a gene/loci. To indicate a genome with multiple
copies/alleles of the same gene, format the fasta header as follows:

`>lcl|genome_identifier|allele_identifier`

or

`>genome_identifier|allele_identifier`

The subtype file is a tab-delimited text file containing two colums:

`genome_identifier<\tab>subtype`

ABOUT PHYLOTYPER
================

### Ancestral State Reconstruction

### Steps in the Phylotyper pipeline

1.  Identify subtype gene loci in input genomes using blastn or blastx
    depending on subtype sequence type. For each gene loci found in the
    input genomes:
2.  Compare input gene against reference genes. If an identical
    reference gene found, report subtype and terminate at this point. If
    no identical sequences found, proceed with Phylotyper analysis.
3.  Align input genes against a pre-aligned set of reference genes using
    the tool MAFFT's `--add` option.
4.  Automatically trim alignment using tool trimAI.
5.  If multiple loci are involved, concatenate individual alignments
    into superalignment.
6.  Generate maximum likelihood phylogenetic tree of aligned genes with
    FastTree.
7.  Run phytools `rerootingMethod` using the phylogenetic tree and
    assigned subtypes.
8.  Identify the subtype with maximum marginal likelihood for the
    unknown gene.

### Adding new subtype schemes

Phylotyper includes functionality to allow you to add your own subtype
schemes. Schemes can be one or more loci (currently no limit defined,
however, we have only tested with two loci). The sequences can be
nucleotide or amino acid. The new subtype pipeline automatically
generates and stores all required input files, so that future subtype
prediction runs only have to reference a subtype name.

#### Reference files

Phylotyper reference files and their locations for available subtype
schemes are defined in a YAML-format file. The default subtype file is
`<package_name>/phylotyper/subtype_index.yaml`. This file can be
specified using the `--index` option. The `subtype_index.yaml` file is
automatically updated when new subtypes are created. In
`subtype_index.yaml`, the `root` field specifies the parent data
directory. If a relative path is provided, the data directory is a
subdirectory of `<package_name>/phylotyper/`. The default is
`<package_name>/phylotyper/data`. Each subtype scheme will have a
subdirectory under this `root` data directory. The `subtypes` field in
`subtype_index.yaml` is a list of all subtypes indexed by name. Under a
subtype the following files/options will be defined:

1.  alignment: the reference sequence alignment
2.  desc: Verbal description of subtype scheme defined in `new`
    arguments
3.  lookup: JSON object containing all reference sequences
4.  nloci: Number of loci in scheme
5.  rate\_matrix: emperical transition matric for the Mk model
6.  search\_database: BLAST database for searching genomes
7.  seq: sequence type (nt|aa)
8.  subtype: tab-delimited file listing genomes (column 1) and subtypes
    (column 2)

#### Multiple loci

Multiple loci can be used in a Phylotyper subtype scheme. For example,
schemes `stx1`, `stx2`, `wz` all use two genes. The individual loci will
be BLAST'd and aligned independently. The individual loci alignments are
concatenated to form a superalignment. The superalignment is used as
input into the phylogenetic tree building step.

#### Parameterization of the transition matrix

Phylotyper offers flexibility in the parameterization of evolution model
used in the ASR step. A component of the underlying ASR framework, is an
Mk or markov model of subtype evolution, for which an emperical
transition rate matrix is estimated from the data. The transition matrix
is used to calculate the expected number of subtype state changes given
a distance in the phylogenetic tree. Different model parameterizations
can be defined for the transition rate matrix. The simplest
parameterization available in Phylotyper is the equal rates model; all
subtypes have the same forward and reverse rate. The most complex
parameterization available in Phylotyper is the symmetric model, wherein
each forward and reverse rate for a given pair of subtypes are assigned
a separate parameter. Frequently, the number of subtypes makes the
symmetric model too computationally prohibitive (it is unavailable for
schemes with \>10 subtypes). To offer more flexible models in these
situations with reduced numbers of free parameters, two custom
parameterization approaches were developed. The custom approaches both
use a binning strategy that attempts to identify sets of subtypes that
would have similar rates and assign them a single parameter as a set.

1.  Small-distance binning: The rationale is to select the closest
    subtypes in the phylogenetic tree as free parameters. Subtypes with
    large distances are assigned a common parameter. Maximum
    inter-patristic distances are collected for all subtypes and
    modelled as a normal mixture distribution. The smallest normal
    distribution is selected (based on mean) and all pairs of subtypes
    belonging to this distribution are set as a free parameters. Other
    subtypes are assigned a common parameter.

2.  Iterateive binning: The rationale is to approximate the transition
    rates individually for each subtype pair, cluster subtypes by the
    approximate rates and then assign each cluster a separate parameter.
    Each subtype pair (forward and reverse directions) are set as a free
    parameter while fixing all other parameters. The Mk model estimation
    is run and transition rate is recorded. The collected transition
    rates are clustered using `Mclust`. Each cluster is assigned a
    separate parameter.

Each of these model parameterizations; equal, symmetric and the two
custom models are tested and evaluated in new subtype pipeline. The
parameterization that has highest accuracy (based on a leave-one-out
cross-validation analysis) is selected. In the case of ties, the model
with the fewest parameters is given precidence (the symmetric model is
not tested when the number of subtypes is over 10)

### Evaluating new subtype schemes

#### Check 1. Clades with relatively small inter-patristic distance have the same subtype

In large cases, clades in the phylogenetic tree with distinct subtypes
may indicate a subtype that is not correlated with the phylogeny. In
limited cases, it might indicate annotation errors. Phylotyper computes
an inter-patritristic distance threshold for a given subtype scheme. The
threshold is equivalent to 0.4 probability that inter-patristic
distances equal to or greater then threshold belong to the same-subtype
distribution (basically collects all inter-patristic distances for nodes
with the same subtype and then uses R's Mclust to model the
distribution). Subtrees with a max inter-patristic distance less than
this threshold that have distinct subtypes are flagged. You will be
notified if this situation is detected. An updated subtype input file
called `phylotyper_proposed_subtypes.csv` will be generated with
proposed corrections and written to the output directory.

#### Check 2. The predictive performance is above a minimum threshold

All new subtype schemes are subject to a leave-one-out cross-validation
analysis. You will be notified if the F1-score (a equally weighted
average of precision and recall) is below 0.9.

CONTACT
=======

Matt Whiteside matthew.whiteside@phac-aspc.gc.ca

TODO
====
