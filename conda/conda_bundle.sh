conda config --add channels compbiocore
conda config --add channels bioconda
conda config --add channels defaults

conda metapackage phylotyper-bundle 0.1.0 --dependencies r r-ggplot2 r-optparse r-fitdistrplus r-mclust r-ape r-rocr r-rcolorbrewer r-devtools r-phangorn r-phytools r-igraph bioconductor-biostrings mafft blast fasttree trimal 
