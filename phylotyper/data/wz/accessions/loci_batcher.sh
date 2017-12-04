#!/bin/bash

######################################################
## Name: loci_batcher.sh
## Desc: Helper script to run phylotyper loci mode 
##       on batches of genomes in a directory
## Author: Matt Whiteside
## Date: Dec 1, 2017
##
## Usage: loci_batcher.sh inputdir outputdir
######################################################

INPUTDIR=$1
OUTPUTDIR=$2

files=($INPUTDIR/*.fasta)
g=2

for((i=0; i < ${#files[@]}; i+=g))
do
  part=( "${files[@]:i:g}" )
  python -m phylotyper loci wz $OUTPUTDIR "$part[*]"
  cat $OUTPUTDIR/subtype_gene_predictions.tsv >> $OUTPUTDIR/batch_subtype_gene_predictions.tsv
done

