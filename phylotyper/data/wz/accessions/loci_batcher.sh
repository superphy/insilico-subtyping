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

files=($INPUTDIR/*.ffn)
g=100

for((i=0; i < ${#files[@]}; i+=g))
do
  part=( "${files[@]:i:g}" )
  python -m phylotyper loci wz $OUTPUTDIR "${part[@]}"
  cat $OUTPUTDIR/subtype_gene_predictions.tsv >> $OUTPUTDIR/all_subtype_gene_predictions.tsv
done


for file in $OUTPUTDIR/*_loci2.fasta; do mv "$file" "${file/_loci2.fasta/_wzx.fasta}"; done
for file in $OUTPUTDIR/*_loci1.fasta; do mv "$file" "${file/_loci1.fasta/_wzy.fasta}"; done

