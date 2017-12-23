#!/usr/bin/env bash

######################################################
## batch_run.sh
##
## Run phylotyper in batches of 10, concatenating
## results for each batch into single output results
## file
##
## 
##
#####################################################

# Command-line args
gene=$1
output=$2
infile=$3

if [ -z "$gene" ]; then
	echo "Missing gene argument"
	exit 1
fi

if [ -z "$output" ] || [ ! -d "$output" ]; then
	echo "Missing/invalid output directory argument"
	exit 1
fi

if [ -z "$infile" ] || [ ! -f "$infile" ]; then
	echo "Missing/invalid input file argument"
	exit 1
fi

if [ -e ptbatch_failed.txt ]; then
  rm ptbatch_failed.txt
fi

# Batch result file
batchresfile="$output/batch_subtype_predictions.tsv"

# Load inputs
declare -a inputs
readarray -t inputs < $infile

# Iterate over inputs in chunks
firstheader=true
g=30
for((i=0; i < ${#inputs[@]}; i+=g))
do
  part=( "${inputs[@]:i:g}" )
  nparts="${#part[@]}"
  # Create tmp directory
  tmpdir=`mktemp -d /tmp/ptbatchXXXX`
  
  echo "Iteration: ${i}"
  # echo "Inputs in this batch: ${part[*]}"
  # echo "Output directory: ${tmpdir}"

  # Run phylotyper batch
  resfile="$tmpdir/subtype_predictions.tsv"
  cmd="python -m phylotyper genome $gene $tmpdir "${part[*]}""
  #echo $cmd
  eval $cmd

  # Check output
  success=true
  if [ ! -f "$resfile" ]; then
  	success=false
  else
  	lines=`wc -l $resfile | cut -f1 -d' '`
    if [ "$lines" -le "$nparts" ]; then
    	success=false
    fi
  fi
  
  if [ "$success" = true ]; then
  	# Save output
  	echo "SUCCESS"
  	if [ "$firstheader" = true ]; then
  	  cat $resfile > $batchresfile
  	  firstheader=false
  	else
  	  tail -n +2 $resfile >> $batchresfile
  	fi 
  else
  	# Failed
  	echo "FAILED"
  	printf "%s\n" ${part[@]} >> ptbatch_failed.txt
  fi

  #rm -rf $tmpdir
  echo "DONE"
  echo ""

done





