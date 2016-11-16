#!/usr/bin/env python

"""Loci Detection

Script for finding gene loci in input genomes

Examples:
    To run loci detection:

        $ python search.py genome.fasta --db phylotyper_blast_db [--fasta loci1.fasta [loci2.fasta ...]] 

"""


import argparse
import logging

from config import PhylotyperOptions
from genome.loci import LociSearch

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


logger = None

if __name__ == "__main__":
    """Run loci search

    Parses command-line arguments and calls LociSearch()
    functions

    """

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('phylotyper.genome.main')

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help='Phylotyper config options file')
    parser.add_argument('output', help='Output loci sequences')
    parser.add_argument('db', help='Blast database')
    parser.add_argument('genome', nargs='*', help='Input genome DNA sequence')
    parser.add_argument('--fasta', nargs='*', help='Input reference protein/gene sequences', )
    parser.add_argument('--aa', action='store_true', help='Amino acid sequence database')
    parser.add_argument('--list', help='Optionally specify file listing file paths to input genomes')
   
    options = parser.parse_args()

    # Parse .ini config file
    config = PhylotyperOptions(options.config)

    # Initialization
    seqtype = 'prot' if options.aa else 'nucl'
    detector = LociSearch(config, options.db, options.fasta, seqtype)

    # Search for instances of loci in genomes
    genomes = options.genome

    if options.list:
        with open(options.list) as f:
            genomes = [x.strip('\n') for x in f.readlines()]

    for genome in genomes:
        detector.search(genome, options.output, append=True)
        

   
  

