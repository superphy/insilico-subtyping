#!/usr/bin/env python

"""Concatenate fasta records in multiple fasta files

"""

import argparse
import logging
from Bio import SeqIO
from collections import defaultdict


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


logger = None

if __name__ == "__main__":
    """Concatenate fasta entries in multiple fasta files

    """

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('data.stx.translate')

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument('output', help='Output fasta file')
    parser.add_argument('input', nargs='*', help='Input fasta files')
    
    options = parser.parse_args()
    
    # Load
    seqs = defaultdict(list)
    
    for f in options.input:
        fasta = SeqIO.parse(f, 'fasta')
        for rec in fasta:
            seqs[rec.description].append(str(rec.seq))
    with open(options.output, 'w') as outfh:       
        for h,sarr in seqs.iteritems():
            outfh.write('>{}\n{}\n'.format(h,''.join(sarr)))



