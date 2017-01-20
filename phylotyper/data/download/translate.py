#!/usr/bin/env python

"""Translate sequences

"""

import argparse
import logging
import re
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC



__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


logger = None

if __name__ == "__main__":
    """Translate sequences

    """

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('data.stx.translate')

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Input fasta file')
    parser.add_argument('output', help='Output fasta file')
    
    options = parser.parse_args()

    # Translate
    with open(options.output, 'w') as outfh:
        fasta = SeqIO.parse(options.input, 'fasta')
        for rec in fasta:
            sobj = Seq(str(rec.seq).replace('-',''), IUPAC.ambiguous_dna)
            pobj = sobj.translate(table=11, stop_symbol='')

            outfh.write('>{}\n{}\n'.format(rec.description, str(pobj)))



