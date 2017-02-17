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

def translate(sobj, l):
    # Translate genome hits to amino acid sequences
    # Discards failed translations 

    global logger

    try:
        pobj = sobj.translate(table=11, stop_symbol='*')
    except Exception as e:
        logger.warning('Translation failed for {}: {}'.format(l,e))

    nstops = pobj.count('*')

    if nstops > 1:
        try:
            sobj1 = sobj.reverse_complement()
            pobj1 = sobj1.translate(table=11, stop_symbol='*')

            nstops1 = pobj1.count('*')
            if nstops1 > 1:
                logger.warning('Frameshift error for {} ({} & {} stop codons in forward and reverse)'.format(l, nstops, nstops1))
                if nstops1 < nstops:
                    sobj = sobj1
                    pobj = pobj1
            else:
                sobj = sobj1
                pobj = pobj1

        except Exception as e:
            logger.warning('Translation failed for {}: {}'.format(l,e))

    elif not re.search(r'^M', str(pobj)):
        try:
            sobj1 = sobj.reverse_complement()
            pobj1 = sobj1.translate(table=11, stop_symbol='*')

            if re.search(r'^M', str(pobj1)):
                sobj = sobj1
                pobj = pobj1

        except Exception as e:
            logger.warning('Translation failed for {}: {}'.format(l,e))

    elif not re.search(r'\*$', str(pobj)):
        try:
            sobj1 = sobj.reverse_complement()
            pobj1 = sobj1.translate(table=11, stop_symbol='*')

            if re.search(r'\*$', str(pobj)):
                sobj = sobj1
                pobj = pobj1

        except Exception as e:
            logger.warning('Translation failed for {}: {}'.format(l,e))

    return pobj


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
            pobj = translate(sobj, rec.id)

            outfh.write('>{}\n{}\n'.format(rec.description, str(pobj)))


