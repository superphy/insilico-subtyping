#!/usr/bin/env python

"""Merge multiple fasta/subtype files into single non-redundant set


Example:
        $ python combine_phylotyper_inputs.py output_fasta_file output_subtype_file input_fasta_file1 input_subtype_file1 [input_fasta_file2 input_subtype_file2]

"""

import argparse
import logging
import os, sys, inspect
treefolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"..","..","tree")))
if treefolder not in sys.path:
    sys.path.insert(0, treefolder)
from seq import SeqDict

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mdwhitesi@gmail.com"


logging.basicConfig(
    filename='combine_phylotyper_inputs.log',
    level=logging.DEBUG,
    format='%(asctime)s %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
    filemode='w')
logger = logging.getLogger(__name__)


def filesets(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]



if __name__ == "__main__":
    """Merge files

    """

    # Parse command-line args
    parser = argparse.ArgumentParser(description='Combine phylotyper inputs')
    parser.add_argument('fastafile', action="store")
    parser.add_argument('subtypefile',action="store")
    parser.add_argument('inputs', default=[], nargs=argparse.REMAINDER)
    
    args = parser.parse_args()

    if len(args.inputs) == 0:
        msg = 'Missing input file argument.'
        raise Exception(msg)

    if len(args.inputs) % 2 != 0:
        msg = 'Input file arguments must be in sets of 2: fastafile, subtypefile.'
        raise Exception(msg)

    logger.debug("Collapsing identical sequences")
    seqdict = SeqDict(format_name=False)
    
    files = filesets(args.inputs, 2)
    for fs in files:
        fastafile = fs[0]
        subtypefile = fs[1]

        logger.debug("Processing files: {}, {}".format(fastafile, subtypefile))
        seqdict.build(fastafile, subtypefile)

    logger.debug("Identified {} non-redundant sequence/subtype pairs in input files".format(seqdict.num))

    # Output unique set
    seqdict.write(args.fastafile, args.subtypefile)

    subtype_counts = seqdict.subtype_occurences()
    logger.debug('Subtypes encountered:\n{}\n'.format(str(subtype_counts)))




   