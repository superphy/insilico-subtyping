#!/usr/bin/env python

"""Merge multiple fasta/subtype files into single non-redundant set


Example:
        $ python combine_phylotyper_inputs.py output_fasta_file output_subtype_file input_fasta_file1 input_subtype_file1 [input_fasta_file2 input_subtype_file2]

"""

import argparse
import logging
import csv
from Bio import SeqIO
from collections import Counter


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


def collate(uniq, invalid, fastafile, subtypefile):
    """Record non-conflicting, non-redundant sequence subtype pairs

    Results are stored in uniq dictionary

        Args:
            pmids (list): List of ints
            outputfile (str): Filepath to write results to


    """

    
    pending = {}
    
    for row in csv.reader(open(subtypefile,'r'),delimiter='\t'):
        name = row[0].upper()
        subt = row[1].lower()

        if not name in pending:
            # Never seen before
            pending[name] = subt

        elif pending[name] == subt:
            # Non-conflicting duplicates
            logger.info("Duplicate instances of {}, subtype:{} in subtype file: {}".format(name,subt,subtypefile))

        else:
            # Conflict
            logger.warning("Duplicate instances of {} with conflicting subtypes: {} vs {} in subtype file: {}. Removing entry.".format(name,subt,
                    pending[name],subtypefile))
            invalid.add(name)


    with open(fastafile, 'r') as fh:
        for seqrecord in SeqIO.parse(fh, 'fasta'):
            
            name = seqrecord.id.upper()
            seq = seqrecord.seq.upper()

            if not name in invalid:

                if name in pending:
                    # Has matching subtype
                    
                    if name in uniq:
                        # Duplicate in fasta file
                        if uniq[name][0] == seq and uniq[name][1] == pending[name]:
                             logger.info("Duplicate instance of {} in fasta file: {}".format(name, fastafile))
                        else:
                            logger.warning("Duplicate instances of {} with conflicting sequence/subtype in fasta file: {}".format(name, fastafile))
                            invalid.add(name)
                            uniq.pop(name, None)

                    else:
                        # novel sequence with subtype
                        uniq[name] = (seq,pending[name])
                        

                else:
                    # No matching subtype
                    logger.info("No subtype for sequence {} in fasta file: {}".format(name, fastafile))

    

    return None


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


    files = filesets(args.inputs, 2)
    uniq = {}
    invalid = set([])
    for fs in files:
        fastafile = fs[0]
        subtypefile = fs[1]

        logger.debug("Processing files: {}, {}".format(fastafile, subtypefile))
        collate(uniq, invalid, fastafile, subtypefile)

    logger.debug("Identified {} non-redundant sequence/subtype pairs in input files".format(len(uniq)))

    subtype_counts = Counter()
    with open(args.fastafile, 'w') as fastah, open(args.subtypefile, 'w') as subth:
        for name in uniq:
            seq = uniq[name][0]
            subt = uniq[name][1]

            fastah.write('>{}\n{}\n'.format(name, seq))
            subth.write('{}\t{}\n'.format(name, subt))

            subtype_counts[subt] += 1


    logger.debug('Subtypes encountered:\n{}\n'.format(str(subtype_counts)))




   