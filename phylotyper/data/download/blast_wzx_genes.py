#!/usr/bin/env python

"""Blast wzx gene sequences against Genbank


Example:
        $ python blast_wzx_genes.py .

"""

import argparse
import csv
import logging
import os
import re
from Bio import SeqIO


from blast import Blast
from utils import DownloadUtils, SubtypeParser


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mdwhitesi@gmail.com"


logging.basicConfig(
    filename='blast_wzx_genes.log',
    level=logging.DEBUG,
    format='%(asctime)s %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
    filemode='w')
logger = logging.getLogger(__name__)

if __name__ == "__main__":
    """Run

    """

    # Parse command-line args
    parser = argparse.ArgumentParser(description='Download and store NCBI blast results')
    parser.add_argument('fasta_file', action="store")
    parser.add_argument('output_directory', action="store")

    args = parser.parse_args()
    tmp_fasta_file = os.path.join(args.output_directory, 'tmp.fasta')

    blast = Blast("Escherichia coli")
    blast.run(args.fasta_file, tmp_fasta_file)

    # Initialize serotype parser
    opattern = r"([o0](?:\d+|nt|r|untypeable|\?|\-)[a-z]?)\:(?:h\d+[a-z]?)(?:\:|\b)"
    parser = SubtypeParser([re.compile(opattern, flags=re.IGNORECASE)], source_fields=['organism','strain','serotype','serovar'],
        annotation_fields=['source','serotype','organism','serovar'])

    # Initialize Download object
    dutil = DownloadUtils(args.output_directory, 'Escherichia coli', ['wzx'], parser, None)

    # Perform Downloads
    dutil.download_by_accession(tmp_fasta_file, fasta_format=True)

    # Parse genbank files for H-types
    dutil.parse_subtype()

    # Generate final output
    invalid = set([])
    subtypes = {}
    for row in csv.reader(open(dutil.subtypefile,'r'),delimiter='\t'):
        name = row[0]
        subt = row[1]

        if not name in subtypes:
            # Never seen before
            subtypes[name] = subt

        elif subtypes[name] == subt:
            # Non-conflicting duplicates
            logger.info("Duplicate instances of {}, subtype:{} in subtype file: {}".format(name,subt,dutil.subtypefile))

        else:
            # Conflict
            logger.warning("Duplicate instances of {} with conflicting subtypes: {} vs {} in subtype file: {}. Removing entry.".format(name,subt,
                    subtypes[name],dutil.subtypefile))
            invalid.add(name)



    seqs = {}
    invalid2 = set([])
    with open(tmp_fasta_file, 'r') as infh:
        for seqrecord in SeqIO.parse(infh, 'fasta'):
            
            name = seqrecord.id
            seq = seqrecord.seq

            if not name in invalid:

                if name in subtypes:
                    # Has matching subtype
                    
                    if name in seqs:
                        # Duplicate in fasta file
                        if seqs[name] == seq:
                             logger.info("Duplicate instances of {} in fasta file: {}".format(name, dutil.fastafile))
                        else:
                            logger.warning("Duplicate instances of {} with conflicting sequences in fasta file: {}".format(name, dutil.fastafile))
                            invalid2.add(name)

                    else:
                        # novel sequence with subtype
                        seqs[name] = seq
                        
                else:
                    # No matching subtype
                    logger.info("No subtype for sequence {} in fasta file: {}".format(name, dutil.fastafile))

            

    with open(dutil.fastafile, 'w') as outfh, open(dutil.subtypefile, 'w') as subtfh:
        for name in seqs:

            outfh.write('>{}\n{}\n'.format(name, seqs[name]))
            subtfh.write('{}\t{}\n'.format(name, subtypes[name]))
   



   