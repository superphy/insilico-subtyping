#!/usr/bin/env python

"""Download fliC gene sequences from Genbank


Example:
        $ python download_flic_genes.py .

"""

import argparse
import logging
import re

from utils import DownloadUtils, SubtypeParser, GeneFilter

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mdwhitesi@gmail.com"


logging.basicConfig(
    filename='download_flic_genes.log',
    level=logging.DEBUG,
    format='%(asctime)s %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
    filemode='w')

if __name__ == "__main__":
    """Run pipeline downloading sequences

    """

    # Parse command-line args
    parser = argparse.ArgumentParser(description='Download and store NCBI genes sequences')
    parser.add_argument('output_directory', action="store")
    
    args = parser.parse_args()

    # Initialize gene filter
    seq_tests = [lambda x: len(x) > 1400]
    gfilter = GeneFilter(sequence_tests=seq_tests)

    # Initialize serotype parser
    hpattern = r"[o0](?:\d+|nt|r|untypeable|\?|\-)[a-z]?\:(h\d+[a-z]?)(?:\:|\b)"
    parser = SubtypeParser([re.compile(hpattern, flags=re.IGNORECASE)], source_fields=['organism','strain','serotype'],
        annotation_fields=['source','serotype','organism'])

    # Initialize Download object
    dutil = DownloadUtils(args.output_directory, 'Escherichia coli', ['fliC'], parser, gfilter, 
        gene_regex=[re.compile(r'\bfliC\b', flags=re.IGNORECASE), re.compile(r'^flagellin$',flags=re.IGNORECASE)])

    # Perform Downloads
    #dutil.download_genes()

    # Parse genbank files for H-types
    dutil.parse_subtype()



   