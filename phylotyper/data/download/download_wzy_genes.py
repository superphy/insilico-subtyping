#!/usr/bin/env python

"""Download Wzy gene sequences from Genbank


Example:
        $ python download_wzy_genes.py .

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
    filename='download_wzy_genes.log',
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
    parser.add_argument('acc_file', action="store")
    
    args = parser.parse_args()

    # Initialize gene filter
    seq_tests = [lambda x: len(x) > 800]
    gfilter = GeneFilter(sequence_tests=seq_tests)

    # Initialize Subtype parser
    opattern = r"(?:\b|serogroup\:)([o0]x?\d+(?:[a-z]{1,2})?)(?:\b|\:)"
    parser = SubtypeParser([re.compile(opattern, flags=re.IGNORECASE)], 
        source_fields=['organism','strain','serotype','serovar','note'],
        annotation_fields=['source','serotype','organism','serovar'])

    # Initialize Download object
    dutil = DownloadUtils(args.output_directory, 'Escherichia coli', ['wzy','O-antigene polymerase','O-antigen polymerase'], parser, gfilter)

    # Download
    dutil.download_by_accession(args.acc_file)

    # Parse genbank files for known intimin types
    dutil.parse()




   