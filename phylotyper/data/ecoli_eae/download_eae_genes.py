#!/usr/bin/env python

"""Download Intimin (eae) gene sequences from Genbank


Example:
        $ python download_eae_genes.py .

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
    filename='download_eae_genes.log',
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
    parser.add_argument('pmid_file', action="store")
    
    args = parser.parse_args()

    # Initialize gene filter
    seq_tests = [lambda x: len(x) > 2500]
    gfilter = GeneFilter(sequence_tests=seq_tests)

    # Initialize Subtype parser
    subtype_names = '((?:alpha|beta|gamma|delta|epsilon|zeta|eta|theta|jota|kappa|lambda|mu|nu|xi|omicron|pi|rho|sigma)(?:[\-_\s]?\d)?)'
    pattern1 = "(?:intimin|eae)[-_\s]%s" % subtype_names
    pattern2 = "%s[-_\s](?:intimin|eae)" % subtype_names

    parser = SubtypeParser([re.compile(pattern1),re.compile(pattern2)])


    # Initialize Download object
    dutil = DownloadUtils(args.output_directory, 'Escherichia coli', ['eae','Intimin'], parser, gfilter)

    # Perform Download
    dutil.download_pubmed_genes(pmidfile=args.pmid_file)

    # Parse genbank files for known intimin types
    dutil.parse()


   