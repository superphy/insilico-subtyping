#!/usr/bin/env python

"""Download Intimin (eae) gene sequences from Genbank



Example:
        $ python download_eae_genes.py

"""

import argparse
import logging
import re

from utils import DownloadUtils, SubtypeParser

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
    # parser.add_argument('--dbuser', action="store")
    # parser.add_argument('--dbhost', action="store")
    # parser.add_argument('--dbport', action="store")
    # parser.add_argument('--dbpass', action="store")
    # parser.add_argument('--blastdbcmd', action="store", default='blastdbcmd')
    # parser.add_argument('--diamondcmd', action="store", default='diamond')
    # parser.add_argument('--db_dir', action="store")
    # parser.add_argument('--skip_gi', action="store_true")
    # parser.add_argument('--skip_ftp', action="store_true")
    # parser.add_argument('--update', action="store_true")
    
    args = parser.parse_args()

    # Initialize Subtype parser
    subtype_names = '((?:alpha|beta|gamma|delta|epsilon|zeta|eta|theta|jota|kappa|lambda|mu|nu|xi|omicron|pi|rho|sigma)(?:[\-_\s]?\d)?)'
    pattern1 = "(?:intimin|eae)[-_\s]%s" % subtype_names
    pattern2 = "%s[-_\s](?:intimin|eae)" % subtype_names

    parser = SubtypeParser([re.compile(pattern1),re.compile(pattern2)])


    # Initialize Download object
    dutil = DownloadUtils(args.output_directory, 'Escherichia coli', ['eae','Intimin'], parser)

    # Perform Download
    #dutil.download()

    # Parse genbank files for known intimin types
    dutil.parse()

    # Trim output alignment



   