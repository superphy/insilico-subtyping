#!/usr/bin/env python

"""Download Intimin (eae) gene sequences from Genbank



Example:
        $ python download_eae_genes.py

"""

import argparse
import logging

from utils import DownloadUtils

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

    # Initialize Download object

    dutil = DownloadUtils(args.output_directory, 'Escherichia coli', ['eae','Intimin'])

    dutil.download()

    # Perform Download



   