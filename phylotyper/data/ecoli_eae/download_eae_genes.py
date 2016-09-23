#!/usr/bin/env python

"""Download Intimin (eae) gene sequences from Genbank


Example:
        $ python download_eae_genes.py .

"""

import argparse
import csv
import logging
import os
import re
from collections import Counter

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
    parser.add_argument('acc_file', action="store")
    
    args = parser.parse_args()

    # Initialize gene filter
    seq_tests = [lambda x: len(x) > 2500]
    gfilter = GeneFilter(sequence_tests=seq_tests)

    # Initialize Subtype parser
    subtype_names = '((?:alpha|beta|gamma|delta|epsilon|zeta|eta|theta|jota|iota|kappa|lambda|mu|nu|xi|omicron|pi|rho|sigma)(?:[\-_\s]?\d)?)'
    pattern1 = "(?:intimin|eae)[-_\s]%s" % subtype_names
    pattern2 = "%s[-_\s](?:intimin|eae)" % subtype_names
    pattern3 = "intimin[-_\s]type[-_\s]%s" % subtype_names

    parser = SubtypeParser([re.compile(pattern1),re.compile(pattern2),re.compile(pattern3)])


    # Initialize Download object
    dutil = DownloadUtils(args.output_directory, 'Escherichia coli', ['eae','Intimin'], parser, gfilter)

    # Perform Downloads
    dutil.download_pubmed_genes(pmidfile=args.pmid_file)

    dutil.download_by_accession(args.acc_file)

    # Parse genbank files for known intimin types
    dutil.parse()

    # Fixes, some subtypes are mislabelled
    replacements = [('AF530553.1_allele1', 'iota-2'), ('AF449420.1_allele1', 'theta-2'), 
    ('AF449419.1_allele1', 'theta-2'),('AF449415.1_allele1','theta-2'),('AF449414.1_allele1', 'theta-2'),
    ('AB334564.1_allele1','theta-2'),('AB334563.1_allele1', 'theta-2'),('AJ308551.1_allele1','iota-1'),
    ('AF449418.1_allele1', 'theta-2'),('AJ875027.1_allele1','kappa-1'),('AF253561.1_allele1','theta-2')
    ]
    count = Counter()
    subtypefile = os.path.join(args.output_directory, 'eae.txt')
    newsubtypefile = os.path.join(args.output_directory, 'tmp.txt')
    with open(newsubtypefile, 'w') as w:
        for x in csv.reader(open(subtypefile,'r'), delimiter='\t'):
            st = x[1]
            al = x[0]
            for k, v in replacements:
                if k in al:
                    st = v

            # append digit if one is missing
            st = re.sub(r'(\D)$',r'\1-1',st)
            w.write('\t'.join((x[0],st)) + '\n')
            count[st] += 1

    os.rename(newsubtypefile, subtypefile)

    print(count)



   