#!/usr/bin/env python

"""Syncing subtype tab-delim file and fasta files.

"""

import argparse
import csv
import logging
import re
from Bio import SeqIO
from datetime import datetime
from shutil import copyfile


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


logger = None


def backup_file(file):
    """Create timestamp'd backup of a file

    Args:
        file (str): filepath

    Returns:
        backupfile(str)

    """

    current_time = datetime.now()
    time_stamp =  current_time.strftime("%b-%d-%y-%H.%M.%S")
    backupfile = file +'.bkp_'+ time_stamp
    copyfile(file, backupfile)

    return(backupfile)


if __name__ == "__main__":
    """Filter subtypefile and fastafiles for accessions matching ids

        A backup (.bkp) of each file will be made and then file will be altered in-place.

    """

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('data.download.filter')

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument('id', action='store', help='ID list')
    parser.add_argument('subtype', action='store', help='Subtypes for input sequences')
    parser.add_argument('fasta', nargs='+', help='Fasta input(s)')
    
    options = parser.parse_args()

    # Load IDs
    with open(options.id, 'r') as fh:
        content = fh.readlines()
    ids = [x.strip() for x in content]

    # Parse IDs?
    if re.match(r'^\d+\-.+__', ids[0]):
        ids2=list()
        for i in ids:
            m = re.search(r'^\d+\-(.+)__', i)
            ids2.append(m.group(1))
        ids = ids2
    lookup = set(ids)
    
    # Backup files
    subtypefile = backup_file(options.subtype)
    fastafiles = []
    for f in options.fasta:
        fastafiles.append(backup_file(f))

     # Filter subtypefile
    with open(options.subtype, 'w') as outfh:
        for row in csv.reader(open(subtypefile,'r'),delimiter='\t'):
            name = row[0]
            subt = row[1]

            if name in lookup:
                outfh.write("{}\t{}\n".format(name, subt))

    # Filter fastafiles
    for i,f in enumerate(options.fasta):
        bf = fastafiles[i]

        with open(f, 'w') as outfh:
            fasta = SeqIO.parse(bf, 'fasta')
            for rec in fasta:
                hasallele = re.search('^(.+)\|allele\d+$', rec.id)
                if hasallele:
                    acc = hasallele.group(1)
                else:
                    acc = rec.id

                if acc in lookup:
                    outfh.write('>{}\n{}\n'.format(rec.description, rec.seq))







