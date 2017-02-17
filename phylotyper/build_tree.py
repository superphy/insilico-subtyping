#!/usr/bin/env python

"""Build tree from alignment

Script for building phylogenetic tree from alignment

"""


import argparse
import logging
import os
import tempfile

from config import PhylotyperOptions
from tree.fasttree import FastTreeWrapper
from tree.seq import LociConcat

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


logger = None

if __name__ == "__main__":
    """Run loci search

    Parses command-line arguments and calls LociSearch()
    functions

    """

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('phylotyper.build_tree')

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', action='store', help='Phylotyper config options file')
    parser.add_argument('output', help='Newick output file')
    parser.add_argument('inputs', nargs='+', help='Alignment fasta files')
    parser.add_argument('--aa', action='store_true', help='Amino acid sequences')
   
    options = parser.parse_args()

    # Parse .ini config file
    # Location of config file is defined by ENV variable PHYLOTYPER_CONFIG or by --config (overrides previous)
    if options.config:
        config_file = options.config
    elif os.environ.get('PHYLOTYPER_CONFIG'):
        config_file = os.environ.get('PHYLOTYPER_CONFIG')
    else:
        msg = 'Missing config file argument.\nMust provide Phylotyper config file using' \
            ' enviroment variable PHYLOTYPER_CONFIG or command-line argument --config.'
        raise Exception(msg)
    config = PhylotyperOptions(config_file)

    ft = FastTreeWrapper(config)
    nt = not options.aa

    if len(options.inputs) > 1:
        with tempfile.NamedTemporaryFile() as tmpfh:

            concat = LociConcat()
            concat.safe_collapse(options.inputs, tmpfh.name)
   
            ft.build(tmpfh.name, options.output, nt=nt)

    else:
        ft.build(options.inputs[0], options.output, nt=nt)
        

   
  

