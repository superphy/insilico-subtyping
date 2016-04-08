#!/usr/bin/env python

"""Phylotyper program

Script for running various phylotyper functions



Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

.. _Google Python Style Guide:
   http://google.github.io/styleguide/pyguide.html
"""

import logging
import argparse
from config import PhylotyperOptions
from tree.fasttree import FastTreeWrapper


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mwhiteside@canada.ca"




if __name__ == "__main__":
	"""Run phylotyper function

    """
    
	logging.basicConfig(level=logging.DEBUG)

    # Parse command-line args
    # Phylotyper functions are broken up into commands
    # Each command has its own options and subparser
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(help='commands')

	# Alignment command
	aln_parser = subparsers.add_parser('aln', help='Build sequence alignment')
	aln_parser.add_argument('config', action='store', help='Phylotyper config options file')
	aln_parser.add_argument('input', action='store', help='Fasta input')
	aln_parser.add_argument('output', action='store', help='Alignment output')
	aln_parser.add_argument('--nt', action='store_true', help='Nucleotide sequences')

	aln_parser.set_defaults(which='aln')

	# Tree command
	tree_parser = subparsers.add_parser('tree', help='Build phylogenetic tree')
	tree_parser.add_argument('config', action='store', help='Phylotyper config options file')
	tree_parser.add_argument('input', action='store', help='Alignment input')
	tree_parser.add_argument('output', action='store', help='Phylogenetic tree output')
	tree_parser.add_argument('--nt', action='store_true', help='Nucleotide sequences')
	tree_parser.set_defaults(which='tree')

	options = parser.parse_args()

	config = PhylotyperOptions(options.config)

	if options.which == 'aln':
		pass
	elif options.which == 'tree':
		# Build tree
		
		tree = FastTreeWrapper(config)
		tree.build(options.input, options.output)

	else:
		raise Exception("Unrecognized command")

  