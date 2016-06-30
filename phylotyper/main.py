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
from tree.seqaligner import SeqAligner


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"



def align_sequences(input, output, config):
	"""Build MSA

	Args:
	  input (str): Fasta file
	  output (str): Output file for MSA
	  config (obj): PhylotyperOptions object

    """

	aln = SeqAligner(config)
	aln.align(input, output)
	


def build_tree(input, output, nt, fast, config):
	"""Build phylogenetic tree

	Args:
	  input (str): Fasta file
	  output (str): Output file for newick tree
	  nt (bool): True when nucleotide sequences
	  fast (bool): True to prioritize speed over accuracy
	  config (obj): PhylotyperOptions object

    """

	tree = FastTreeWrapper(config)
	tree.build(input, output, nt, fast)
	

def predict_subtypes(inputs):
	pass


def subtype_pipeline(inputs):
	"""Run phylotyper pipeline

    """

    # Validate inputs


    # Align

    # Compute tree

    # Predict subtypes
	pass

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
	aln_parser.add_alignment('add', action='store', help='Existing Alignment file')
	aln_parser.add_alignment('reference', action='store', help='Phylotyper reference alignment file')
	aln_parser.set_defaults(which='aln')

	# Tree command
	tree_parser = subparsers.add_parser('tree', help='Build phylogenetic tree')
	tree_parser.add_argument('config', action='store', help='Phylotyper config options file')
	tree_parser.add_argument('input', action='store', help='Alignment input')
	tree_parser.add_argument('output', action='store', help='Phylogenetic tree output')
	tree_parser.add_argument('--nt', action='store_true', help='Nucleotide sequences')
	tree_parser.set_defaults(which='tree')

	# Subtype command
	subtype_parser = subparsers.add_parser('subtype', help='Predict subtype')
	subtype_parser.add_argument('config', action='store', help='Phylotyper config options file')
	subtype_parser.add_argument('sequences', action='store', help='Reference sequences for tree')
	subtype_parser.add_argument('subtype', action='store', help='Reference subtypes')
	subtype_parser.add_argument('input', action='store', help='Fasta input for unknowns')
	subtype_parser.add_argument('output', action='store', help='Subtype predictions')
	subtype_parser.add_argument('--nt', action='store_true', help='Nucleotide sequences')
	subtype_parser.set_defaults(which='subtype')

	options = parser.parse_args()

	config = PhylotyperOptions(options.config)

	if options.which == 'aln':
		# Build alignment

		# Run
		align_sequences(options.input, options.output, config)
		
	elif options.which == 'tree':
		# Build tree

		# Nucleotide sequences
		nt = options.nt

		# Fast mode of tree calculation
		fast = False
		
		# Run
		build_tree(options.input, options.output, nt, fast, config)

	elif options.which == 'subtype':
		# Compute subtype

		# Nucleotide sequences
		nt = options.nt

		# Fast mode of tree calculation
		fast = False
		
		

	else:
		raise Exception("Unrecognized command")

  

