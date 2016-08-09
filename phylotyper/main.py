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
import os
from Bio import SeqIO
from config import PhylotyperOptions
from builtin_subtypes import SubtypeConfig
from phylotyper import Phylotyper
from tree.fasttree import FastTreeWrapper
from tree.seqaligner import SeqAligner


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"



def align_all_sequences(input, output, config):
	"""Build MSA

	Args:
	  input (str): Fasta file
	  output (str): Output file for MSA
	  config (obj): PhylotyperOptions object

    """

	aln = SeqAligner(config)
	aln.align(input, output)


def align_new_sequences(input, alignment, output, config):
	"""Add new sequences to existing MSA using
	profile alignment

	Args:
	  input (str): Fasta file
	  alignment (str): Aligned fasta file
	  output (str): Output file for MSA
	  config (obj): PhylotyperOptions object

    """

	aln = SeqAligner(config)
	aln.add(input, alignment, output)
	


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
	

def predict_subtypes(options, config):
	"""Phylotyper subtype prediction


	"""

	# Define files
	subtfile = os.path.abspath(options['subtype'])
	treefile = os.path.abspath(os.path.join(options['output_directory'], 'combined.tree'))

	pt = Phylotyper()

	pt.subtype(treefile, subtfile, config)



def subtype_pipeline(options, config):
	"""Run phylotyper pipeline

    """

    # Define files
	alnfile = os.path.join(options['output_directory'], 'combined.aln')
	oldalnfile = options['alignment']
	treefile = os.path.join(options['output_directory'], 'combined.tree')

	# Rename sequences with unique ids
	uniquify_sequences(options)

    # Align
	align_new_sequences(options['input'], oldalnfile, alnfile, config)

	# Compute tree
	nt = options['seq'] == 'nt'
	build_tree(alnfile, treefile, nt, options['fast'], config)

    # Predict subtypes
	predict_subtypes(options, config)


def uniquify_sequences(options):
	"""Create temporary sequence names that are unique

    """

    # Define files
	output_file = os.path.join(options['output_directory'], 'input.fasta')
	mapping_file = os.path.join(options['output_directory'], 'mapping.txt')
	input_file = options['input']
	
	i = 1
	pre = 'pt_'
	fasta_sequences = SeqIO.parse(open(input_file),'fasta')
	with open(output_file, 'w') as out, open(mapping_file, 'w') as mapping:
	    for fasta in fasta_sequences:
	        name, sequence = fasta.description, str(fasta.seq)
	        newname = '%s%i' % (pre, i)
	        i += 1
	        mapping.write('%s\t%s\n' % (name, newname))
	        out.write('>%s\n%s\n' % (newname, sequence))

	options['input'] = output_file
	options['user_input'] = input_file



if __name__ == "__main__":
	"""Run phylotyper function

    """

   	subtype_config_file = 'builtin_subtypes.yaml'
   
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
	aln_parser.add_argument('add', action='store', help='Existing Alignment file')
	aln_parser.add_argument('reference', action='store', help='Phylotyper reference alignment file')
	aln_parser.set_defaults(which='aln')

	# Tree command
	tree_parser = subparsers.add_parser('tree', help='Build phylogenetic tree')
	tree_parser.add_argument('config', action='store', help='Phylotyper config options file')
	tree_parser.add_argument('input', action='store', help='Alignment input')
	tree_parser.add_argument('output', action='store', help='Phylogenetic tree output')
	tree_parser.add_argument('--nt', action='store_true', help='Nucleotide sequences')
	tree_parser.set_defaults(which='tree')

	# User-supplied subtype command
	subtype_parser = subparsers.add_parser('custom', help='Predict subtype for user-defined subtype scheme')
	subtype_parser.add_argument('config', action='store', help='Phylotyper config options file')
	subtype_parser.add_argument('sequences', action='store', help='Fasta sequences of aligned reference genes for tree')
	subtype_parser.add_argument('subtype', action='store', help='Reference gene subtypes')
	subtype_parser.add_argument('input', action='store', help='Fasta input for unknowns')
	subtype_parser.add_argument('output', action='store', help='Directory for Subtype predictions')
	subtype_parser.add_argument('--nt', action='store_true', help='Nucleotide sequences')
	subtype_parser.set_defaults(which='custom')

	# Builtin subtype command
	subtype_parser = subparsers.add_parser('subtype', help='Predict subtype for scheme provided in phylotyper')
	subtype_parser.add_argument('config', action='store', help='Phylotyper config options file')
	subtype_parser.add_argument('gene', action='store', help='Subtype gene name')
	subtype_parser.add_argument('input', action='store', help='Fasta input for unknowns')
	subtype_parser.add_argument('output', action='store', help='Directory for subtype predictions')
	subtype_parser.add_argument('--nt', action='store_true', help='Nucleotide sequences')
	subtype_parser.set_defaults(which='subtype')


	options = parser.parse_args()

	config = PhylotyperOptions(options.config)

	if options.which == 'aln':
		# Build alignment

		# Run
		align_all_sequences(options.input, options.output, config)
		
	elif options.which == 'tree':
		# Build tree

		# Nucleotide sequences
		nt = options.nt

		# Fast mode of tree calculation
		fast = False
		
		# Run
		build_tree(options.input, options.output, nt, fast, config)

	elif options.which == 'subtype':
		# Compute subtype for builtin scheme

		# Check arguments

		# Nucleotide sequences
		nt = options.nt

		# Fast mode of tree calculation
		fast = False

		# Check input file exists
		if not os.path.isfile(options.input):
			msg = 'Invalid/missing input file argument.'
			raise Exception(msg)

		# Check output directory exists, if not create it if possible
		if not os.path.exists(options.output):
			os.makedirs(options.output)

		# Load requested subtype data files
		stConfig = SubtypeConfig(subtype_config_file)
		scheme = options.gene

		subtype_options = stConfig.get_subtype_config(scheme)
		subtype_options['input'] = options.input
		subtype_options['output_directory'] = options.output
		subtype_options['fast'] = False

		if options.nt and (subtype_options['seq'] != 'nt'):
			msg = 'Sequence type of input does not match Phylotyper gene sequences for %s' % (scheme)
			raise Exception(msg)

		# Run pipeline
		subtype_pipeline(subtype_options, config)






		
		

	else:
		raise Exception("Unrecognized command")

  

