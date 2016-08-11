#!/usr/bin/env python

"""Phylotyper program main

Script for running various phylotyper functions

Examples:
	To run subtyping routine on Stx1 genes:

        $ python main.py subtype ../phylotyper_example.ini ecoli_stx1 test/ecoli_stx1.ffn output/test/

"""


import argparse
import csv
import logging
import os
import pprint
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


logger = None

def align_all_sequences(input, output, config):
	"""Build MSA

	Args:
		input (str): Fasta file
		output (str): Output file for MSA
		config (obj): PhylotyperOptions object

    """

	logger.debug('Performing full alignment')

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

	logger.debug('Aligning genes to existing alignment')

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

	logger.debug('Building fasttree')

	tree = FastTreeWrapper(config)
	tree.build(input, output, nt, fast)
	

def predict_subtypes(options, config):
	"""Calls subtype method

	Wrapper around the Phylotyper.subtype method. Identifies
	predictions above significance cutoff. Writes results to 
	file.

	Args:
        options (dict): user defined settings from __main__
        config (obj): PhylotyperConfig with .ini file settings

	"""

	logger.debug('Running phylotyper')

	# Define files
	subtfile = options['subtype']
	treefile = options['tree_file']
	assfile = options['result_file']

	pt = Phylotyper(config)

	assignment_dict = pt.subtype(treefile, subtfile, options['output_directory'])

	# Load mapping
	rev_mapping = {}
	for x in csv.reader(open(options['mapping_file'],'r'), delimiter='\t'):
		rev_mapping[x[1]] = x[0]

	# Write assignments
	cutoff = float(config.get('phylotyper','prediction_threshold'))
	logger.debug('Using posterior probability cutoff: %f' % (cutoff))

	with open(assfile, 'w') as outf:
		outf.write('#genome\tsubtype\tprobability\tphylotyper_assignment\n')
		for genome_id,subtup in assignment_dict.items():

			pp = subtup[1]
			subt = ','.join(subtup[0])

			pred = 'undetermined'
			print pp
			if pp > cutoff:
				pred = subt
			
			outf.write('\t'.join((
				rev_mapping[genome_id],
				','.join(subt),
				str(pp),
				pred +'\n'
			)))

	None


def subtype_pipeline(options, config):
	"""Run phylotyper pipeline 

	Runs individual steps in phylotyper pipeline
	for internal subtype scheme
        
        Args:
            options (dict): user defined settings from __main__
            config (obj): PhylotyperConfig with .ini file settings

    """

    # Define files
	alnfile = os.path.join(options['output_directory'], 'combined.aln')
	oldalnfile = options['alignment']
	treefile = os.path.join(options['output_directory'], 'combined.tree')
	options['tree_file'] = treefile
	options['result_file'] = os.path.join(options['output_directory'], 'subtype_predictions.txt')

	# Rename sequences with unique ids, change input files
	uniquify_sequences(options)

	logger.info('Settings:\n%s' % (pprint.pformat(options)))
	logger.info('Config:\n%s' % (config.pformat()))

    # Align
	align_new_sequences(options['input'], oldalnfile, alnfile, config)

	# Compute tree
	nt = options['seq'] == 'nt'
	build_tree(alnfile, treefile, nt, options['fast'], config)

    # Predict subtypes & write to file
	predict_subtypes(options, config)


def uniquify_sequences(options):
	"""Create temporary sequence names that are unique

	Autogenerate temporary name placeholders and rewrite
	inputs with these names.

	Sets filename key-values in options dict for new temp input file
	and mapping file.
        
        Args:
            options (dict): user defined settings from __main__

    """

	logger.debug('Renaming sequences')

    # Define files
	output_file = os.path.join(options['output_directory'], 'input.fasta')
	mapping_file = os.path.join(options['output_directory'], 'mapping.txt')
	input_file = options['input']

	i = 1
	pre = 'pt_'
	fasta_sequences = SeqIO.parse(open(input_file),'fasta')
	with open(output_file, 'w') as outf, open(mapping_file, 'w') as mapf:
	    for fasta in fasta_sequences:
	        name, sequence = fasta.description, str(fasta.seq)
	        newname = '%s%i' % (pre, i)
	        i += 1
	        mapf.write('%s\t%s\n' % (name, newname))
	        outf.write('>%s\n%s\n' % (newname, sequence))

	options['input'] = output_file
	options['user_input'] = input_file
	options['mapping_file'] = mapping_file



if __name__ == "__main__":
	"""Run phylotyper functions

	Parses command-line arguments and calls appropriate
	functions

    """

   	subtype_config_file = 'builtin_subtypes.yaml'
   
	logging.basicConfig(level=logging.DEBUG)
	logger = logging.getLogger('phylotyper.main')

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

	# Parse .ini config file
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
		subtype_options['input'] = os.path.abspath(options.input)
		subtype_options['output_directory'] = os.path.abspath(options.output)
		subtype_options['fast'] = False

		if options.nt and (subtype_options['seq'] != 'nt'):
			msg = 'Sequence type of input does not match Phylotyper gene sequences for %s' % (scheme)
			raise Exception(msg)

		# Run pipeline
		subtype_pipeline(subtype_options, config)


	else:
		raise Exception("Unrecognized command")

  

