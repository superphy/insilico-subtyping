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
from collections import Counter

from config import PhylotyperOptions
from builtin_subtypes import SubtypeConfig
from phylotyper import Phylotyper
from tree.fasttree import FastTreeWrapper
from tree.seqaligner import SeqAligner
from tree.seq import SeqDict


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


logger = None

def align_all_sequences(input, output, summary, config):
    """Build MSA

    Args:
        input (str): Fasta file
        output (str): Output file for MSA
        summary (str): Trimming summary file
        config (obj): PhylotyperOptions object

    """

    logger.debug('Performing full alignment')

    aln = SeqAligner(config)
    aln.align(input, output)
    aln.trim(output, output, trimming_summary_file=summary)


def align_new_sequences(input, alignment, summary, output, config):
    """Add new sequences to existing MSA using
    profile alignment

    Args:
        input (str): Fasta file
        alignment (str): Aligned fasta file
        summary (str): Trimming summary file
        output (str): Output file for MSA
        config (obj): PhylotyperOptions object

    """

    logger.debug('Aligning genes to existing alignment')

    aln = SeqAligner(config)
    aln.add(input, alignment, output)
    aln.trim(output, output, trimming_summary_file=summary)


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
    

def predict_subtypes(options, config, assignments):
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

    pt = Phylotyper(config)

    assignment_dict = pt.subtype(treefile, subtfile, options['output_directory'], options['posterior_probability_plot'])

    # Load mapping
    rev_mapping = {}
    for x in csv.reader(open(options['mapping_file'],'r'), delimiter='\t'):
        rev_mapping[x[1]] = x[0]

    # Compute assignments
    cutoff = float(config.get('phylotyper','prediction_threshold'))
    logger.debug('Using posterior probability cutoff: %f' % (cutoff))

    for genome_id,subtup in assignment_dict.items():

        pp = subtup[1]
        subt = ','.join(subtup[0])

        pred = 'undetermined'
        print pp
        if pp > cutoff:
            pred = subt
        
        assignments.append([
            rev_mapping[genome_id],
            subt,
            str(pp),
            pred
        ])

    
    return(assignments)

def results_header():
    # Returns header string
    return '#genome\tsubtype\tprobability\tphylotyper_assignment\n'

def print_results(options, assignments):
    """Writes results to file.

    Args:
        options (dict): user defined settings from __main__
        assignments (list): list of lists with individual results

    """

    logger.debug('Writing results')

    assfile = options['result_file']
    mode = 'w'
    with open(assfile, mode) as outf:
        
        outf.write(results_header())
        for row in assignments:
            outf.write('\t'.join(row)+'\n')

    None


def subtype_pipeline(options, config, write_results=True):
    """Run phylotyper pipeline 

    Runs individual steps in phylotyper pipeline
    for internal subtype scheme
        
        Args:
            options (dict): user defined settings from __main__
            config (obj): PhylotyperConfig with .ini file settings
            write_results (bool): A switch so that in dispatch/multi mode, 
                individual results returned not printed to file.

    """

    # Define files
    alnfile = os.path.join(options['output_directory'], 'combined.aln')
    options['profile_alignment'] = alnfile
    oldalnfile = options['alignment']
    treefile = os.path.join(options['output_directory'], 'combined.tree')
    options['tree_file'] = treefile
    options['result_file'] = os.path.join(options['output_directory'], 'subtype_predictions.csv')
    if not options['noplots']:
        options['posterior_probability_plot'] = os.path.join(options['output_directory'], 'posterior_probability_tree.png')
    summary = os.path.join(options['output_directory'], 'alignment_trimming_summary.html')

    logger.info('Settings:\n%s' % (pprint.pformat(options)))
    logger.info('Config:\n%s' % (config.pformat()))

    assignments = [] # Phylotyper subtype assignments

    # Check if sequences match known subtyped sequences
    # Otherwise. rename sequences with unique ids, change input files
    num_remaining = prep_sequences(options, assignments)

    # Run phylotyper on remaining untyped input sequences
    if num_remaining > 0:    

        # Align
        align_new_sequences(options['input'], oldalnfile, summary, alnfile, config)

        # Compute tree
        nt = options['seq'] == 'nt'
        build_tree(alnfile, treefile, nt, options['fast'], config)

        # Predict subtypes & write to file
        predict_subtypes(options, config, assignments)

    if write_results:
        print_results(options, assignments)

    return assignments


def evaluate_subtypes(options, config):
    """Examine correlation of subtype in phylogenetic tree

    Wrapper around the Phylotyper.subtype method. Identifies
    predictions above significance cutoff. Writes results to 
    file.

    Args:
        options (dict): user defined settings from __main__
        config (obj): PhylotyperConfig with .ini file settings

    """

    logger.debug('Running phylotyper performance tests')

    # Define files
    subtfile = options['subtype']
    treefile = options['tree_file']

    pt = Phylotyper(config)

    pt.evaluate(treefile, subtfile, options['output_directory'])

    None


def build_pipeline(options, config):
    """Create and evaluate new reference alignment for subtyping

    User provides new reference set for subtyping. Build and
    refine alignment. Evaluate predictive ability of tree.

    Args:
        options (dict): user defined settings from __main__
        config (obj): PhylotyperConfig with .ini file settings

    """

    alnfile = options['alignment']
    tmpfile = os.path.join(options['output_directory'], 'tmp.fasta')
    treefile = options['tree_file'] = os.path.join(options['output_directory'], 'test.tree')
    summary = os.path.join(options['output_directory'], 'alignment_trimming_summary.html')

    logger.info('Settings:\n%s' % (pprint.pformat(options)))
    logger.info('Config:\n%s' % (config.pformat()))

    # Check sequence IDs
    check_gene_names(options)

    # Remove identical sequences
    logger.debug('Collapsing identical sequences')
    seqdict = SeqDict()
    seqdict.build(options['input'], options['subtype_orig'])
    # Output unique set
    seqdict.write(tmpfile, options['subtype'])
    # Save lookup object
    seqdict.store(options['lookup'])

    # Align
    align_all_sequences(tmpfile, alnfile, summary, config)

    # Compute tree
    nt = options['seq'] == 'nt'
    #build_tree(alnfile, treefile, nt, options['fast'], config)

    # Run evaluation
    #evaluate_subtypes(options, config)


def prep_sequences(options, identified):
    """Create temporary sequence names that are unique

    Autogenerate temporary name placeholders and rewrite
    inputs with these names.  Also flags sequences that are identical
    to subtyped sequences in the reference set.  These do not need
    to be run with phylotyper. Subtype assignment will be transfered
    from identical sequences.

    Sets filename key-values in options dict for new temp input file
    and mapping file.
        
    Args:
        options (dict): user defined settings from __main__
        identified (list): Inputs that match reference sequences will be appended to this list

    Returns:
        Integer indicated number of remaining unsubtyped sequences

    """

    logger.debug('Searching for identical sequences with known subtype')
    logger.debug('Renaming sequences')

    # Define files
    output_file = os.path.join(options['output_directory'], 'input.fasta')
    mapping_file = os.path.join(options['output_directory'], 'mapping.txt')
    input_file = options['input']

    # Load lookup object
    lookup = SeqDict()
    lookup.load(options['lookup'])

    i = 0
    pre = 'pt_'
    fasta_sequences = SeqIO.parse(open(input_file, 'r'),'fasta')
    with open(output_file, 'w') as outf, open(mapping_file, 'w') as mapf:
        for fasta in fasta_sequences:
            name, sequence = fasta.description, str(fasta.seq)
            found = lookup.find(sequence)

            if found:
                subt = found['subtype']
                hit = found['name']
                identified.append([
                    name,
                    subt,
                    'identical to {}'.format(hit),
                    subt
                    ])

            else:
                i += 1
                newname = '%s%i' % (pre, i)
                mapf.write('%s\t%s\n' % (name, newname))
                outf.write('>%s\n%s\n' % (newname, sequence))

    options['input'] = output_file
    options['user_input'] = input_file
    options['mapping_file'] = mapping_file

    return(i)


def check_gene_names(options):
    """Make sure gene names Fasttree and Mafft safe
        
        Args:
            options (dict): user defined settings from __main__

    """

    logger.debug('Checking gene sequence names')

  
    # Define files
    subtype_file = options['subtype_orig']
    input_file = options['input']

    # Check fasta file
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    msg1 = 'Invalid fasta file: '
    msg2 = 'Invalid subtype file: '

    uniq = Counter()
    reserved = set(':(), ') # Newick reserve characters


    for fasta in fasta_sequences:
        name = fasta.id
        desc = fasta.description[0:20]

        if len(fasta.description) > 300:
            raise Exception(msg1+"{} header is too long".format(desc))

        if uniq[name] > 0:
            raise Exception(msg1+"{} id in header is not unique".format(desc))
        uniq[name] += 1

        if any((c in reserved) for c in name):
            raise Exception(msg1+"invalid character in header {}".format(desc))

    fasta_sequences.close()

    # Check subtype file
    for row in csv.reader(open(subtype_file,'r'),delimiter='\t'):
        name = row[0]
        subt = row[1]

        if any((c in reserved) for c in subt):
            raise Exception(msg+"invalid character in subtype {}".format(subt))

        if len(subt) > 120:
            raise Exception(msg2+"{} subtype name is too long".format(subt))

        if uniq[name] != 1:
            raise Exception(msg2+"unknown gene {}".format(name))

        uniq[name] += 1


    for name in uniq:
        if uniq[name] != 2:
            raise Exception(msg2+"missing gene {}".format(name))


    return True


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

    # New subtype command
    new_parser = subparsers.add_parser('new', help='Build new subtype resource')
    new_parser.add_argument('config', action='store', help='Phylotyper config options file')
    new_parser.add_argument('input', action='store', help='Fasta input')
    new_parser.add_argument('subtypein', action='store', help='Input reference gene subtypes')
    new_parser.add_argument('subtypeout', action='store', help='Ouput reference gene subtypes')
    new_parser.add_argument('alignment', action='store', help='Reference gene alignment output')
    new_parser.add_argument('lookup', action='store', help='Reference gene sequence storage')
    new_parser.add_argument('output', action='store', help='Directory for evaluation result files')
    new_parser.add_argument('--aa', action='store_true', help='Amino acid sequences')
    new_parser.set_defaults(which='new')

    # # User-supplied subtype command
    # subtype_parser = subparsers.add_parser('custom', help='Predict subtype for user-defined subtype scheme')
    # subtype_parser.add_argument('config', action='store', help='Phylotyper config options file')
    # subtype_parser.add_argument('sequences', action='store', help='Fasta sequences of aligned reference genes for tree')
    # subtype_parser.add_argument('subtype', action='store', help='Reference gene subtypes')
    # subtype_parser.add_argument('input', action='store', help='Fasta input for unknowns')
    # subtype_parser.add_argument('output', action='store', help='Directory for Subtype predictions')
    # subtype_parser.add_argument('--nt', action='store_true', help='Nucleotide sequences')
    # subtype_parser.set_defaults(which='predict')

    # Builtin subtype command
    subtype_parser = subparsers.add_parser('subtype', help='Predict subtype for scheme provided in phylotyper')
    subtype_parser.add_argument('config', action='store', help='Phylotyper config options file')
    subtype_parser.add_argument('gene', action='store', help='Subtype gene name')
    subtype_parser.add_argument('input', action='store', help='Fasta input for unknowns')
    subtype_parser.add_argument('output', action='store', help='Directory for subtype predictions')
    subtype_parser.add_argument('--aa', action='store_true', help='Amino acid sequences')
    subtype_parser.add_argument('--noplots', action='store_true', help='Do not generate tree image file')
    subtype_parser.set_defaults(which='subtype')

    options = parser.parse_args()

    # Parse .ini config file
    config = PhylotyperOptions(options.config)

    if options.which == 'new':
        # Build & evaluate new subtype alignment

        # Check arguments

        # Fast mode of tree calculation
        fast = False

        # Check input file exists
        if not os.path.isfile(options.input):
            msg = 'Invalid/missing input file argument.'
            raise Exception(msg)

        # Check subtype file exists
        if not os.path.isfile(options.subtypein):
            msg = 'Invalid/missing subtype file argument.'
            raise Exception(msg)

        # Check output directory exists, if not create it if possible
        if not os.path.exists(options.output):
            os.makedirs(options.output)

        outdir = os.path.abspath(options.output)

        # Save options
        subtype_options = {}
        subtype_options['input'] = os.path.abspath(options.input)
        subtype_options['alignment'] = os.path.abspath(options.alignment)
        subtype_options['subtype_orig'] = os.path.abspath(options.subtypein)
        subtype_options['subtype'] = os.path.abspath(options.subtypeout)
        subtype_options['lookup'] = os.path.abspath(options.lookup)
        subtype_options['output_directory'] = outdir
        subtype_options['fast'] = False

        if options.aa:
            subtype_options['seq'] = 'aa'
        else:
            subtype_options['seq'] = 'nt'

        # Run pipeline
        build_pipeline(subtype_options, config)
    

    elif options.which == 'subtype':
        # Compute subtype for builtin scheme

        # Check arguments

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
        subtype_options['noplots'] = False

        if options.noplots:
            subtype_options['noplots'] = True

        if options.aa and (subtype_options['seq'] != 'aa'):
            msg = 'Sequence type of input does not match Phylotyper gene sequences for %s' % (scheme)
            raise Exception(msg)

        # Run pipeline
        subtype_pipeline(subtype_options, config)


    else:
        raise Exception("Unrecognized command")

  

