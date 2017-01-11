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
import re
from Bio import SeqIO
from collections import Counter, defaultdict

from config import PhylotyperOptions
from genome.loci import LociSearch
from subtypes_index import SubtypeConfig
from phylotyper import Phylotyper
from tree.fasttree import FastTreeWrapper
from tree.seqaligner import SeqAligner
from tree.seq import SeqDict, AlleleID, LociConcat


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


logger = logging.getLogger('phylotyper.main')

def align_all_sequences(inputs, outputs, superaln_output, summary, config):
    """Build MSA

    Args:
        inputs (list): List of Fasta files
        outputs (list): Output files for MSA
        superaln_output (str): Output file for concatenated alignment
        summary (str): Trimming summary file
        config (obj): PhylotyperOptions object

    """

    logger.debug('Performing full alignment')

    aln = SeqAligner(config)
    aln.malign(inputs, outputs, superaln_output)
    aln.trim(superaln_output, superaln_output, trimming_summary_file=summary)


def align_new_sequences(inputs, alignments, trim_summary, output, config):
    """Add new sequences to existing MSA using
    profile alignment

    Args:
        inputs (list): List of Fasta files
        alignment (str): List of aligned fasta files
        trim_summary (str): trimming summary file
        output (str): Output files for MSA
        config (obj): PhylotyperOptions object

    """

    logger.debug('Aligning genes to existing alignment')

    aln = SeqAligner(config)
    aln.madd(inputs, alignments, output)
    aln.trim(output, output, trimming_summary_file=trim_summary)


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
    

def predict_subtypes(treefile, subtypefile, plotfile, options, config):
    """Calls subtype method

    Wrapper around the Phylotyper.subtype method. Identifies
    predictions above significance cutoff. Writes results to 
    file.

    Args:
        treefile (str): Filepath to newick tree
        subtypefile (str): Filepath to tab-delim subtype assignments for leaves in tree
        plotfile (str|False): If provided, filepath to probability plot
        options (dict): user defined settings from __main__
        config (obj): PhylotyperConfig with .ini file settings

    Returns:
        dict with keys:
            subtype
            probability
            phylotyper_assignment

    """

    logger.debug('Running phylotyper')

    pt = Phylotyper(config)

    assignment_dict = pt.subtype(treefile, subtypefile, options['rate_matrix'], plotfile)

    # Compute assignments
    cutoff = float(config.get('phylotyper','prediction_threshold'))
    logger.debug('Using posterior probability cutoff: %f' % (cutoff))

    results = {}
    for genome_id,subtup in assignment_dict.items():

        pp = subtup[1]
        subt = ','.join(subtup[0])

        pred = 'non-significant/undetermined'
        if pp > cutoff:
            pred = subt

        results[genome_id] = {
            'subtype': subt,
            'probability': str(pp),
            'phylotyper_assignment': pred,
        }
        
    return(results)


def results_header():
    # Returns output fields
    return ['genome','tree_label','subtype','probability','phylotyper_assignment','loci']
    

def subtype_pipeline(options, config):
    """Run phylotyper pipeline 

    Runs individual steps in phylotyper pipeline
    for pre-built subtype scheme
        
        Args:
            options (dict): user defined settings from __main__
            config (obj): PhylotyperConfig with .ini file settings
            write_results (bool): A switch so that in dispatch/multi mode, 
                individual results returned not printed to file.

    """
    
    # Define files
    genome_search = False
    if 'genomes' in options:
        # Inputs are genome files
        # Need to define fasta files for each loci
        # found in the blast search
        locifiles = []
        for i in xrange(options['nloci']):
            file_name = 'search_results'
            locifiles.append(os.path.join(options['output_directory'], '{}.locus{}'.format(file_name, i)))
        options['input'] = locifiles
        genome_search = True

    refalnfiles = options['alignment']
    subtypefile = options['subtype']
    options['result'] = os.path.join(options['output_directory'], 'subtype_predictions.csv')

    logger.info('Settings:\n%s' % (pprint.pformat(options)))
    logger.info('Config:\n%s' % (config.pformat()))

    if genome_search:
        # Identify loci in input genomes
        seqtype = 'prot' if options['seq'] == 'aa' else 'nucl'
        detector = LociSearch(config, options['search_database'], None, seqtype, options['nloci'])

        for genome in options['genomes']:
            detector.search(genome, options['input'], append=True)

    # Predict subtypes
    with open(options['result'], 'w') as resfile:
        assignments = csv.DictWriter(resfile, fieldnames=results_header(), delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        assignments.writeheader()

        # Check if sequences match known subtyped sequences
        remaining = identical_sequences(options, assignments)

        # Run phylotyper on remaining untyped input sequences
        for genome in remaining:
            allele = 1

            for alleleset in remaining[genome]:

                # Setup input/output files
                infiles = []
                loci = 1

                # Only one gene for one genome submitted
                # don't need to split genes into individual files
                tree_label = "{}-allele{}".format(genome,allele)
                filename = tree_label.replace('|','_')
                
                # Only one gene, don't need to distinguish alleles
                if len(remaining[genome]) == 1:
                    tree_label = genome

                # Only one gene for one genome, don't need to split up input
                if options['ngenomes'] == 1 and len(remaining[genome]) == 1: 
                        infiles = options['input']
                else:
                    for s in alleleset.seqlist():
                        infile = os.path.join(options['output_directory'], "{}_loci{}_step2_alignment_input.fasta".format(filename, loci))
                        loci += 1
                        with open(infile, 'w') as outfh:
                            outfh.write('>{}\n{}\n'.format(tree_label, s))
                        infiles.append(infile)

                trimfile = os.path.join(options['output_directory'], "{}_step3_alignment_trimming_summary.html".format(filename))
                alnfile = os.path.join(options['output_directory'], "{}_step4_profile_alignment_output.fasta".format(filename))
                treefile = os.path.join(options['output_directory'], "{}_step5_subtype_tree.newick".format(filename))
                plotfile = os.path.join(options['output_directory'], "{}_step5_posterior_probability_tree.png".format(filename))

                # Run alignment on each locus
                align_new_sequences(infiles, refalnfiles, trimfile, alnfile, config)

                # Compute tree
                nt = options['seq'] == 'nt'
                build_tree(alnfile, treefile, nt, options['fast'], config)

                # Predict subtypes & write to file
                results = predict_subtypes(treefile, subtypefile, plotfile, options, config)
                print results
                print tree_label
                if not tree_label in results:
                    raise Exception("Phylotyper failed to complete")

                this_results = results[tree_label]
                this_results['genome'] = genome
                this_results['tree_label'] = tree_label
                this_results['loci'] = alleleset.iddump()
                assignments.writerow(this_results)

                allele += 1



def evaluate_subtypes(options, config, seqdict):
    """Examine correlation of subtype in phylogenetic tree

    Wrapper around the Phylotyper.subtype method. Identifies
    predictions above significance cutoff. Writes results to 
    file.

    Args:
        options (dict): user defined settings from __main__
        config (obj): PhylotyperConfig with .ini file settings
        seqdict (obj): SeqDict object

    """

    logger.debug('Running phylotyper new scheme evaluation')

    # Define files
    subtfile = options['subtype']
    treefile = options['tree']
    ratematfile = options['rate_matrix'] 

    pt = Phylotyper(config)

    pt.evaluate(treefile, subtfile, ratematfile, options['output_directory'], seqdict.accession_map())

    None


def build_pipeline(options, config):
    """Create and evaluate new reference alignment for subtyping

    User provides new reference set for subtyping. Build and
    refine alignment. Evaluate predictive ability of tree.

    Args:
        options (dict): user defined settings from __main__
        config (obj): PhylotyperConfig with .ini file settings

    """

    alnfiles = options['alignment']
    tmpfiles = []
    for i in xrange(options['nloci']):
        alnfile = os.path.basename(options['alignment'][i])
        tmpfiles.append(os.path.join(options['output_directory'], '{}.tmp{}'.format(alnfile, i)))

    tmpfile = os.path.join(options['output_directory'], 'tmpsuperaln.fasta')
    treefile = options['tree'] = os.path.join(options['output_directory'], 'test.tree')
    summary = os.path.join(options['output_directory'], 'alignment_trimming_summary.html')

    logger.info('Settings:\n%s' % (pprint.pformat(options)))
    logger.info('Config:\n%s' % (config.pformat()))

    # Check sequence IDs
    check_gene_names(options)

    # Create blast database for searching genomes
    seqtype = 'prot' if options['seq'] == 'aa' else 'nucl'
    LociSearch(config, options['search_database'], options['input'], seqtype, options['nloci'])

    # Remove identical sequences, 
    logger.debug('Collapsing identical sequences')
    seqdict = SeqDict(options['nloci'])
    seqdict.build(options['input'], options['subtype_orig'])
    # Output unique set
    seqdict.write(tmpfiles, options['subtype'])
    # Save lookup object
    seqdict.store(options['lookup'])

    # Align
    align_all_sequences(tmpfiles, alnfiles, tmpfile, summary, config)

    # Compute tree
    nt = options['seq'] == 'nt'
    build_tree(tmpfile, treefile, nt, options['fast'], config)

    # Run evaluation
    evaluate_subtypes(options, config, seqdict)


def identical_sequences(options, identified):
    """Looks for exact matches

    Sequences that are identical to subtyped sequences in 
    the reference set do not need to be run with phylotyper. 
    Subtype assignment will be transfered from identical sequences.

    Sets filename values in options dictionary for downstream 
    analyses.
        
    Args:
        options (dict): user defined settings from __main__
        identified (list): Inputs that match reference sequences will be appended to this list

    Returns:
        Integer indicated number of remaining unsubtyped sequences

    """

    logger.debug('Searching for identical sequences with known subtype')
   
    # Load lookup object
    lookup = SeqDict()
    lookup.load(options['lookup'])

    # Load all allele combinations for each loci
    # in all genomes
    remaining = defaultdict(list)
    concat = LociConcat()
    supersequences = concat.load(options['input'])

    for genome, typing_sequences in supersequences.iteritems():

        for alleleset in typing_sequences.iteralleles():

            sequence = alleleset.seqlist()
            found = lookup.find(sequence)

            if found:
                # Matches reference sequence

                subt = found['subtype']
                hit = found['name']
                identified.writerow({
                    'genome': genome,
                    'tree_label': 'not applicable',
                    'subtype': subt,
                    'probability': 'identical to {}'.format(hit),
                    'phylotyper_assignment': subt,
                    'loci': alleleset.iddump() # Use the power of json to ecode this fasta header string
                })

            else:
                # Need to run through phylotyper
                remaining[genome].append(alleleset)
            
    return(remaining)


def check_gene_names(options):
    """Make sure gene names Fasttree and Mafft safe
        
        Args:
            options (dict): user defined settings from __main__

    """

    logger.debug('Checking gene sequence names')

    # Define files
    subtype_file = options['subtype_orig']
    input_files = iter(options['input'])

    # Check first fasta file
    i = 1
    input_file = next(input_files)
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    
    uniq = Counter()
    genomes = Counter()
    
    for fasta in fasta_sequences:
        name = fasta.id
        desc = fasta.description[0:20]

        allele = validate_fasta_header(name)

        genomes[allele.genome] = i
        if uniq[str(allele)] > 0:
            raise Exception("allele id in {} in file {} is not unique (allele: {})".format(desc, input_file, str(allele)))
        uniq[str(allele)] += 1

    fasta_sequences.close()
    
    # Check subtype file
    i += 1
    reserved = set(':(), ') # Newick reserved characters
    for row in csv.reader(open(subtype_file,'r'),delimiter='\t'):
        name = row[0]
        subt = row[1]

        if any((c in reserved) for c in subt):
            raise Exception("invalid character in subtype {}".format(subt))

        if len(subt) > 20:
            raise Exception("{} subtype name is too long".format(subt))

        if genomes[name] != 1:
            raise Exception("unknown genome {} in subtype file".format(name))
        genomes[name] = i

    for name in genomes:
        if genomes[name] != 2:
            raise Exception("missing genome {} in subtype file".format(name))

    # Check remaining input files
    for input_file in input_files:
        fasta_sequences = SeqIO.parse(open(input_file),'fasta')
        uniq = Counter()
        i += i
        for fasta in fasta_sequences:
            allele = validate_fasta_header(fasta.id)

            if uniq[str(allele)] > 0:
                raise Exception("allele id in {} in file {} is not unique (allele: {})".format(desc, input_file, str(allele)))
            uniq[str(allele)] += 1

            if not allele.genome in genomes:
                raise Exception("unknown genome {} in file {}".format(allele.genome, input_file))
            genomes[allele.genome] = i

        for genome in genomes:
            if genomes[genome] != i:
                raise Exception("missing genome entry {} in file {}".format(genome, input_file))


    return True


def validate_fasta_header(name):
    """Check that fasta header is suitable for phylotyper downstream
    applications. 
        
        Args:
            id (str): Fasta ID check

        Returns:
            AlleleID namedtuple

        Raises:
            Exception if format unsuitable

    """


    reserved = set(':(), ') # Newick reserved characters

    if len(name) > 40:
        raise Exception("{} id in fasta header is too long".format(name))

    if any((c in reserved) for c in name):
        raise Exception("invalid character in fasta header id {}".format(name))

    return LociConcat.parse_fasta_id(name)


# Modify filename, not extension
def modify_filename(filename, uid, altext=None):
    name, ext = os.path.splitext(filename)
    if not altext is None:
        ext = altext
    return "{name}_{uid}{ext}".format(name=name, uid=uid, ext=ext)


if __name__ == "__main__":
    """Run phylotyper functions

    Parses command-line arguments and calls appropriate
    functions

    """

    logging.basicConfig(level=logging.INFO)
   
    # Parse command-line args
    # Phylotyper functions are broken up into commands
    # Each command has its own options and subparser
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='commands')

    # New subtype command
    new_parser = subparsers.add_parser('new', help='Build new subtype scheme')
    new_parser.add_argument('gene', action='store', help='Subtype name')
    new_parser.add_argument('subtype', action='store', help='Subtypes for reference sequences')
    new_parser.add_argument('results', action='store', help='Directory for evaluation result files')
    new_parser.add_argument('ref', nargs='+', metavar='reference', help='Fasta input(s) for reference sequences')
    new_parser.add_argument('--aa', action='store_true', help='Amino acid')
    new_parser.add_argument('--index', help='Specify non-default location of YAML-formatted file index for pre-built subtype schemes')
    new_parser.add_argument('--config', action='store', help='Phylotyper config options file')
    new_parser.set_defaults(which='new')

    # Builtin subtype command with genome as input
    genome_parser = subparsers.add_parser('genome', help='Predict subtype for scheme provided in phylotyper for genome input')
    genome_parser.add_argument('gene', action='store', help='Subtype gene name')
    genome_parser.add_argument('output', action='store', help='Directory for subtype predictions')
    genome_parser.add_argument('inputs', nargs='+', help='Fasta input for genomes')
    genome_parser.add_argument('--index', help='Specify non-default location of YAML-formatted file index for pre-built subtype schemes')
    genome_parser.add_argument('--aa', action='store_true', help='Amino acid sequences')
    genome_parser.add_argument('--noplots', action='store_true', help='Do not generate tree image file')
    genome_parser.add_argument('--config', action='store', help='Phylotyper config options file')
    genome_parser.set_defaults(which='genome')

    # Builtin subtype command with gene as input
    subtype_parser = subparsers.add_parser('subtype', help='Predict subtype for scheme provided in phylotyper')
    subtype_parser.add_argument('gene', action='store', help='Subtype gene name')
    subtype_parser.add_argument('input', action='store', help='Fasta input for unknowns')
    subtype_parser.add_argument('output', action='store', help='Directory for subtype predictions')
    subtype_parser.add_argument('--index', help='Specify non-default location of YAML-formatted file index for pre-built subtype schemes')
    subtype_parser.add_argument('--aa', action='store_true', help='Amino acid sequences')
    subtype_parser.add_argument('--noplots', action='store_true', help='Do not generate tree image file')
    subtype_parser.add_argument('--config', action='store', help='Phylotyper config options file')
    subtype_parser.set_defaults(which='subtype')

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

    # Default index location
    subtype_config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'subtypes_index.yaml')

    # Non-default location
    if options.index:
        if not os.path.isfile(options.index):
            msg = 'Invalid/missing index file option.'
            raise Exception(msg)
        subtype_config_file = options.index
    stConfig = SubtypeConfig(subtype_config_file)
    
    if options.which == 'new':
        # Build & evaluate new subtype alignment

        # Check arguments

        # Fast mode of tree calculation
        fast = False

        # Check input files exists
        n_loci = 0
        for f in options.ref:
            if not os.path.isfile(f):
                msg = 'Invalid/missing input file argument.'
                raise Exception(msg)
            n_loci += 1

        # Check subtype file exists
        if not os.path.isfile(options.subtype):
            msg = 'Invalid/missing subtype file argument.'
            raise Exception(msg)

        # Check output directory exists, if not create it if possible
        if not os.path.exists(options.results):
            os.makedirs(options.results)
        outdir = os.path.abspath(options.results)

        # Create subtype directory & file names
        subtype_options = stConfig.create_subtype(options.gene, n_loci, options.aa)

        # Save additional build options inputted by user
        subtype_options['input'] = [os.path.abspath(f) for f in options.ref]
        subtype_options['subtype_orig'] = os.path.abspath(options.subtype)
        subtype_options['output_directory'] = outdir
        subtype_options['nloci'] = n_loci

        subtype_options['fast'] = False

        # Run pipeline
        build_pipeline(subtype_options, config)

        # Update subtype YAML file
        stConfig.save()

    elif options.which == 'genome':
        # Compute subtype for builtin scheme
        # Genome input

        # Check arguments

        # Fast mode of tree calculation
        fast = False

        # Check input file exists
        n_genomes = 0
        for f in options.inputs:
            if not os.path.isfile(f):
                msg = 'Invalid/missing input file argument.'
                raise Exception(msg)
            n_genomes += 1

        # Check output directory exists, if not create it if possible
        if not os.path.exists(options.output):
            os.makedirs(options.output)

        # Load requested subtype data files
        scheme = options.gene
        subtype_options = stConfig.get_subtype_config(scheme)

        # Add pipeline options
        subtype_options['genomes'] = [os.path.abspath(f) for f in options.inputs]
        subtype_options['output_directory'] = os.path.abspath(options.output)
        subtype_options['ngenomes'] = n_genomes
        subtype_options['fast'] = False
        subtype_options['noplots'] = False

        if options.noplots:
            subtype_options['noplots'] = True

        if options.aa and (subtype_options['seq'] != 'aa'):
            msg = 'Sequence type of input does not match Phylotyper gene sequences for %s' % (scheme)
            raise Exception(msg)

        # Run pipeline
        subtype_pipeline(subtype_options, config)

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
        scheme = options.gene
        subtype_options = stConfig.get_subtype_config(scheme)

        # Add pipeline options
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

  

