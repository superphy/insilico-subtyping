#!/usr/bin/env python

"""Dispatcher for running phylotyper on multiple inputs

Script for iterating over multiple loci in a fasta file and running
phylotyper

Examples:
    To run subtyping routine on Stx1 genes:

        $ python dispatc.py gene 

"""


import argparse
import logging
import os
import re
from Bio import SeqIO

import main
from config import PhylotyperOptions
from builtin_subtypes import SubtypeConfig

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


logger = None

if __name__ == "__main__":
    """Run phylotyper on each input

    Parses command-line arguments. Phylotyper is run on each
    sequence in an input fasta file. Input files and output
    file names are based on the fasta header


    """

    # Parse command-line args
    parser = argparse.ArgumentParser()
    
    # Builtin subtype command
    parser.add_argument('config', action='store', help='Phylotyper config options file')
    parser.add_argument('gene', action='store', help='Subtype gene name')
    parser.add_argument('input', action='store', help='Fasta input for unknowns')
    parser.add_argument('output', action='store', help='Directory for subtype predictions')
    parser.add_argument('--aa', action='store_true', help='Amino acid sequences')
    parser.add_argument('--noplots', action='store_true', help='Do not generate tree image file')

    options = parser.parse_args()

    # Check input file exists
    if not os.path.isfile(options.input):
        msg = 'Invalid/missing input file argument.'
        raise Exception(msg)

    # Check output directory exists, if not create it if possible
    if not os.path.exists(options.output):
        os.makedirs(options.output)

   
    # Initialize common options
    subtype_config_file = 'builtin_subtypes.yaml'
   
    # logging
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('phylotyper.dispatch')
    main.logger = logger

    # Config settings
    config = PhylotyperOptions(options.config)

    # Fast mode of tree calculation
    fast = False

    # Load requested subtype data files
    stConfig = SubtypeConfig(subtype_config_file)
    scheme = options.gene

    subtype_options = stConfig.get_subtype_config(scheme)
    subtype_options['fast'] = False
    subtype_options['noplots'] = False

    if options.noplots:
        subtype_options['noplots'] = True

    if options.aa and (subtype_options['seq'] != 'aa'):
        msg = 'Sequence type of input does not match Phylotyper gene sequences for %s' % (scheme)
        raise Exception(msg)

    # Root directory for outputs
    root_dir = os.path.abspath(options.output)

    # Iterate through fasta inputs, run phylotyper on each
    is_unique = set()
    fasta_sequences = SeqIO.parse(open(options.input, 'r'),'fasta')
    result_file = os.path.join(root_dir, 'merged_subtype_predictions.csv')

    with open(result_file, 'w') as resfh:
        
        # Print header to result file
        resfh.write(main.results_header())
       
        for record in fasta_sequences:

            # Get unique descriptor for sequence

            # Try processing standard NCBI format
            result = re.search(r'^\w{3}\|(\S+)', record.id)
            uid = result.group(1)

            if not uid:
                # Format not recognized, use entire fasta id
                uid = record.id
                
            uid = uid.replace('|', '-')
            if uid in is_unique:
                raise Exception('ID {} in fasta header {} is not unique.'.format(uid, record.id))

            is_unique.add(uid)

            input_file = os.path.join(root_dir, 'single_record.fasta')
            with open(input_file, 'w') as outfh:
                fasta = '\n>{}\n{}'.format(uid, record.seq)
                outfh.write(fasta)
                outfh.close()

            subtype_options['input'] = input_file
            subtype_options['output_directory'] = root_dir

            # Run pipeline
            assignment = main.subtype_pipeline(subtype_options, config, False)

            # Rename output files
            rename = []
            if not subtype_options['noplots']:
                rename.append('posterior_probability_plot')

            for f in rename:
                if os.path.exists(subtype_options[f]):
                    filename = os.path.basename(subtype_options[f])
                    newfilename = os.path.join(root_dir, '{}_{}'.format(uid, filename))
                    logger.debug('Moving file {} to {}'.format(subtype_options[f], newfilename))
                    os.rename(subtype_options[f], newfilename)
                else:
                    logger.debug('File {} not found... skipping'.format(subtype_options[f]))

            # Save results to file
            for row in assignment:
                resfh.write('\t'.join(row)+'\n')

       


  

