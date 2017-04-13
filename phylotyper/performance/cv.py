#!/usr/bin/env python

"""Cross-validation test main

Script for running cross-valiation performance test

Examples:
    To run subtyping routine on Stx1 genes:

        $ python main.py subtype ../phylotyper_example.ini ecoli_stx1 test/ecoli_stx1.ffn output/test/

"""


import argparse
import csv
import inspect
import logging
import os
import random
import shutil
import sys
import tempfile
from Bio import SeqIO

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

from config import PhylotyperOptions
from tbh import parse_fasta_id, Blaster


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2017, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


logger = logging.getLogger('phylotyper.performance.cv')

def write(fastafile, queryfile, subjectfile, query):

     with open(queryfile, 'w') as qfh, open(subjectfile, 'w') as sfh:
        fasta = SeqIO.parse(fastafile, 'fasta')

        for f in fasta:
            if not f.id in query:
                sfh.write(">{id}\n{seq}\n".format(id=f.description, seq=str(f.seq).replace('-','')))
            else:
                qfh.write(">{id}\n{seq}\n".format(id=f.description, seq=str(f.seq).replace('-','')))


def loo(config, options):
    """Leave-One-Out Cross Validation

    Args:
        config (PhylotyperConfig): Instance of PhylotyperConfig object
        options (namedtuple): Command-line arguments

    """

    # Load subtypes
    subtypes = load_subtypes(options)[0]

    # Load genomes
    genomes = set()
    fasta = SeqIO.parse(options.genes, 'fasta')
    for record in fasta:
        genomes.add(record.id)
        allele = parse_fasta_id(record.id)
        subt = subtypes[allele.genome]
    
        if not subt:
            raise Exception("No subtype for genome {}".format(allele.genome))

    # Leave each genome out
    results = []
    n = 0
    failed_blasts = []
    for g in genomes:
        allele = parse_fasta_id(g)
        truevalue = subtypes[allele.genome]

        dirname = tempfile.mkdtemp(prefix='ptcv')
        logger.debug('Temp directory: {}'.format(dirname))

        dbfile = os.path.join(dirname, 'blastdb')
        subjectfile = os.path.join(dirname, 'subject.fasta')
        queryfile = os.path.join(dirname, 'query.fasta')

        # Test input
        write(options.genes, queryfile, subjectfile, [g])

        # Make blast db
        blastyper = Blaster(config, dbfile, 'nucl')
        blastyper.make_db(subjectfile, options.subtype)

        # Run test
        result, tophit = blastyper.search(queryfile)
        logger.info('Subtype prediction for {} is {}'.format(g, result))

        results.append((truevalue, result))
        
        if result:
            if not result == truevalue:
                failed_blasts.append(tophit)

        shutil.rmtree(dirname)
        n += 1

    stats = compute_stats(results)

    outputfile = os.path.join(options.results, 'performance_results.csv')
    with open(outputfile, 'w') as csvfile:
        csvwriter =  csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerow(['class', 'tp', 'fp', 'fn', 'precision', 'recall', 'fscore', 'support'])

        for s in stats:
            csvwriter.writerow(s)

    logger.info('FAILED RESULTS:')
    for f in failed_blasts:
        logger.info(f)


def kfcv(config, options, k=10):
    """K-folds Cross Validation

    Args:
        config (PhylotyperConfig): Instance of PhylotyperConfig object
        options (namedtuple): Command-line arguments
        k (int): Number of partitions for training set

    """

    # Load subtypes
    subtypes, subtype_classes = load_subtypes(options)

    # Load genomes
    genomes = set()
    fasta = SeqIO.parse(options.genes, 'fasta')
    for record in fasta:
        genomes.add(record.id)
        allele = parse_fasta_id(record.id)
        subt = subtypes[allele.genome]
    
        if not subt:
            raise Exception("No subtype for genome {}".format(allele.genome))

    # Result objects
    results = { k: [] for k in subtype_classes }
    maxresults = []
    ng = len(genomes)
    niter = 100
    
    kfold = ng // k

    if kfold < 5:
        raise Exception("Training set too small")
    
    for i in xrange(niter):

        testg = random.sample(genomes, kfold)

        dirname = tempfile.mkdtemp(prefix='ptcv')
        logger.debug('Temp directory: {}'.format(dirname))

        dbfile = os.path.join(dirname, 'blastdb')
        subjectfile = os.path.join(dirname, 'subject.fasta')
        queryfile = os.path.join(dirname, 'query.fasta')

        # Test input
        write(options.genes, queryfile, subjectfile, testg)

        # Make blast db
        blastyper = Blaster(config, dbfile, 'nucl')
        blastyper.make_db(subjectfile, options.subtype)

        # Run test
        hits = blastyper.multisearch(queryfile)
        
        # Append results
        for g in testg:
            allele = parse_fasta_id(g)
            truevalue = subtypes[allele.genome]

            if g in hits:
                
                besthit = None
                for h in hits[g]:

                    if besthit:
                        if h.bitscore > besthit.bitscore:
                            besthit = h
                    else:
                        besthit = h

                    label = 0
                    if truevalue == h.sprop:
                        label = 1

                    results[truevalue].append([h.pident, label])

                label = 0
                if truevalue == besthit.sprop:
                    label = 1
                maxresults.append([besthit.pident, label])

            else:
                # Blast returned no results,
                # so get a FN right out of the gates
                results[truevalue].append([0, 1])
                maxresults.append([0, 1])
                logger.debug('No BLAST result for {}'.format(g))

        i+=1
        logger.debug('Iteration {}'.format(i))


    outputfile = os.path.join(options.results, 'performance_result_tbh_kfcv.csv')
    with open(outputfile, 'w') as csvfile:
        csvwriter =  csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)

        scores = [ r[0] for r in maxresults ]
        labels = [ r[1] for r in maxresults ]
        csvwriter.writerow(scores)
        csvwriter.writerow(labels)
        
        for st in results:
            scores = [ r[0] for r in results[st] ]
            labels = [ r[1] for r in results[st] ]
            csvwriter.writerow(st)
            csvwriter.writerow(scores)
            csvwriter.writerow(labels)

    

def lco(config, options, pident_cutoff=90.0):
    """Leave-Class-Out Cross Validation

    Args:
        config (PhylotyperConfig): Instance of PhylotyperConfig object
        options (namedtuple): Command-line arguments

    """

    # Load subtypes
    subtypes, subtype_classes = load_subtypes(options)

    # Load genomes
    genomes = set()
    fasta = SeqIO.parse(options.genes, 'fasta')
    for record in fasta:
        genomes.add(record.id)
        allele = parse_fasta_id(record.id)
        subt = subtypes[allele.genome]
    
        if not subt:
            raise Exception("No subtype for genome {}".format(allele.genome))

    # Leave each subtype class out
    results = []

    print subtype_classes

    for c in subtype_classes:
        logger.debug('Running {}'.format(c))

        testg = []
        for g in genomes:
            allele = parse_fasta_id(g)
            truevalue = subtypes[allele.genome]

            if truevalue == c:
                testg.append(g)
        logger.debug('Test set: {}'.format(testg))

        dirname = tempfile.mkdtemp(prefix='ptcv')
        logger.debug('Temp directory: {}'.format(dirname))

        dbfile = os.path.join(dirname, 'blastdb')
        subjectfile = os.path.join(dirname, 'subject.fasta')
        queryfile = os.path.join(dirname, 'query.fasta')

        # Test input
        write(options.genes, queryfile, subjectfile, testg)

        # Make blast db
        blastyper = Blaster(config, dbfile, 'nucl')
        blastyper.make_db(subjectfile, options.subtype)

        # Run test
        hits = blastyper.multisearch(queryfile)

        fp=0
        for g in testg:
            allele = parse_fasta_id(g)
            truevalue = subtypes[allele.genome]

            if g in hits:
            
                besthit = None
                for h in hits[g]:
                    if besthit:
                        if h.bitscore > besthit.bitscore:
                            besthit = h
                    else:
                        besthit = h

                print besthit
                if besthit.pident > pident_cutoff:
                    fp +=1

        results.append([c, fp, len(testg)])

    logger.debug('Iteration {}'.format(c))

    outputfile = os.path.join(options.results, 'performance_result_tbh_lco.csv')
    with open(outputfile, 'w') as csvfile:
        csvwriter =  csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        
        for r in results:
            csvwriter.writerow(r)
        


def compute_stats(results):
    # Compute precision etc

    stats = []
    beta = float(1)

    classes = set()
    n = len(results)
    for r in results:
        if r[0]:
            classes.add(r[0])
        if r[1]:
            classes.add(r[1])

    for c in classes:

        print c

        tp, fp, fn = 0, 0, 0
        for r in results:
            if r[0] == c:
                # Positive

                if r[1] == c:
                    # TP
                    tp += 1
                else:
                    # FN
                    fn += 1

            else:
                # Negative

                if r[1] == c:
                    # FP
                    fp += 1

        tp = float(tp)
        fp = float(fp)
        fn = float(fn)
        n = tp+fn
        p = tp+fp

        if n == 1:
            fp == 0
            recall = 0
        else:
            recall = tp/n

        if p == 0:
            precision = 0
            fscore = 0

        else:
            precision = tp/p
            
            denom = (beta**2*precision)+recall
            fscore = (1+beta**2)*precision*recall/denom
    
        stats.append([c, tp, fp, fn, precision, recall, fscore, n])

    totals = [0, 0, 0, 0, 0, 0, 0]
    for r in stats:
        totals = [ sum(x) for x in zip(r[1:], totals) ]
    
    totals.insert(0, 'total')
    totals[4] = float(totals[1]) / (float(totals[1]) + float(totals[2]))
    totals[5] = float(totals[1]) / (float(totals[1]) + float(totals[3]))
    denom = (beta**2*totals[4])+totals[5]
    totals[6] = (1+beta**2)*totals[4]*totals[5]/ denom

    stats.append(totals)

    return stats


def load_subtypes(options):
    # Load subtypes

    subtypes = {}
    classes = set()

    for row in csv.reader(open(options.subtype,'r'), delimiter='\t'):
        name = row[0]
        subt = row[1]

        subtypes[name] = subt
        classes.add(subt)

    return (subtypes, classes)

        

if __name__ == "__main__":
    """Run test

    Parses command-line arguments and calls appropriate
    functions

    """

    logging.basicConfig(level=logging.DEBUG)
   
    # Parse command-line args
    # Phylotyper functions are broken up into commands
    # Each command has its own options and subparser
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='commands')

    # Leave one out command
    loo_parser = subparsers.add_parser('loo', help='Leave-one-out cross-validation')
    loo_parser.add_argument('genes', action='store', help='Fasta input')
    loo_parser.add_argument('subtype', action='store', help='Subtypes for fasta input')
    loo_parser.add_argument('results', action='store', help='Directory for evaluation result files')
    loo_parser.add_argument('--config', action='store', help='Phylotyper config options file')
    loo_parser.set_defaults(which='loo')

    # Leave subtype out command
    lco_parser = subparsers.add_parser('lco', help='Leave-class-out cross-validation')
    lco_parser.add_argument('genes', action='store', help='Fasta input')
    lco_parser.add_argument('subtype', action='store', help='Subtypes for fasta input')
    lco_parser.add_argument('results', action='store', help='Directory for evaluation result files')
    lco_parser.add_argument('--config', action='store', help='Phylotyper config options file')
    lco_parser.set_defaults(which='lco')

    # KFCV command
    kfcv_parser = subparsers.add_parser('kfcv', help='Leave-class-out cross-validation')
    kfcv_parser.add_argument('genes', action='store', help='Fasta input')
    kfcv_parser.add_argument('subtype', action='store', help='Subtypes for fasta input')
    kfcv_parser.add_argument('results', action='store', help='Directory for evaluation result files')
    kfcv_parser.add_argument('--config', action='store', help='Phylotyper config options file')
    kfcv_parser.set_defaults(which='kfcv')

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

    if options.which == 'loo':
        # Run Leave-One-Out CV

        # Check arguments

        # Check input files exists
        if not os.path.isfile(options.genes):
            msg = 'Invalid/missing genes file argument.'
            raise Exception(msg)
           
        # Check subtype file exists
        if not os.path.isfile(options.subtype):
            msg = 'Invalid/missing subtype file argument.'
            raise Exception(msg)

        # Check output directory exists, if not create it if possible
        if not os.path.exists(options.results):
            os.makedirs(options.results)
        outdir = os.path.abspath(options.results)
        options.results = outdir
        
        loo(config, options)

    elif options.which == 'kfcv':
        # Run K-folds CV

        # Check arguments

        # Check input files exists
        if not os.path.isfile(options.genes):
            msg = 'Invalid/missing genes file argument.'
            raise Exception(msg)
           
        # Check subtype file exists
        if not os.path.isfile(options.subtype):
            msg = 'Invalid/missing subtype file argument.'
            raise Exception(msg)

        # Check output directory exists, if not create it if possible
        if not os.path.exists(options.results):
            os.makedirs(options.results)
        outdir = os.path.abspath(options.results)
        options.results = outdir
        
        kfcv(config, options)

    elif options.which == 'lco':
        # Run Leave-Class-Out CV

        # Check arguments

        # Check input files exists
        if not os.path.isfile(options.genes):
            msg = 'Invalid/missing genes file argument.'
            raise Exception(msg)
           
        # Check subtype file exists
        if not os.path.isfile(options.subtype):
            msg = 'Invalid/missing subtype file argument.'
            raise Exception(msg)

        # Check output directory exists, if not create it if possible
        if not os.path.exists(options.results):
            os.makedirs(options.results)
        outdir = os.path.abspath(options.results)
        options.results = outdir
        
        lco(config, options)

    else:
        raise Exception("Unrecognized command")

  

