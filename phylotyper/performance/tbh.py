#!/usr/bin/env python

"""Top-BLAST hit detection

Compare tree-based method to straight-forward sequence-based approach

"""

from __future__ import division

import csv
import logging
import re
import tempfile
from Bio import SeqIO
from collections import namedtuple
from subprocess import check_output, CalledProcessError, STDOUT

def parse_fasta_id(id):
        """Extract Genome and gene identifier fasta header

        From the ID component of a fasta header, determine
        for a particular genome copy, whether there is multiple 
        alleles of a gene in the genome.  Format is:

            >[lcl|]genome|allele

            Args:
                id (str): Fasta ID

            Returns:
                AlleleID namedtuple

        """

        # Remove database id
        id = re.sub(r'^(?:lcl\|)|(?:gi\|)', '', id)

        # Extract genome and gene Ids if they exist
        parts = re.split(r'\|', id, maxsplit=1)
        
        return AlleleID(*parts)


# Factory for AlleleID class
# Stores genome and allele ID components
class AlleleID(namedtuple('AlleleID', ['genome', 'allele'])):
    __slots__ = ()
    def __new__(cls, genome, allele=None):
        # add default values
        return super(AlleleID, cls).__new__(cls, genome, allele)
    def __str__(self):
        # String
        allele = '|'+str(self.allele) if self.allele else ''
        return '{}{}'.format(self.genome, allele)

# Factory for BlastRow class
# Stores outfmt 6 blast results
class BlastRow(namedtuple('BlastRow', ['qseqid', 'sseqid', 'pident',
    'qlen', 'length', 'mismatch', 'gapopen', 'evalue', 'bitscore'])):
    __slots__ = ()
    def __new__(cls, qseqid='Working', sseqid=None, pident=0, qlen=0, length=0, 
        mismatch=0, gapopen=0, evalue=1, bitscore=0):
        # add default values
        return super(BlastRow, cls).__new__(cls, qseqid, sseqid, float(pident), int(qlen), int(length), 
            int(mismatch), int(gapopen), float(evalue), float(bitscore))


class Blaster(object):
    """Blast-based detection of alleles

    """
    
    def __init__(self, config, database, sequence_type):
        """Constructor

        Either a pre-existing blast database or one or more loci fasta files are required as arguments.
        Each input fasta file must contain one or more examples of a single loci. Multiple loci are input
        as multiple fasta files.

        Args:
            config (PhylotyperConfig): Instance of PhylotyperConfig object
            database (str): Blast database location
            sequence_type (str): nucl or prot

        Returns:
            None

        """

        self.logger = logging.getLogger('phylotyper.performance.tbh.Blaster')

        if sequence_type == 'nucl' or sequence_type == 'prot':
            self.seqtype = sequence_type
        else:
            raise Exception('Invalid sequence_type argument {}'.format(sequence_type))

        self._alignment_coverage = 0.95
        self._percent_identity = 0.90
        self._evalue = 1e-5

        if self.seqtype == 'prot':
            self._percent_identity = 0.80

        # Default options
        self._blast_options = {
            'evalue': self._evalue,
            #'evalue': 1,
            'db': database,
            'outfmt': "\"6 qseqid sseqid pident qlen length mismatch gapopen evalue bitscore\""
        }

        filt = 'dust' if self.seqtype == 'nucl' else 'seg'
        self._blast_specific_options = {
            filt: 'no'
        }

        self._makeblastdb_options = {
            'input_type': 'fasta',
            'dbtype': self.seqtype,
            'out': database,
        }

        self._makeblastdb_exe = config.get('external', 'makeblastdb')
        blastexe = 'blastn' if self.seqtype == 'nucl' else 'blastx'
        self._blast_exe = config.get('external', blastexe)
        self._blastdbcmd_exe = config.get('external', 'blastdbcmd')

      
    # Prepare
    def make_db(self, fastafile, subtypefile):
        """Make blast database with associated properties

        Args:
            fastafile (str): Filepath to fasta sequences.
            subtypefile (str)

        Returns:
            None

        """

        # Load subtypes into memory
        subtypes = {}
        reserved = set(':#') # Newick reserved characters
        for row in csv.reader(open(subtypefile,'r'), delimiter='\t'):
            name = row[0]
            subt = row[1]

            if any((c in reserved) for c in subt):
                raise Exception("invalid character in subtype {}".format(subt))

            if any((c in reserved) for c in name):
                raise Exception("invalid character in genome name {}".format(name))

            if len(subt) > 20:
                raise Exception("{} subtype name is too long".format(subt))

            subtypes[name] = subt

        # Build properties into blast name
        with tempfile.NamedTemporaryFile() as tmpfh:

            fasta = SeqIO.parse(fastafile, 'fasta')
            for record in fasta:
                allele = parse_fasta_id(record.id)
                subt = subtypes[allele.genome]
                if not subt:
                    raise Exception("No subtype for genome {}".format(allele.genome))
                newheader = self.encode(record.id, subt)
                tmpfh.write(">{}\n{}\n".format(newheader, record.seq))

            tmpfh.flush()

            opts = self._makeblastdb_options
            opts['in'] = tmpfh.name

            cmd = "{blastcmd} -in {in} -out {out} -dbtype {dbtype} -input_type {input_type} ".format(
                blastcmd=self._makeblastdb_exe, **opts)

            self.logger.debug('Creating new blast database with input file: {}'.format(','.join(fasta)))
            self.logger.debug('Running Blast command-line command:\n{}'.format(cmd))

            try:
                check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
            except CalledProcessError as e:
                msg = "Blaster failed: {} (return code: {}).".format(e.output, e.returncode)                                                                                                   
                raise Exception(msg)


    def encode(self, string, prop):
        """Encode property in blast db name

        """
        return "{name}::{prop}".format(name=string, prop=prop)


    def decode(self, string):
        """Extract property from blast db name

        """
        parts = string.split('::',1)
        return parts[1]


    def search(self, genome):
        """Run Blast with genome as input

        Args:
            genome (str): Filepath to input fasta sequence

        Raises Exception if blast search fails

        """

        with tempfile.NamedTemporaryFile() as tmpfh:
            opts = self._blast_options
            opts['query'] = genome
            opts['out'] = tmpfh.name

            cmd = "{blastcmd} -evalue {evalue} -outfmt {outfmt} -db {db} -query {query} -out {out}".format(
                blastcmd=self._blast_exe, **opts)
            if self._blast_specific_options:
                cmd += ' '+' '.join(['-'+k+' '+v for k,v in self._blast_specific_options.items()])

            self.logger.debug('Running Blast')
            self.logger.debug('Running Blast command-line command:\n{}'.format(cmd))
            try:
                check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
            except CalledProcessError as e:
                msg = "Blaster failed: {} (return code: {}).".format(e.output, e.returncode)                                                                                                   
                raise Exception(msg)

            self.logger.debug('Parsing Blast output')
            
            tophit = BlastRow()
            for row in csv.reader(open(tmpfh.name,'r'), delimiter='\t'):

                blastresult = BlastRow(*row)

                if blastresult.bitscore > tophit.bitscore:
                    tophit = blastresult

            self.logger.debug('Best match: {}'.format(tophit))

            # Does tophit pass the criteria
            if ( tophit.evalue < self._evalue and
                tophit.pident > self._percent_identity and
                ( float(tophit.length) / float(tophit.qlen) ) > self._alignment_coverage):
                self.logger.debug('Match passes criteria')
                return self.decode(tophit.sseqid), tophit

            else:
                return None, None
               
    