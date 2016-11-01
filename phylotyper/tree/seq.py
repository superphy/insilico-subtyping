#!/usr/bin/env python

"""Bio Sequences Utilities

Tools for manipulating / validating input DNA and protein sequences

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

import csv
import hashlib
import json
import re
from Bio import SeqIO
from collections import Counter

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mwhiteside@canada.ca"



class SeqDict(object):
    """Library of non-redundant set of sequences for quick lookup



    """

    def __init__(self):
        """Constructor

        Args:
            None

        """

        self._seqdict = {}
        self._genenum = 0
        self._hash_algorithm = 'md5'


    @property
    def seqs(self):
        """dict: ."""

        return self._seqdict

    @property
    def num(self):
        """Number sequences"""

        return self._genenum

   
    def build(self, fasta_file, subtype_file):
        """Populate SeqDict with entries in fasta and subtype files

        Args:
            fasta_file (str): Filepath to input sequence file
            subtype_file (str): Filepath to input subtype

        Returns:
            None

        """

        # Check subtype file
        subtypes = {}

        for row in csv.reader(open(subtype_file,'r'),delimiter='\t'):
            name = row[0].upper()
            subt = row[1].lower()

            subtypes[name] = subt
            
        
        fasta = SeqIO.parse(fasta_file, 'fasta')
        for record in fasta:
            name = record.id.upper()
            seq = str(record.seq).upper()
            seq.replace('-','')
            this_subt = subtypes[name]

            matched = self.find(seq)
            if matched:
                # Existing identical sequence
                if not matched['subtype'] == this_subt:
                    raise Exception('Identical sequence {} has different subtype assignment than {} (expected: {}, assigned: {})'.format(
                        name, ','.join(matched['accessions']), matched['subtype'], this_subt))

                # Record gene/genome
                matched['accessions'].append(name)

            else:
                self.add(seq, name, this_subt)

        fasta.close()


    def find(self, seq):
        """Search for identical sequence in dictionary
        
        Args:
            seq (str): sequence

        Returns:
            None if no identical sequence found -or- SeqDict entry
            
        """

        
        searchstr = self.digest(seq)

        if searchstr in self.seqs:
            hits = self.seqs[searchstr]

            if seq in hits:
                return hits[seq]
            else:
                return None
            
        else:
            return None

   
    def add(self, seq, name, subt):
        """Add new entry to SecDict
        
        Args:
            seq (str): sequence
            name (str): accession
            subt (str): subtype

        Returns:
            None
            
        """

        searchstr = self.digest(seq)

        if not searchstr in self.seqs:
            self.seqs[searchstr] = {}

        self._genenum += 1

        self.seqs[searchstr][seq] = {
            'subtype': subt,
            'name': self.format_name(name, subt),
            'accessions': [name]
        }


    def format_name(self, name, subtype):
        """Standard gene naming"""

        # Try to isolate NCBI accession
        nm = re.sub(r'^(?:\d+\-)?([A-Z0-9\.]+)_.+$', r'\1', name.upper())
        genename = '{}-{}__{}'.format(self._genenum, nm, subtype.lower())

        return genename



    def store(self, reffile):
        """Save SeqDict to file"""
        with open(reffile, 'w') as rfh:
            json.dump(self.seqs,rfh)

        return None


    def load(self, reffile):
        """Load SeqDict from file"""
        with open(reffile, 'r') as rfh:
            self._seqdict = json.load(rfh)

        # Check format    
        for seqkey,seqs in self._seqdict.iteritems():
            for seq,seqentry in seqs.iteritems():
                self._genenum += 1

                for k in ['name','subtype','accessions']:
                    if not k in seqentry:
                        raise Exception('Improperly formated SeqDict object')

        return None


    def digest(self, seq):
        """Compute quick lookup hash value from sequence"""

        h = hashlib.new(self._hash_algorithm)
        h.update(seq)
        dig = h.hexdigest()

        return dig


    def write(self, fasta_filepath, subtype_filepath):
        """Write unique sequences in fasta format
        Args:
            fasta_filepath (str): Filepath to output fasta file
            subtype_filepath (str): Filepath to output subtype file
            
        Returns:
            None
        """

        with open(fasta_filepath, 'w') as f, open(subtype_filepath, 'w') as s:
            for seqkey,seqs in self.seqs.iteritems():
                for seq,seqentry in seqs.iteritems():

                    f.write(">{}\n{}\n".format(seqentry['name'], seq))
                    s.write("{}\t{}\n".format(seqentry['name'], seqentry['subtype']))


        return None


    def subtype_occurences(self):
        """Count occurences of each subtype in SeqDict

        Returns:
            Counter object

        """

        subtype_counts = Counter()

        for seqkey,seqs in self.seqs.iteritems():
            for seq,seqentry in seqs.iteritems():

                subtype_counts[seqentry['subtype']] += 1


    
        return subtype_counts



        


