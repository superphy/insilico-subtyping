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
from collections import Counter, defaultdict

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mwhiteside@canada.ca"



class SeqDict(object):
    """Library of non-redundant set of sequences for quick lookup



    """

    def __init__(self, n_loci=1):
        """Constructor

        Args:
            None

        """

        self._seqdict = {}
        self._genenum = 0
        self._hash_algorithm = 'md5'
        self._nloci = n_loci


    @property
    def seqs(self):
        """dict: ."""

        return self._seqdict

    @property
    def num(self):
        """Number sequences"""

        return self._genenum

   
    def build(self, fasta_files, subtype_file):
        """Populate SeqDict with entries in fasta and subtype files

        Args:
            fasta_file (str|list): Filepath to input sequence file or list of such files
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
            
        
        if isinstance(fasta_files, str):
            # Create list
            fasta_files = [fasta_files]

        if len(fasta_files) != self._nloci:
            raise Exception("Missing fasta file. {} fasta files provided for {} number of loci.".format(len(fasta_files), self._nloci))

        concat = LociConcat()
        sequences = concat.collect(fasta_files)
       
        for name,seqs in sequences.iteritems():
            name = name.upper()
            seq = ''.join(seqs).upper()
            seq.replace('-','')
            this_subt = subtypes[name]

            


    def find(self, seq):
        """Search for identical sequence in dictionary
        
        Args:
            seq (str|list): sequence or list of sequences

        Returns:
            None if no identical sequence found -or- SeqDict entry
            
        """

        if not isinstance(seq, str):
            # Concatenate list
            seq = ''.join(seq)
        
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
            seq (list|str): list of sequences or single sequence if one loci
            name (str): accession
            subt (str): subtype

        Returns:
            None
            
        """
        keystr = seq
        if self._nloci > 1:
            if isinstance(seq, str):
                raise Exception('Invalid seq parameter. Multiple loci need to be passed as list')
            keystr = ''.join(seq)

        else:
            if not isinstance(seq, str):
                keystr = seq[0]

        matched = self.find(seq)
        if matched:
            # Existing identical sequence
            # Append accession
            if not matched['subtype'] == subt:
                raise Exception('Identical sequence {} has different subtype assignment than {} (expected: {}, assigned: {})'.format(
                    name, ','.join(matched['accessions']), matched['subtype'], subt))

            # Record gene/genome
            matched['accessions'].append(name)

        else:
            # Add new sequence
            searchstr = self.digest(keystr)
            if not searchstr in self.seqs:
                self.seqs[searchstr] = {}

            self._genenum += 1

            self.seqs[searchstr][keystr] = {
                'subtype': subt,
                'name': self.format_name(name, subt),
                'accessions': [name]
            }

            if self._nloci > 1:
                self.seqs[searchstr][keystr]['loci'] = seq


    def accession_map(self, name):
        """Return dict of SeqDict names mapped to lists of original accessions"""

        acc_map = {}
        for seqkey,seqs in self._seqdict.iteritems():
            for seq,seqentry in seqs.iteritems():
                acc_map[seqentry['name']] = seqentry['accessions']

        acc_map


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
        keys =  ['name','subtype','accessions']
        if self._nloci > 1:
            keys.append('loci')
        for seqkey,seqs in self._seqdict.iteritems():
            for seq,seqentry in seqs.iteritems():
                self._genenum += 1

                for k in keys:
                    if not k in seqentry:
                        raise Exception('Improperly formated SeqDict object')

                if self._nloci > 1:
                    if len(seqentry['loci']) != self._nloci:
                        raise Exception('Improperly formated SeqDict object')

        return None


    def digest(self, seq):
        """Compute quick lookup hash value from sequence"""

        h = hashlib.new(self._hash_algorithm)
        h.update(seq)
        dig = h.hexdigest()

        return dig


    def write(self, fasta_filepaths, subtype_filepath):
        """Write unique sequences in fasta format
        Args:
            fasta_filepaths (list): List of filepaths to output fasta files
            subtype_filepath (str): Filepath to output subtype file
            
        Returns:
            None
        """

        if self._nloci == 1 and isinstance(fasta_filepaths, str):
            fasta_filepaths = [fasta_filepaths]

        if len(fasta_filepaths) != self._nloci:
            raise Exception("Missing fasta file. {} fasta files provided for {} loci.".format(len(fasta_filepaths), self._nloci))


        # Open all needed files
        try:
            fhs = [ open(filename, 'w') for filename in fasta_filepaths ]

            with open(subtype_filepath, 'w') as s:
                for seqkey,seqs in self.seqs.iteritems():
                    for seq,seqentry in seqs.iteritems():

                        if self._nloci == 1:
                            f = fhs[0]
                            f.write(">{}\n{}\n".format(seqentry['name'], seq))
                            
                        else:
                            for f,loci_seq in zip(fhs, seqentry['loci']):
                                f.write(">{}\n{}\n".format(seqentry['name'], loci_seq))
                               
                    s.write("{}\t{}\n".format(seqentry['name'], seqentry['subtype']))


        except IOError as e:
            for f in fhs:
                f.close()

            raise e


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



class LociConcat(object):
    """Concatenate loci from multiple fasta input files into supersequences.

    Occaisonly multiple loci are used to predict traits.  These can be concatenated
    to form supersequences.

    """
    
    def collapse(self, inputs, fasta_filepath=None):
        """Concatenate loci from mulitple fasta files

        Returns dict indexed by fasta ID containg supersequences.
        Optionally, writes to fasta file.

        Args:
            inputs (list): Filepaths to fasta sequences.
            fasta_filepath (str)[Optional]: If provided, write supersequences to this file
        
        Returns:
            dict

        Raise:
            Exception when missing loci in a file

        """

        i = 0
        uniq = Counter()
        sequences = defaultdict(str)

        for ff in inputs:
            fasta = SeqIO.parse(ff, 'fasta')
            for record in fasta:
                name = record.id

                uniq[name] += 1
                sequences[name] += str(record.seq)


            i += 1
            for name in uniq:
                if uniq[name] != i:
                    raise Exception('Missing gene {} in file {}'.format(name, ff))


        if fasta_filepath:
            with open(fasta_filepath, 'w') as f:
                for name,seq in sequences.iteritems():
                    f.write(">{}\n{}\n".format(name, seq))
                
        return sequences


    def collect(self, inputs):
        """Concatenate loci from mulitple fasta files

        Returns dict indexed by fasta ID containg list of loci sequences
       
        Args:
            inputs (list): Filepaths to fasta sequences.
        
        Returns:
            dict

        Raise:
            Exception when missing loci in a file

        """

        i = 0
        uniq = Counter()
        sequences = defaultdict(list)

        for ff in inputs:
            fasta = SeqIO.parse(ff, 'fasta')
            for record in fasta:
                name = record.id

                uniq[name] += 1
                sequences[name].append(str(record.seq))


            i += 1
            for name in uniq:
                if uniq[name] != i:
                    raise Exception('Missing gene {} in file {}'.format(name, ff))
 
        return sequences




