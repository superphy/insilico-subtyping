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
from Bio import SeqIO
from collections import Counter

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mwhiteside@canada.ca"


class SeqDict(object):
    """Generate non-redundant set of sequences



    """

    def __init__(self):
        """Constructor

        Args:
            None

        """

        self._seqdict = {}
        self._namedict = {}


    @property
    def seqs(self):
        """dict: ."""

        return self._seqdict

    @property
    def names(self):
        """dict: ."""

        return self._namedict

   
    def load(self, fasta_file, subtype_file):
        """Populate SeqDict


        Args:
            fasta_file (str): Filepath to input sequence file
            subtype_file (str): Filepath to input subtype

        Returns:
            None

        """

        # Check subtype file
        alleles = Counter()
        subtypes = {}
        assigned = {}

        for row in csv.reader(open(subtype_file,'r'),delimiter='\t'):
            name = row[0]
            subt = row[1]

            subtypes[name] = subt
            
        
        fasta = SeqIO.parse(fasta_file, 'fasta')
        for record in fasta:
            name = record.id
            seq = str(record.seq)

            this_subt = subtypes[name]

            if seq in self.seqs:
                if not self.seqs[seq] == this_subt:
                    raise Exception('Identical sequence {} has different subtype assignment than {} (expected: {}, assigned: {})'.format(
                        name, ','.join(assigned[seq]), self.seqs[seq], this_subt))

                # Have multiple copies of sequence, generate new meta-name based on subtype
                alleles[this_subt] += 1
                self.names[seq] = '{}_allele{}'.format(this_subt, alleles[this_subt])

            else:
                self.seqs[seq] = this_subt
                self.names[seq] = name

            # track ids with this sequence
            if seq in assigned:
                assigned[seq].append(name)
            else:
                assigned[seq] = [name]

        fasta.close()


    def find(self, sequence):
        """Find matching sequence in SeqDict

        Args:
            sequence (str): Bio sequence
            

        Returns:
            name (str) or None

        """

        if sequence in self.seqs:
            return self.seqs[sequence]
        else:
            return None


    def write(self, fasta_filepath, subtype_filepath):
        """Write unique sequences in fasta format

        Args:
            fasta_filepath (str): Filepath to output fasta file
            subtype_filepath (str): Filepath to output subtype file
            

        Returns:
            None

        """

        with open(fasta_filepath, 'w') as f, open(subtype_filepath, 'w') as s:
            for seq,subt in self.seqs.iteritems():
                name = self.names[seq]

                f.write(">{}\n{}\n".format(name, seq))
                s.write("{}\t{}\n".format(name, subt))






    




        


