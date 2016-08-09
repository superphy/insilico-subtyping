#!/usr/bin/env python

"""Multiple Sequence Alignment utilities.

Tools for building and manipulating alignments

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

from subprocess import check_output, CalledProcessError, STDOUT

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mwhiteside@canada.ca"


class SeqAligner(object):
    """Run MSA program.

    Wrapper around the multiple sequence alignment program

    """

    def __init__(self, config):
        """Constructor

        Args:
            config (PhylotyperConfig): Instance of PhylotyperConfig object

        """

        self._mafft = config.get('external', 'mafft')
        self._mafft_args = {
            'auto': '--auto',
            'add': '--add'
        }


    @property
    def aligner(self):
        """str: Path to executable."""

        return self._mafft

    @property
    def aligner_args(self):
        """str: Path to executable."""

        return self._mafft_args['auto']

    @property
    def aligner_add_args(self):
        """str: Path to executable."""

        return self._mafft_args['add']

   
    def align(self, fasta_file, alignment_file):
        """Run aligner program

        Args:
            fasta_file (str): Filepath to input fasta file
            alignment_file (str): Filepath for output from alignment


        Returns:
            None

        """

        cmd_args = self.aligner_args
        cmd = "{} {} {} > {}".format(self.aligner, cmd_args, fasta_file, alignment_file)

        try:
            check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
        except CalledProcessError as e:
            msg = "Aligner failed: {} (return code: {}).".format(e.output, e.returncode)                                                                                                   
            raise Exception(msg)

        None

    def add(self, fasta_file, existing_alignment_file, output_alignment_file):
        """Add new sequences to existing alignment

        Args:
            fasta_file (str): Filepath to input fasta file
            existing_alignment_file (str): Filepath to aliged fasta file
            output_alignment_file (str): Filepath for output from alignment

        Returns:
            None

        """

        cmd_args = self.aligner_add_args
        cmd = "{} {} {} {} > {}".format(self.aligner, cmd_args, fasta_file, existing_alignment_file, 
            output_alignment_file)

        try:
            check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
        except CalledProcessError as e:
            msg = "Aligner --add failed: {} (return code: {}).".format(e.output, e.returncode)                                                                                                   
            raise Exception(msg)

        None


