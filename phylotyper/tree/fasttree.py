#!/usr/bin/env python

"""Phylogenetic tree utilities.

Tools for building and manipulating phylogenetic trees

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


class FastTreeWrapper(object):
    """Run FastTree ML Tree program.

    Wrapper around the FastTree maximum likelihood phylogenetic tree building
    program


    """

    def __init__(self, config):
        """Constructor

        Args:
            config (PhylotyperConfig): Instance of PhylotyperConfig object

        """

        self._fasttree = config.get('external', 'fasttree')
        self._fasttree_args = {
            'nt': '-gtr -nt',
            'nt_fast': '-gtr -nt -nosupport -fastest -mlnni 4',
            'aa': '-gtr'
        }


    @property
    def fasttree(self):
        """str: Path to FastTree executable."""

        return self._fasttree

   
    def build(self, alignment_file, tree_file, nt=False, fast=False):
        """Run fasttree


        Args:
            alignment_file (str): Filepath to input alignment file
            tree_file (str): Filepath for output from fasttree
            ft (bool[Default=False]): Set to true for nucleotide sequences
            fast (bool[Default=False]): Set to true to run fasttree in mode optimized for speed

        Returns:
            None

        """


        cmd_args = None
        if nt and not fast:
            cmd_args = self._fasttree_args['nt']
        elif nt and fast:
            cmd_args = self._fasttree_args['nt_fast']
        elif not nt and not fast:
            cmd_args = self._fasttree_args['aa']


        if not cmd_args:
            raise Exception("No fasttree mode defined for argument settings nt: {}, fast: {}".format(nt, fast))

        cmd = "{} {} {} > {}".format(self.fasttree, cmd_args, alignment_file, tree_file)

        try:
            check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
        except CalledProcessError as e:
            msg = "FastTree failed: {} (return code: {}).".format(e.output, e.returncode)                                                                                                   
            raise Exception(msg)

        None


