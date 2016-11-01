#!/usr/bin/env python

"""Loci Detection

Given set of reference genes, find loci in input genome sequences

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

import logging
import re
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from os.path import splitext, basename
from subprocess import check_output, CalledProcessError, STDOUT


class LociSearch(object):
    """Blast loci sequence database against genome sequences


    """
    
    def __init__(self, config, database, inputs):
        """Constructor

        Either a pre-existing blast database or one or more loci fasta files are required as arguments.
        Each input fasta file must contain one or more examples of a single loci. Multiple loci are input
        as multiple fasta files.

        Args:
            config (PhylotyperConfig): Instance of PhylotyperConfig object
            database (str): Blast database location
            inputs (tuple): Filepaths to fasta sequences.

        Returns:
            None

        """

        self.logger = logging.getLogger('phylotyper.genome.loci.LociSearch')

        self._db_title = 'Phylotyper Blast Database v1'
        self._db_prefix = 'ptloci'
        self._alignment_coverage = 0.95

        # Default options
        self._blast_options = {
            'perc_identity': 90,
            'evalue': 0.00001,
            'db': database,
            'outfmt': 5
        }

        self._makeblastdb_options = {
            'input_type': 'fasta',
            'dbtype': 'nucl',
            'out': database,
            'title': self._db_title
        }

        self._blastdbcmd_options = {
            'db': database
        }

        self._blastdbcmd_options2 = {
            'db': database,
            'outfmt': '%t\t%l',
            'entry': 'all'
        }

        self._makeblastdb_exe = config.get('external', 'makeblastdb')
        self._blastn_exe = config.get('external', 'blastn')
        self._blastdbcmd_exe = config.get('external', 'blastdbcmd')


        if inputs:
            self.make_db(inputs)
        
        self.load_db(database)
      

    # Prepare
    def make_db(self, inputs):
        """Make blast database

        Args:
            database (str): Blast database location
            inputs (tuple): Filepaths to fasta sequences.

        Returns:
            None

        """

        # Create naming system that identifies loci sets in the database
        with tempfile.NamedTemporaryFile() as tmpfh:

            # Concatenate all the fasta sequences
            # Assign a loci set identifier
            loci_id = 1
            for ff in inputs:
                fasta = SeqIO.parse(ff, 'fasta')
                for record in fasta:
                    newheader = '{}{}|{}'.format(self._db_prefix, loci_id,record.description)
                    tmpfh.write(">{}\n{}\n".format(newheader, record.seq))

                loci_id += 1

            tmpfh.flush()

            opts = self._makeblastdb_options
            opts['in'] = tmpfh.name

            cmd = "{blastcmd} -in {in} -out {out} -dbtype {dbtype} -input_type {input_type} -title \"{title}\"".format(
                blastcmd=self._makeblastdb_exe, **opts)

            self.logger.debug('Creating new blast database with input files: {}'.format(','.join(inputs)))
            self.logger.debug('Running Blast command-line command:\n{}'.format(cmd))

            try:
                check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
            except CalledProcessError as e:
                msg = "LociSearch failed: {} (return code: {}).".format(e.output, e.returncode)                                                                                                   
                raise Exception(msg)

        
    def load_db(self, database):
        """Check that this blast database was created by LociSearch

        Args:
            database (str): Blast database location

        Raises Exception

        """

        opts = self._blastdbcmd_options
        cmd = "{blastcmd} -info -db {db}".format(blastcmd=self._blastdbcmd_exe, **opts)

        self.logger.debug('Checking blast database')
        self.logger.debug('Running Blast command-line command:\n{}'.format(cmd))

        try:
            info = check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
        except CalledProcessError as e:
            msg = "LociSearch failed: {} (return code: {}).".format(e.output, e.returncode)                                                                                                   
            raise Exception(msg)

        database_header = r'Database: {}'.format(self._db_title)
        if not re.match(database_header, info):
            msg = "Not a Phylotyper Blast database: {}".format(database)                                                                                                   
            raise Exception(msg)
        else:
            self.logger.debug('Verified Phylotyper Blast database')

        # Load subject lengths
        opts = self._blastdbcmd_options2
        cmd = "{blastcmd} -db {db} -entry {entry} -outfmt \"{outfmt}\"".format(blastcmd=self._blastdbcmd_exe, **opts)

        self.logger.debug('Loading sequence lengths')
        self.logger.debug('Running Blast command-line command:\n{}'.format(cmd))

        try:
            lengths = check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
        except CalledProcessError as e:
            msg = "LociSearch failed: {} (return code: {}).".format(e.output, e.returncode)                                                                                                   
            raise Exception(msg)

        self._sequence_lengths = {}
        for line in lengths.splitlines():
            [desc, length] = line.split('\t',2)
            
            if not re.match(r'^{}'.format(self._db_prefix), desc):
                raise Exception('Improperly formatted sequence title in database: {}'.format(desc))

            self._sequence_lengths[desc] = int(length)


        None


    def search(self, genome, output, append=True, fasta_prefix=None):
        """Run Blast with genome as input

        Args:
            genome (str): Filepath to input fasta sequence
            output (str): Filepath for fasta output
            append (bool): False = overwrite output file
            fasta_prefix (str): All fasta headers in output will start with this

        Raises Exception

        """

        def locad(contig, start, stop):
            # Loci address
            s = sorted([start, stop])
            return '{}:{}-{}'.format(contig, *s)

        loci = {}

        if not fasta_prefix:
            # If an output fasta header ID prefix is not provided
            # try to identify in the headers from the input fasta file itself
            fasta = SeqIO.parse(genome, 'fasta')
            record = fasta.next()
            result = re.search(r'^(\w{3}\|\S+\|)', record.id)

            if result:
                fasta_prefix = result.group(1)

            else:
                # If the fasta ID does not match convention, use the filename
                fp = splitext(basename(genome))[0]
                fasta_prefix = 'lcl|{}|'.format(fp)

        self.logger.debug('Output fasta headers will use <{}> as the prefix'.format(fasta_prefix))

        with tempfile.NamedTemporaryFile() as tmpfh:
            opts = self._blast_options
            opts['query'] = genome
            opts['out'] = tmpfh.name

            cmd = "{blastcmd} -perc_identity {perc_identity} -evalue {evalue} -outfmt {outfmt} -db {db} -query {query} -out {out}".format(
                blastcmd=self._blastn_exe, **opts)

            self.logger.debug('Running Blast')
            self.logger.debug('Running Blast command-line command:\n{}'.format(cmd))
            try:
                check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
            except CalledProcessError as e:
                msg = "LociSearch failed: {} (return code: {}).".format(e.output, e.returncode)                                                                                                   
                raise Exception(msg)

            self.logger.debug('Parsing Blast output')
            
            with open(tmpfh.name, 'r') as blast_handle:
                blast_records = NCBIXML.parse(blast_handle)

                # Iterate through results
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                      
                        result = re.search(r' (?P<subject>{}(?P<lociset>\d+)\|.+$)'.format(self._db_prefix), alignment.title)
                        subject = result.group('subject')
                        loci_set = result.group('lociset')
                        if not loci_set in loci:
                            loci[loci_set] = {}

                        for hsp in alignment.hsps:
                            
                            # Is majority of subject aligned?
                            alen = hsp.align_length
                            slen = self._sequence_lengths[subject]
                            coverage = alen / slen

                            if coverage > self._alignment_coverage:
                                # Found one, save it if its a new loci

                                addr = locad(blast_record.query, hsp.query_start, hsp.query_end)
                                # self.logger.debug('Hit {} location: {}, length: {}, coverage: {}, identities:{}'.format(
                                #     subject, addr, alen, int(coverage*100), hsp.identities))

                                if not addr in loci[loci_set]: 
                                    # New
                                    seq = Seq(hsp.query, IUPAC.unambiguous_dna)
                                    if hsp.sbjct_start > hsp.sbjct_end:
                                        # For some reason strandedness is missing
                                        seq = seq.reverse_complement()
                                    loci[loci_set][addr] = seq
                
                for ls in loci:
                    self.logger.debug("{} hits found for loci {}".format(len(loci[ls]), ls))

                # Concatenate loci sets
                loci_sequences = self.concatenate(loci.keys(), loci)

                # Output
                mode = 'w'
                if append:
                    mode = 'a'

                nloci = 1
                with open(output, mode) as outfh:
                    for s in loci_sequences:
                        outfh.write("\n>{}loci{}\n{}".format(fasta_prefix, nloci, s))


    def concatenate(self, loci_sets, loci):
        # Generate all combinations of loci set hits

        ls = loci_sets.pop()

        if not loci_sets:
            return loci[ls].values()

        else:
            seqs = self.concatenate(loci_sets, loci)
            newseqs = []
            for v in loci[ls].values():
                for s in seqs:
                    newseqs.append(v+s)

            return newseqs






