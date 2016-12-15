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

from __future__ import division

import logging
import re
import tempfile
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from os.path import splitext, basename
from subprocess import check_output, CalledProcessError, STDOUT


class LociSearch(object):
    """Blast loci sequence database against genome sequences


    """
    
    def __init__(self, config, database, inputs, sequence_type):
        """Constructor

        Either a pre-existing blast database or one or more loci fasta files are required as arguments.
        Each input fasta file must contain one or more examples of a single loci. Multiple loci are input
        as multiple fasta files.

        Args:
            config (PhylotyperConfig): Instance of PhylotyperConfig object
            database (str): Blast database location
            inputs (list): Filepaths to fasta sequences.
            sequence_type (str): nucl or prot

        Returns:
            None

        """

        self.logger = logging.getLogger('phylotyper.genome.loci.LociSearch')

        if sequence_type == 'nucl' or sequence_type == 'prot':
            self.seqtype = sequence_type
        else:
            raise Exception('Invalid sequence_type argument {}'.format(sequence_type))

        self._db_title = 'Phylotyper Blast Database v1'
        self._db_prefix = 'ptloci'
        self._alignment_coverage = 0.95
        self._percent_identity = 0.90

        if self.seqtype == 'prot':
            self._percent_identity = 0.80

        self._database = database

        # Default options
        self._blast_options = {
            'evalue': 0.00001,
            #'evalue': 1,
            'db': database,
            'outfmt': 5
        }

        filt = 'dust' if self.seqtype == 'nucl' else 'seg'
        self._blast_specific_options = {
            filt: 'no'
        }

        self._makeblastdb_options = {
            'input_type': 'fasta',
            'dbtype': self.seqtype,
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
        blastexe = 'blastn' if self.seqtype == 'nucl' else 'blastx'
        self._blast_exe = config.get('external', blastexe)
        self._blastdbcmd_exe = config.get('external', 'blastdbcmd')


        if inputs:
            self.make_db(inputs)
        
        self.load_db(database)
      

    # Prepare
    def make_db(self, inputs):
        """Make blast database

        Args:
            inputs (tuple): Filepaths to fasta sequences.

        Returns:
            None

        """

        # Create naming system that identifies loci sets in the database
        with tempfile.NamedTemporaryFile() as tmpfh:

            # Assign a loci set identifier
            loci_id = 1
            for ff in inputs:
                fasta = SeqIO.parse(ff, 'fasta')
                for record in fasta:
                    newheader = '{}{}|{}'.format(self._db_prefix, loci_id, record.description)
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

        Raises Exception if blast search fails

        """

        loci = {}
        locations = defaultdict(list)

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

            cmd = "{blastcmd} -evalue {evalue} -outfmt {outfmt} -db {db} -query {query} -out {out}".format(
                blastcmd=self._blast_exe, **opts)
            if self._blast_specific_options:
                cmd += ' '+' '.join(['-'+k+' '+v for k,v in self._blast_specific_options.items()])

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

                            alen = hsp.align_length

                            # Is the precent identity above thresold
                            # Note: I want identically functional proteins, not the broader ortholog/homolog families
                            pident = hsp.identities / alen
                            
                            # Is majority of subject aligned?
                            slen = self._sequence_lengths[subject]
                            coverage = alen / slen

                            addr = self.locad(blast_record.query, hsp.query_start, hsp.query_end)
                            self.logger.debug('Hit {} location: {}, length: {}, subject length:{}, coverage: {}, identities:{}, percent identity:{}'.format(
                                subject, addr, alen, slen, int(coverage*100), hsp.identities, int(pident*100)))

                            if coverage > self._alignment_coverage and pident > self._percent_identity:
                                # Found one, save it if its a new loci

                                self.logger.debug('Found loci {}'.format(addr))

                                if not addr in loci[loci_set]:

                                    # New
                                    seq = hsp.query
                                    if self.seqtype == 'nucl':
                                        seq = Seq(hsp.query, IUPAC.unambiguous_dna)
                                        if hsp.sbjct_start > hsp.sbjct_end:
                                            # For some reason strandedness is missing
                                            seq = seq.reverse_complement()
                                    else:
                                        nstops = seq.count('*')
                                        if nstops > 1:
                                            self.logger.warning('Query hit {} has multiple stop codons: {}'.format(blast_record.query, nstops))
                                        seq = str(seq).replace('*','')

                                    loci[loci_set][addr] = seq
                                    locations[blast_record.query].append(sorted((hsp.query_start, hsp.query_end)))
                
                for ls in loci:
                    loci_hits = ','.join(loci[ls].keys())
                    self.logger.debug("{} hits found for loci {}\n\t({})".format(len(loci[ls]), ls, loci_hits))

                # Find non-overlapping hits
                nonoverlapping = self.non_overlapping(locations)

                # Remove overlapping hits
                new_loci = {}
                for ls in loci.keys():
                    new_loci[ls] =  { addr: loci[ls][addr] for addr in loci[ls].keys() if addr in nonoverlapping }

                for ls in new_loci:
                    loci_hits = ','.join(new_loci[ls].keys())
                    self.logger.info("{} hits found for loci {} after filtering overlapping hits\n\t({})".format(len(new_loci[ls]), ls, loci_hits))

                # Concatenate loci sets
                loci_sequences = self.concatenate(new_loci.keys(), new_loci)

                # Output
                mode = 'w'
                if append:
                    mode = 'a'

                nloci = 1
                with open(output, mode) as outfh:
                    for s in loci_sequences:
                        outfh.write("\n>{}loci{} {}\n{}".format(fasta_prefix, nloci, s[0], s[1]))
                        nloci += 1


    def concatenate(self, loci_sets, loci):
        # Generate all combinations of loci set hits

        ls = loci_sets.pop()

        if not loci_sets:
            return loci[ls].iteritems()

        else:
            seqs = self.concatenate(loci_sets, loci)
            newseqs = []
            for k,v in loci[ls].iteritems():
                for s in seqs:
                    newseqs.append((k+','+s[0],v+s[1]))

            return newseqs


    def locad(self, contig, start, stop):
            # Loci address
            s = sorted([start, stop])
            return '{}:{}-{}'.format(contig, *s)


    def non_overlapping(self, locations):
        # Remove hits that are contained within other hits

        def position_sort(x, y):
            if x[0] == y[0]:
                return y[1] - x[1]
            else:
                return x[0] - y[0]

        final = []

        for q in locations.keys():
        
            starts = sorted(locations[q], cmp=position_sort)

            # Save the first one - can't be overlapping
            s = 0
            addr = self.locad(q, starts[s][0], starts[s][1])
            final.append(addr)

            mx = len(starts)
            while s < (mx-1):
                # Find the next non-overlapping in the list 
                l = [i for i in range(s+1,mx) if starts[i][1] > starts[s][1]]
                if l:
                    s = l[0]
                    addr = self.locad(q, starts[s][0], starts[s][1])
                    final.append(addr)
                else:
                    # Exit loop, no more non-overlapping
                    s = mx
                
        return final


    def translate(self, loci):
        # Translate genome hits to amino acid sequences
        # Discards failed translations 

        for ls in loci.keys():
            for l in loci[ls].keys():
                sobj = loci[ls][l]
                try:
                    pobj = sobj.translate(table=11, stop_symbol='*')
                except Exception as e:
                    self.logger.warning('Translation failed for {}: {}'.format(l,e))

                nstops = pobj.count('*')

                if nstops > 1:
                    try:
                        sobj1 = sobj.reverse_complement()
                        pobj1 = sobj1.translate(table=11, stop_symbol='*')

                        nstops1 = pobj1.count('*')
                        if nstops1 > 1:
                            self.logger.warning('Frameshift error for {} ({} & {} stop codons in forward and reverse)'.format(l, nstops, nstops1))
                            if nstops1 < nstops:
                                sobj = sobj1
                                pobj = pobj1
                        else:
                            sobj = sobj1
                            pobj = pobj1

                    except Exception as e:
                        self.logger.warning('Translation failed for {}: {}'.format(l,e))


                loci[ls][l] = pobj



