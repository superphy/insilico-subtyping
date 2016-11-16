#!/usr/bin/env python

"""Download PMIDs for references related to pubmed journals


Example:
        $ python download_pmids.py outputfile

"""

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import re
import logging

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mdwhitesi@gmail.com"


class Blast(object):
    """Retrieve sequences using BLAST

    """

    def __init__(self, organism, dna=True):
        """Constructor
        
            Args:
        
        """

        if not dna:
            raise Exception('Protein blast options not defined')

        self._entrez_query = '{}[Organism]'.format(
            "\""+organism+"\""
        )

        self._qblast_args = {
            'program': 'blastn',
            'database': 'nt',
            'megablast': True,
            'expect': 1e-10,
            'perc_ident': 85,
            'alignments': 100,
            'descriptions': 100,
            'hitlist_size': 50,
            'entrez_query': self._entrez_query
        }

        self._query_coverage = float(0.95)

        self.logger = logging.getLogger(__name__)


    def run(self, inputfile, outputfile):
        """Perform blast job, one sequence at a time

        """
        args = self._qblast_args

        with open(outputfile, 'w') as fh:

            fasta_iterator = SeqIO.parse(open(inputfile, "r"), "fasta")

            for fasta_record in fasta_iterator:

                naln = 0
                nhits = 0
                qlen = float(len(fasta_record.seq))
                result_handle = NCBIWWW.qblast(sequence=fasta_record.seq, **args)
                blast_records = NCBIXML.parse(result_handle)

                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        # Only look at the best alignment
                        # Assumes hsps are sorted by expect value
                        hsp = alignment.hsps[0]

                        # Is majority of query aligned?
                        alen = hsp.query_end - hsp.query_start + 1
                        self.logger.debug('Alignment length: {}, coverage: {}, idents: {}'.format(alen, (alen / qlen), hsp.identities))
                        if (alen / qlen) > self._query_coverage:

                            # Is it identical? I want distinct sequences
                            if hsp.identities < alen:
                                sheader = alignment.title
                                sseq = hsp.sbjct

                                ids = self.parse_header(sheader)
                                if ids:
                                    fh.write(">{}\n{}\n".format(ids[1], sseq))
                                    nhits += 1
                                else:
                                    self.logger.info('Unrecognized blast header: {}'.format(sheader))

                            else:
                                self.logger.info('Identical sequences')
                        else:
                            self.logger.info('Alignment too short')

                        naln += 1

                self.logger.debug("{} sequences found for query {} out of {} alignments".format(nhits, fasta_record.id, naln))


                        
    def parse_header(self, header):
        """Extract genbank ID from header (if it exists)

        """

        match = re.search(r'^gi\|\w+\|(\w+)\|([\w\.]+)\|\s', header)

        if match:
            return(match.groups())
        else:
            return None



   