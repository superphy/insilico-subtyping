"""NCBI Gene Sequence Download Utils

Classes for downloading genes from NCBI

Examples:
	To retrieve eae genes in E.coli:

        $ python main.py subtype ../phylotyper_example.ini ecoli_stx1 test/ecoli_stx1.ffn output/test/

"""


from Bio import Entrez, SeqIO
import logging
import os
import time

class GeneFilter(object):
	"""

	"""

	None

class SubtypeParser(object):
	"""Extract subtypes assignments from Genbank records

	Needs apriori list of valid subtype regular expressions and searches for
	instances of these in definition, keyword fields of the top-level record, 
	and in CDS feature fields allele, gene, product

	"""

	def __init__(self, valid_subtypes):
        """Constructor
        
        	Args:
            
   

        """

        self._subtype_regexp = valid_subtypes


    def find_subtype(self, input):
    	"""Find instances of subtype in input

    	"""

    	
	

class DownloadUtils(object):
    """Use Biopython's Eutils to retrieve gene sequences from NCBI

    """

    def __init__(self, output_dir, organism, gene_names, email="superphy.info@gmail.com"):
        """Constructor

        	Args:
            
   

        """

        # Set email for eutils interface
        Entrez.email = email
       
    	# Logger
        self.logger = logging.getLogger(__name__)

        # Output
        filename = gene_names[0]
        self._outdir = output_dir
        self._gbfile = os.path.join(self._outdir, filename+'.gb')
        self._fastafile = os.path.join(self._outdir, filename+'.fasta')
        self._subtypefile = os.path.join(self._outdir, self._outdir, filename+'.txt'

        # Search string
        gene_names = ("\"%s\"[Gene]" %x for x in gene_names)
        self._search_string = "(" + " OR ".join(gene_names) + ") AND \""+organism+"\"[Organism]"
        
        self._db = "nucleotide"

        self._batch_size = 500



    @property
    def query_string(self):
        return self._search_string

    @property
    def output_directory(self):
        return self._outdir


    def download(self):
    	"""Perform download of genes

    	"""

    	# Create directory to store genbank files
    	gbout = open(self._gbfile, 'w')

    	search = Entrez.esearch(db=self._db, term=self._search_string,
        	retmax=1, usehistory="y")

    	search_results = Entrez.read(search)
        search.close()

    	count = int(search_results["Count"])
    	webenv = search_results["WebEnv"]
    	query_key = search_results["QueryKey"]


        self.logger.info("%i genes sequences in NCBI matching query %s" % (count, self._search_string))

        # Download in batches
        self.logger.info("Starting download...")

        for start in range(0, count, self._batch_size):
            end = min(count, start+self._batch_size)
           
            attempt = 1
            while attempt <= 3:
                try:
                    fetch = Entrez.efetch(db=self._db, query_key=query_key,
                        retstart=start, retmax=self._batch_size, webenv=webenv, usehistory='y',
                        retmode='text', rettype='gb')

                    gbout.write(fetch.read())
                    fetch.close()

                    attempt = 4
                    
                except Exception, e:
                    self.logger.warning("Received error from server %s" % str(e))
                    self.logger.warning("Attempt %i of 3 for record %i to %i" % (attempt, start+1, end))
                    attempt += 1
                    time.sleep(15)

            self.logger.debug("%i to %i of %i retrieved" % (start+1, end, count))

        gbout.close()
        self.logger.info("download complete.")

        return True


    def parse(self,):
    	"""Extract sequence and subtype for valid genes in each genbank record

    	"""

    	for gb_record in SeqIO.parse(open(self._gbfile,"r"), "genbank"):


