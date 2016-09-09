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
    test
    """

    def __init(self):
        None

class SubtypeParser(object):
    """Extract subtypes assignments from Genbank records

        Needs apriori list of valid subtype regular expressions and searches for
        instances of these in definition, keyword fields of the top-level record, 
        and in CDS feature fields allele, gene, product, note

    """

    def __init__(self, valid_subtypes, 
        record_fields=['definition','keyword'],
        feature_fields=['allele','gene','product','note']):
        """Constructor
        
            Args:
        
        """

        self._subtype_regexp = valid_subtypes
        self._feature_fields = feature_fields
        self._record_fields = record_fields



    def search_feature(self, feature):
        """Find instances of subtype in genbank feature

        """

        matches = []

        for f in feature_fields:
            if f in feature.qualifiers:

                for re in self._subtype_regexp:
                    matches.append(re.findall(feature.qualifiers[f]))
               

        return None


class DownloadUtils(object):
    """Use Biopython's Eutils to retrieve gene sequences from NCBI

    """

    def __init__(self, output_dir, organism, gene_names, subtype_parser, email="superphy.info@gmail.com"):
        """Constructor

            Args:
            
   

        """

        # Set email for eutils interface
        Entrez.email = email
       
        # Logger
        self.logger = logging.getLogger(__name__)

        # Output
        self._genes = gene_names
        filename = gene_names[0]
        self._outdir = output_dir
        self._gbfile = os.path.join(self._outdir, filename+'.gb')
        self._fastafile = os.path.join(self._outdir, filename+'.fasta')
        self._subtypefile = os.path.join(self._outdir, self._outdir, filename+'.txt')

        # Search string
        self._search_string = "(" + " OR ".join(("\"%s\"[Gene]" %x for x in gene_names)) + ") AND \""+organism+"\"[Organism]"
        
        self._db = "nucleotide"

        self._batch_size = 500

        # SubtypeParser object
        self._subtype_parser = subtype_parser



    @property
    def query_string(self):
        return self._search_string

    @property
    def output_directory(self):
        return self._outdir

    @property
    def subtype_parser(self):
        return self._subtype_parser


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


    def parse(self):
        """Extract sequence and subtype for valid genes in each genbank record

        """

        for gb_record in SeqIO.parse(open(self._gbfile,"r"), "genbank"):

            # Find the CDS matching the gene
            for (index, feature) in enumerate(gb_record.features):

                if feature.type == 'CDS':

                    # Check if this CDS is what we are looking for
                    matched = False
                    if 'gene' in feature.qualifiers and feature.qualifiers['gene'][0] in self._genes:
                        matched = True
                    elif 'product' in feature.qualifiers and feature.qualifiers['product'][0] in self._genes:
                        matched = True

                    # Found
                    if matched:

                        # Get sequence
                        dnaseq = feature.extract(gb_record.seq)
                        protseq = None

                        if feature.qualifiers['translation']:
                            protseq = feature.qualifiers['translation'][0]
                        else:
                            protseq = dnaseq.translate(table=11, to_stop=True)

                        dlen = len(dnaseq)
                        plen = len(protseq)

                        # Try to find subtype
                        subtype = self.subtype_parser.search_feature(feature)




