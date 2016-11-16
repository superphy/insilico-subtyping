"""NCBI Gene Sequence Download Utils

Classes for downloading genes from NCBI

Examples:
    To retrieve eae genes in E.coli:

        $ python main.py subtype ../phylotyper_example.ini ecoli_stx1 test/ecoli_stx1.ffn output/test/

"""


from Bio import Entrez, SeqIO
from collections import Counter
import logging
import os
import re
import time


class GeneFilter(object):
    """Remove spurious genes based on user-supplied criteria

        User defines boolean test functions that examine
        sequence or genbank record. This class runs those tests

    """

    def __init__(self, sequence_tests=[], genbank_tests=[]):
        """Constructor
        
            Args:
        
        """

        for f in sequence_tests:
            if not callable(f):
                raise Exception('Attribute sequence_tests must contain a list of callable functions (failed: {}'.format(f))
        self._sequence_tests = sequence_tests

        for f in genbank_tests:
            if not callable(f):
                raise Exception('Attribute genbank_tests must contain a list of callable functions (failed: {}'.format(f))
        self._genbank_tests = genbank_tests


    def test_sequence(self, sequence):
        """Run test functions on sequence

            Args:
                sequence (str) 

        Returns boolean

        """

        for test in self._sequence_tests:
            if not test(sequence):
                return False

        return True

    def test_genbank(self, record):
        """Run test functions on BioPython genbank
        record object

            Args:
                record (Bio.Record) 

        Returns boolean

        """

        for test in self._sequence_tests:
            if not test(record):
                return False

        return True


class SubtypeParser(object):
    """Extract subtypes assignments from Genbank records

        Needs apriori list of valid subtype regular expressions and searches for
        instances of these in definition, keyword fields of the top-level record, 
        and in CDS feature fields allele, gene, product, note

    """

    def __init__(self, valid_subtypes,
        record_fields=['description'],
        feature_fields=['allele','gene','product','note'],
        source_fields=['organism','strain'],
        annotation_fields=['source','organism']):
        """Constructor
        
            Args:
        
        """

        self._subtype_regexp = valid_subtypes
        self.feature_fields = feature_fields
        self.record_fields = record_fields
        self.source_fields = source_fields
        self.annotation_fields = annotation_fields
        self._gotchas = [re.compile(r'ortholog')] # Indicate annotation transfer rather than direct evidence


    def search_feature(self, feature):
        """Find instances of subtype in genbank feature

        """

        matches = set([])

        for f in self.feature_fields:
            if f in feature.qualifiers:

                for regexp in self._subtype_regexp:
                    for v in feature.qualifiers[f]:
                        hits = regexp.findall(v)

                        if len(hits) > 0:
                            ok = True
                            for gotcha in self._gotchas:
                                if gotcha.search(v):
                                    ok = False
                                    break

                            if ok:
                                hits = self.format(hits)
                                for h in hits:
                                    matches.add(h)
               
        return matches


    def search_source(self, feature):
        """Find instances of subtype in genbank source feature

        """

        matches = set([])

        for f in self.source_fields:
            if f in feature.qualifiers:

                for regexp in self._subtype_regexp:
                    for v in feature.qualifiers[f]:
                        hits = regexp.findall(v)
                        hits = self.format(hits)
                        for h in hits:
                            matches.add(h)
               
        return matches
        

    def search_annotations(self, record):
        """Find instances of subtype in genbank annotations

        """

        matches = set([])

        for f in self.annotation_fields:
            if f in record.annotations:
                v = record.annotations[f]

                for regexp in self._subtype_regexp:
                    hits = regexp.findall(v)
                    hits = self.format(hits)
                    for h in hits:
                        matches.add(h)
               
        return matches


    def search_record(self, record):
        """Find instances of subtype in genbank record

        Returns set of matches returned by findall

        """

        matches = set([])

        for f in self.record_fields:
            if f in dir(record):

                for regexp in self._subtype_regexp:
                    hits = regexp.findall(getattr(record, f))
                    hits = self.format(hits)
                    for h in hits:
                        matches.add(h)
               
        return matches


    def format(self, stlist):
        """Apply standard formatting

        Formatting:
            1. to lowercase
            2. whitespace to '-'

        Returns list with formatted strings

        """

        fixed = []
        for s in stlist:
            s = s.lower()
            s = re.sub(r'\s', r'-', s)
            fixed.append(s)

        return fixed




class DownloadUtils(object):
    """Use Biopython's Eutils to retrieve gene sequences from NCBI

    """

    def __init__(self, output_dir, organism, gene_names, subtype_parser, gene_filter, gene_regex=None, dna=True, wgs=False, email="superphy.info@gmail.com"):
        """Constructor

            Args:
            
   

        """

        # Set email for eutils interface
        Entrez.email = email
       
        # Logger
        self.logger = logging.getLogger(__name__)

        if(gene_regex):
            self._genes = gene_regex
        else:
            self._genes = [re.compile(re.escape(x)) for x in gene_names]
        
        filename = gene_names[0]
        self._outdir = output_dir
        self._gbfile = os.path.join(self._outdir, filename+'.gb')
        self._tmpfile = os.path.join(self._outdir, 'tmp.gb')
        self._first_download = True
        self._fastafile = os.path.join(self._outdir, filename+'.fasta')
        self._subtypefile = os.path.join(self._outdir, filename+'.txt')
        self._dna_sequences = dna
        
        # Search string
        self._search_string = "{} AND {}[Organism]".format(
            "(" + " OR ".join(("\"%s\"[Gene]" %x for x in gene_names)) + " OR " +
            " OR ".join(("\"%s\"[Product]" %x for x in gene_names)) +
             ")",
            "\""+organism+"\""
        )

        if not wgs:
            self._search_string += " NOT WGS[Keyword]"
        
        if self._dna_sequences:
            self._db = "nuccore"
            self._dbfrom = "pubmed"
            self._linkname = "pubmed_nuccore"
        else:
            raise Exception('Protein not tested')


        self._batch_size = 500

        # SubtypeParser object
        self._subtype_parser = subtype_parser

        # GeneFilter object
        self._gene_filter = gene_filter



    @property
    def query_string(self):
        return self._search_string

    @property
    def output_directory(self):
        return self._outdir

    @property
    def subtype_parser(self):
        return self._subtype_parser

    @property
    def gene_filter(self):
        return self._gene_filter

    @property
    def subtypefile(self):
        return self._subtypefile

    @property
    def fastafile(self):
        return self._fastafile

    @property
    def genbankfile(self):
        return self._genbankfile


    def download_pubmed_genes(self, pmidfile=None, pmids=None):
        """Perform download of genes linked to pubmed IDs

        Genes must also match search query

            Args:
                pmidfile (str): Filepath containing one PMID ID per line
                    -or- (pmidfile will override pmids)
                pmids (list): List of pmids

            Returns True

        """

        self.prep_download()

        if pmidfile:
            with open(pmidfile, 'r') as f:
                pmids = f.read().splitlines()

        pmid_string = ','.join(pmids)
        links = Entrez.elink(db=self._db, dbfrom=self._dbfrom, linkname=self._linkname, id=pmid_string, cmd='neighbor_history')

        link_results = Entrez.read(links)
        links.close()

        # Only one database, so use first entry
        db_record = link_results[0]
        link_set = db_record['LinkSetDbHistory'][0]
        if not link_set['LinkName'] == self._linkname:
            raise Exception('Elink result LinkSetDB is not expected target: {}'.format(self._linkname))

        webenv = db_record["WebEnv"]
        query_key = link_set["QueryKey"]

        # Limit to genes matching search query
        search = Entrez.esearch(db=self._db, term=self._search_string, webenv=webenv, query_key=query_key, retmax=1, usehistory='y')

        search_results = Entrez.read(search)
        search.close()

        count = int(search_results["Count"])
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        self.logger.info("%i genes sequences in NCBI matching query %s found in pubmed records: %s" % (count, self._search_string,
            pmid_string))

        self.fetch(count, webenv, query_key)
       
        return True


    def download_genes(self):
        """Perform download of genes using query search only

        Returns True

        """

        self.prep_download()

        search = Entrez.esearch(db=self._db, term=self._search_string,
            retmax=1, usehistory="y")

        search_results = Entrez.read(search)
        search.close()

        count = int(search_results["Count"])
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        self.logger.info("%i genes sequences in NCBI matching query %s" % (count, self._search_string))

        self.fetch(count, webenv, query_key)
       
        return True


    def download_by_accession(self, input_file, fasta_format=False):
        """Perform download of genes by ID.

        Returns True

        """

        def batch(sequence, size):
            l = len(accessions)
            for start in range(0, l, size):
                yield sequence[start:min(start + size, l)]

        def process_batch(accessions_batch):
            # get GI for query accessions
            query = " OR ".join(("\"%s\"[Accn]" %x for x in accessions_batch))
            search = Entrez.esearch(db=self._db, term=query,
                retmax=1, usehistory="y")

            search_results = Entrez.read(search)
            search.close()

            count = int(search_results["Count"])
            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]

            self.logger.info("%i genes sequences in NCBI matching query %s" % (count, query))

            self.fetch(count, webenv, query_key)

            return count
           

        self.prep_download()

        with open(input_file, 'r') as f:

            accessions = []
            if fasta_format:
                for fasta in SeqIO.parse(f, "fasta"):
                    accessions.append(fasta.id)
            else:
                accessions = f.read().splitlines()

            batchsize = 100
            b = 1
            accession_batches = batch(accessions, batchsize)
            for acc_batch in accession_batches:
                count = process_batch(acc_batch)
                n = len(acc_batch)
                self.logger.info("Batch %i: %i submitted, %i retrieved." % (b, n, count))
                b += 1
           
       
        return True



    def fetch(self, count, webenv, query_key):
        """Fetch genbank records from ncbi history server

        Results are written to <%genename%>.gb in the output directory
        defined in the constructor

        """

        # Create directory to store genbank files
        gbout = open(self._tmpfile, 'a')

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


    def prep_download(self):
        """Delete existing genbank files

        There are multiple methods to download genes. Each download method appends
        to the genbank file, so the first time a download is called this file must
        be emptied

            Args:
                None

        """

        if self._first_download and os.path.exists(self._tmpfile):
            os.remove(self._tmpfile)

        self._first_download = False


    def remove_duplicates(self):
        """Remove duplicate genbank records

        There are multiple methods to download genes. Each download method appends
        to the genbank file, so method this is needed to remove duplicates.

            Args:
                None

        """

        input_seq_iterator = SeqIO.parse(open(self._tmpfile, "rU"), "genbank")
        unique_seq_iterator = self.unique(input_seq_iterator)

        output_handle = open(self._gbfile, "w")
        SeqIO.write(unique_seq_iterator, output_handle, "genbank")
        output_handle.close()

        os.remove(self._tmpfile)

       
    def unique(self, iterable, key=lambda x: x.id):
        """Generator returns unique genbank records

        """

        seen = set()
        for elem, ekey in ((e, key(e)) for e in iterable):
            if ekey not in seen:
                yield elem
                seen.add(ekey)


    def parse(self):
        """Extract sequence and subtype for valid genes in each genbank record

        """

        # Remove possible duplicates in genbank file
        self.remove_duplicates()

        subtype_counts = Counter()

        # Output files
        with open(self._fastafile, 'w') as ffh, open(self._subtypefile, 'w') as stfh: 

            for gb_record in SeqIO.parse(open(self._gbfile,"r"), "genbank"):

                untyped_features = []
                allele = 0
                record_subtypes = set([])

                # Find the CDS matching the gene
                # Collect general info about organism in source feature
                for (index, feature) in enumerate(gb_record.features):

                    if feature.type == 'source':
                        # Check for subtype in source feature
                        # There should only be one source feature
                        record_subtypes = self.subtype_parser.search_source(feature)

                    if feature.type == 'CDS':

                        # Check if this CDS is what we are looking for
                        field = None
                        matched = False
                        if 'gene' in feature.qualifiers:
                            field = feature.qualifiers['gene'][0].lower()
                            if field and any(x.search(field) for x in self._genes): 
                                matched = True
                        if 'product' in feature.qualifiers:
                            field = feature.qualifiers['product'][0].lower()
                            if field and any(x.search(field) for x in self._genes): 
                                matched = True

                        #self.logger.debug("Found gene? {} {} {}".format(gb_record.id, matched, field))

                        # Check for tags that indicate this annotation was obtained through non-experimental means
                        if 'inference' in feature.qualifiers:
                            matched = False

                        # Found
                        if matched:

                            # Unique naming
                            allele += 1
                            name = "{}_allele{}".format(gb_record.id, allele)
                            seq = None

                            # Get sequence
                            if self._dna_sequences:

                                try:
                                    seq = feature.extract(gb_record.seq)
                                except Exception:
                                    self.logger.info('Error in record: %s, missing sequence in feature %s' % (name, str(feature)))
                                    continue
                            
                            else:

                                if feature.qualifiers['translation']:
                                    seq = feature.qualifiers['translation'][0]
                                else:
                                    dnaseq = None
                                    try:
                                        dnaseq = feature.extract(gb_record.seq)
                                    except Exception:
                                        self.logger.info('Error in record: %s, missing sequence in feature %s' % (name, str(feature)))
                                        continue
                            
                                    seq = dnaseq.translate(table=11, to_stop=True)

                            # Check if sequence passes conditions
                            if not self.gene_filter.test_sequence(seq):
                                self.logger.info('Record filtered %s, based on sequence condition' % (name))
                                continue

                            # Check if feature passes conditions
                            if not self.gene_filter.test_genbank(feature):
                                self.logger.info('Record filtered %s, based on genbank feature condition' % (name))
                                continue

                            # Try to find subtype
                            subtypes = self.subtype_parser.search_feature(feature)

                            # There should only be one unique subtype per feature
                            if len(subtypes) > 1:
                                self.logger.info('Error in record: %s, multiple subtypes %s in feature %s' % (name, subtypes, str(feature)))

                            elif len(subtypes) == 1:
                                # Found a valid subtype
                                st = subtypes.pop()
                                ffh.write(">{}\n{}\n".format(name, seq))
                                stfh.write("{}\t{}\n".format(name, st))
                                subtype_counts[st] += 1 

                            else:
                                # Failed to find subtype info associated with feature
                                # Keep feature info in case there is information in the top-level record or source feature
                                # but need to know if there are multiple alleles
                                untyped_features.append((name, seq))

                if len(untyped_features) == 1 and allele < 2:
                    # Only one allele found in record

                    # Search annotation fields for info on subtype
                    subtypes = self.subtype_parser.search_annotations(gb_record)
                    record_subtypes.update(subtypes)

                    # Search record fields for info on subtype
                    subtypes = self.subtype_parser.search_record(gb_record)
                    record_subtypes.update(subtypes)

                    # There should only be one unique subtype per record (since only one allele)
                    if len(record_subtypes) > 1:
                        self.logger.info('Error in record: %s, multiple subtypes %s in record %s' % (name, subtypes, str(gb_record)))

                    elif len(record_subtypes) == 1:
                        # Found a valid subtype
                        tup = untyped_features[0]
                        st = record_subtypes.pop()
                        ffh.write(">{}\n{}\n".format(tup[0], tup[1]))
                        stfh.write("{}\t{}\n".format(tup[0], st))
                        subtype_counts[st] += 1

                    else:
                        self.logger.debug('Note: no subtype found in record %s' % (str(gb_record)))
                        print(gb_record)

                else:
                    self.logger.debug('Note: %i alleles found in record %s' % (len(untyped_features), str(gb_record)))


        self.logger.debug('Subtypes encountered:\n{}\n'.format(str(subtype_counts)))

        return None


    def parse_subtype(self):
        """Extract subtype in each genbank record

        """

        # Remove possible duplicates in genbank file
        self.remove_duplicates()

        subtype_counts = Counter()

        # Output files
        with open(self._subtypefile, 'w') as stfh: 

            for gb_record in SeqIO.parse(open(self._gbfile,"r"), "genbank"):

                record_subtypes = set([])
                name = gb_record.id

                # Collect general info about organism in source feature and gene's CDS feature
                for (index, feature) in enumerate(gb_record.features):

                    if feature.type == 'source':
                        # Check for subtype in source feature
                        # There should only be one source feature
                        subtypes = self.subtype_parser.search_source(feature)
                        record_subtypes.update(subtypes)

                    if feature.type == 'CDS':

                        # Check if this CDS is what we are looking for
                        field = None
                        matched = False
                        if 'gene' in feature.qualifiers:
                            field = feature.qualifiers['gene'][0].lower()
                            if field and any(x.search(field) for x in self._genes): 
                                matched = True
                        if 'product' in feature.qualifiers:
                            field = feature.qualifiers['product'][0].lower()
                            if field and any(x.search(field) for x in self._genes): 
                                matched = True

                        #self.logger.debug("Found gene? {} {} {}".format(gb_record.id, matched, field))

                        # Check for tags that indicate this annotation was obtained through non-experimental means
                        if 'inference' in feature.qualifiers:
                            matched = False

                        # Found
                        if matched:

                            # Try to find subtype
                            subtypes = self.subtype_parser.search_feature(feature)
                            record_subtypes.update(subtypes)      
                            
                # Search annotation fields for info on subtype
                subtypes = self.subtype_parser.search_annotations(gb_record)
                record_subtypes.update(subtypes)

                # Search record fields for info on subtype
                subtypes = self.subtype_parser.search_record(gb_record)
                record_subtypes.update(subtypes)

                # There should only be one unique subtype per record (since only one allele)
                if len(record_subtypes) > 1:
                    self.logger.info('Error in record: %s, multiple subtypes %s in record %s' % (name, subtypes, str(gb_record)))

                elif len(record_subtypes) == 1:
                    # Found a valid subtype
                    st = record_subtypes.pop()
                    stfh.write("{}\t{}\n".format(name, st))
                    subtype_counts[st] += 1

                else:
                    self.logger.debug('Note: no subtype found in record %s' % (str(gb_record)))

                
        self.logger.debug('Subtypes encountered:\n{}\n'.format(str(subtype_counts)))

        return None
