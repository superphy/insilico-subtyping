#!/usr/bin/env python

"""Download M protein (emm) gene sequences from Genbank

Using NCBI FTP site, download all protein sequences in archive:
  ftp://ftp.ncbi.nlm.nih.gov/ncbi-asn1/protein_fasta/

A list of GI's in gamma proteobacteria is obtained using eutils.
The FTP gzipped files are expanded, cross-referenced in the GI list
 and matched GIs are inserted into a sqlite DB. All data is stored 
 in a single table: 'protein':
  protein:
    gi [int]
    dbsrc [text]
    accession [text]
    version [int]
    description [text]
    fasta_header [text]
    sequence [text]
    lastseen_id [int]

Example:
        $ python example_google.py

"""

import argparse
from Bio import Entrez, SeqIO
from ftplib import FTP
import hashlib
import logging
import multiprocessing
from os import remove, path
import re
from subprocess32 import check_call, check_output, call
import tarfile
import time
from urllib import urlretrieve
from psql_db import DBInterface, DBResource

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mdwhitesi@gmail.com"



def download(dbargs, skip_gi, skip_ftp, db_dir, blastdbcmd, diamondcmd, update):
    """Run GI protein sequence download

    Args:
        dbname, dbuser, etc[str]: Postgres connection parameters
        skip_gi[bool]: Don't retrieve new GI set, use GI in database
        skip_ftp[bool]: Don't download new blastdb from ftp site, use DB in db_dir
        db_dir[str]: Location to store downloaded gzipped blast files
        blastdbcmd[str]: Filepath to Blast+ blastdbcmd executable
        diamondcmd[str]: Filepath to Diamond executable
        update[bool]: Don't erase existings sequence and only insert new sequences in DB

    """
   
    # Setup logger
    FORMAT = '%(asctime)s - %(processName)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.DEBUG, format=FORMAT)
    logger = logging.getLogger(__name__)

    # Perform EUtils Search
    search = Entrez.esearch(db="nuccore", term="emm[Gene Name] AND Streptococcus[Organism]",
        usehistory="y")

    # Retrieve matching IDs in batches
    result = Entrez.read(search)
    search.close()

    # Determine number of GIs
    count = int(result["Count"])

    query = {
        'gi_count': count,
        'query_term': "#%s" % result["QueryKey"],
        'webenv': result["WebEnv"],
        'complete': False
    }

    logger.info("%i GI protein sequences in NCBI" % count)
        
    # Return list of query parameters to retrieve in batches of 100\
    batch_size = 100
    starts = range(0,count, batch_size)
    query['start_positions'] = starts
       
    logger.info("%i Downloads will be performed in %i batches" % length(starts))


def retrieve_ids(query_params, batch_size): 
        """Download part of the GI list using
        start positions

        Operation uses Eutils esearch function

        Args:
            query_params[dict]:
                gi_count[int]: Number of GIs
                webenv[str]: Entrez.webenv input
                query_term[str]: Entrez.term input
                start_positions[list]: start positions in batch

        Returns:
            bool: True if successful

        """

        count = query_params['gi_count']

        for start in query_params['start_positions']:
            end = min(count, start+batch_size)
           
            attempt = 1
            data = None
            while attempt <= 3:
                try:
                    search = Entrez.esearch(db="nucleotide", term=query_params['query_term'],
                        retstart=start, retmax=batch_size, webenv=query_params['webenv'], usehistory='y')

                    data = Entrez.read(search)
                    search.close()

                    attempt = 4
                    
                except Exception, e:
                    logger.warning("Received error from server %s" % str(e))
                    logger.warning("Attempt %i of 3 for record %i to %i" % (attempt, start+1, end))
                    attempt += 1
                    time.sleep(15)

            if not data:
                logger.warning("Eutils search request failed. No results returned for record %i to %i" % (start+1, end))
                return False


            bs = end - start;
            if len(data["IdList"]) != bs:
                logger.warning("Missing IDs in Eutils search request for record %i to %i.  Requested: %i, Got: %i" % (start+1, end, len(data["IdList"]), bs))
                return False
              
            logger.debug("%i to %i of %i retrieved" % (start+1, end, count))

            # Pull out GIs
            dbi.add_gi(map(int, data["IdList"]))

        return True


 

    ncbi = NCBIUtils(logger=logger, db_dir=db_dir, blastdbcmd=blastdbcmd, diamondcmd=diamondcmd)
    if not skip_gi:
        # Download protein GI that match query

        query = ncbi.prepare_query()

        logger.info("Downloading protein GI's matching query")
        logger.info("Query split into %i jobs." % len(query))
        logger.info("Retrieving GI list...")

        download_attempt = 0
        remaining = len(query)
        while(download_attempt < 10):
            # Start
            pool = multiprocessing.Pool(processes=4)
            async_results = []
            j = 0
            for q in query:

                if not q["complete"]:
                    async_results.append(pool.apply_async(gi_iteration, 
                        args=(j, q, dbargs, db_dir), 
                        callback=record_download_result))

                j=j+1

            pool.close()
            # Raise exceptions in processes as soon as they occur
            for r in async_results:
                r.get()
            pool.join()

            this_remaining = sum(q["complete"] == False for q in query)
            if this_remaining:
                # Failed / incomplete jobs found

                if this_remaining < remaining:
                    # Incomplete jobs reduced during last run
                    # Making progress
                    remaining = this_remaining
                    download_attempt += 1
                    time.sleep(60)
                    logger.info("GI download round %i. Remaining Jobs: %i", download_attempt, remaining)

                else:
                    logger.info("GI download round %i made no progress. Aborting.", download_attempt)
                    download_attempt = 10
            else:
                # No remaining jobs, done
                download_attempt = 10

        
        incomplete = [ x for x in query if x["complete"] == False ]
        if incomplete and len(incomplete):
            logger.error("GI download failed.")
            logger.error("The following jobs did not complete:")
            for i in incomplete:
                logger.error(repr(i))

            raise Exception("GI Download failed")
        else:
            logger.info("GI download complete.")

        # Check GI count in db, cleanup obsolete GIs
        dbi = DBInterface(logger=logger, dsn=dbargs)
        with DBResource(dbi):
            dbi.check_gi(ncbi.gi_count)
           


    # Download sequences from ftp site
    if not skip_ftp:
        download_attempt = 0
        query = ncbi.ftp_files()
        remaining = len(query)

        logger.info("Downloading FTP files with NR protein sequences")
        logger.info("Query split into %i jobs." % remaining)
        logger.info("Downloading FTP files...")
        while(download_attempt < 10):
            # Start
            pool2 = multiprocessing.Pool(processes=4)
            async_results = []
            j = 0
        
            for q in query:

                if not q["complete"]:
                    async_results.append(pool2.apply_async(ftp_iteration, 
                        args=(j, q, dbargs, db_dir, blastdbcmd), 
                        callback=record_download_result))

                j=j+1

            pool2.close()
            # Raise exceptions in processes as soon as they occur
            for r in async_results:
                r.get()
            pool2.join()

            this_remaining = sum(q["complete"] == False for q in query)
            if this_remaining:
                # Failed / incomplete jobs found

                if this_remaining < remaining:
                    # Incomplete jobs reduced during last run
                    # Making progress
                    remaining = this_remaining
                    download_attempt += 1
                    time.sleep(60)
                    logger.info("FTP download round %i. Remaining Jobs: %i", download_attempt, remaining)

                else:
                    logger.info("FTP download round %i made no progress. Aborting.", download_attempt)
                    download_attempt = 10
            else:
                # No remaining jobs, done
                download_attempt = 10
            
        
        logger.info("FTP download complete.")

    # Move sequences from blastdb to our db
    dbi = DBInterface(logger=logger, dsn=dbargs)
    with DBResource(dbi):
        ncbi.transfer(dbi);
       
        # Output fasta file
        fasta_file = db_dir + '/nr.ffn'
        dbi.output_fasta(fasta_file)

        # Create diamond DB
        ncbi.build_diamond_db(fasta_file)

        remove(fasta_file)




def record_download_result(job_id):
    """Callback function when parallel job ends.
    Records successful jobs in gi_query and ftp_query 
    by setting flag:
        query[job_id]["complete"] = True
    when job_id > -1.

    Args:
        job_id[int]: successful job ID or -1

    Returns:
        None

    """

    if job_id > -1:
        query[job_id]["complete"] = True


def gi_iteration(job_id, query_params, dbargs, db_dir):
    """Download batch of GI protein sequences
    and insert into database

    This function is called in parallel so required
    objects are created inside function

    Args:
        query_params[dict]:
            gi_count[int]: Number of GIs
            webenv[str]: Entrez.webenv input
            query_term[str]: Entrez.term input
            start_positions[list]: start positions in batch
        dbargs[dict]: Sqlite DB filepath
        db_dir[str]: tmp dir location

    Returns:
        int: job ID or -1 if error occurred

    """

    logger = logging.getLogger(__name__)
   
    process = multiprocessing.current_process()
    fn = "pid{0}.log".format(process.pid)

    if not path.exists(fn):
        fh = logging.FileHandler(fn)
        fh.setLevel(logging.DEBUG)
        fh_format = logging.Formatter('%(asctime)s - %(processName)s - %(levelname)s - %(message)s')
        fh.setFormatter(fh_format)
        logger.addHandler(fh)

    # Initialize objects
    dbi = DBInterface(logger=logger, dsn=dbargs)
    ncbi = NCBIUtils(logger=logger, db_dir=db_dir)

    with DBResource(dbi):
        # Retrieve GI chunk
        gi_result = ncbi.retrieve_gi_batch(query_params, dbi)

        if not gi_result:
            logger.debug("GI Download failed for job %i."%job_id)
            return -1
        
        logger.debug("JOB %i GI dowload succeeded."%job_id)


    return job_id


def ftp_iteration(job_id, query_batch, dbargs, db_dir): 
    """Download batch of FTP fasta files
    and insert sequences into database

    This function is called in parallel so required
    objects are created inside function

    """

    logger = logging.getLogger(__name__)
   
    process = multiprocessing.current_process()
    fn = "pid{0}.log".format(process.pid)

    if not path.exists(fn):
        fh = logging.FileHandler(fn)
        fh.setLevel(logging.DEBUG)
        fh_format = logging.Formatter('%(asctime)s - %(processName)s - %(levelname)s - %(message)s')
        fh.setFormatter(fh_format)
        logger.addHandler(fh)
 
    ncbi = NCBIUtils(logger=logger, db_dir=db_dir)
    ok = ncbi.download_ftp_batch(query_batch)
    if ok:
        logger.debug("JOB %i FTP dowload succeeded." % job_id)
    else:
        logger.debug("FTP Download failed for job %i." % job_id)
        return -1

    return job_id


class NCBIUtils(object):
    """FTP and Eutils access to NCBI

    """

    def __init__(self, db_dir, ftp_host='ftp.ncbi.nlm.nih.gov', ftp_dir='blast/db/', blastdbcmd=None, diamondcmd=None, logger=None):
        """

        Args:
            
   

        """

        # Connection params
        self._ftp_host = ftp_host
        self._ftp_dir = ftp_dir
        self._db_dir = db_dir

        # Blast exe command
        self._bdbcmd = blastdbcmd
        if self._bdbcmd:
            # Check if its exe
            check_output([self._bdbcmd, "-h"])

        # Diamond exe command
        self._dcmd = diamondcmd
        if self._dcmd:
            # Check if its exe
            check_output([self._dcmd, "-h"])

        # Blast exe command
        self._bdbcmd = blastdbcmd
        if self._bdbcmd:
            # Check if its exe
            check_call([self._bdbcmd, "-h"])

        # NCBI config
        self._fasta_pattern = re.compile('^gi\|(\d+)\|(\w+)\|([\w\.]+)\|(?:[\w\.\|]+)*(\s+.+)?$')
        self._gi_batch = 5000
        self._n_gi_batch = 10
        self._gi_batch_size = self._n_gi_batch * self._gi_batch
        self._gi_count = None
        self._search_term = "txid1236[Organism:exp]"

        # Set email for eutils interface
        Entrez.email = "superphy.info@gmail.com"
       
        self.logger = logger or logging.getLogger(__name__)


    @property
    def gi_count(self):
        return self._gi_count

    @property
    def gi_batch_size(self):
        return self._gi_batch_size
    

    def prepare_query(self):
        """Initiate webenv key for query and retrieve total number of GIs
        matching query.

        Break the GI list into batches of start positions.

        Args:
            None

        Returns:
            Generator producing dicts with query parameters

        """

        # Run esearch query, set webenv history
        search = Entrez.esearch(db="protein", term=self._search_term,
            retmax=1, usehistory="y", maxdate=2014)

        result = Entrez.read(search)
        search.close()

        # Determine number of GIs
        count = int(result["Count"])

        ## TEST, REMOVE LATER
        count = 100000
        ##
        
        query = {
            'gi_count': count,
            'query_term': "#%s" % result["QueryKey"],
            'webenv': result["WebEnv"],
            'complete': False
        }

        self.logger.info("%i GI protein sequences in NCBI for gamma-proteobacteria" % count)
        
        self._gi_count = count

        # Return list of query parameters with batch
        # of start positions
        return self.gi_block(query)
       
    
    def gi_block(self, query_params):
        """List of input parameters
        combined with a range of start positions
        for function: retrieve_gi_batch()

        Args:
            dict:
                gi_count[int]: Number of GIs
                webenv[str]: Entrez.webenv input
                query_term[str]: Entrez.term input

        Returns:
            List of dict with the following keys
                gi_count[int]: Number of GIs
                webenv[str]: Entrez.webenv input
                query_term[str]: Entrez.term input
                start_positions[list]: start positions in batch

        """

        param_list = []
        count = query_params['gi_count']
        starts = range(0,count,self._gi_batch)

        n = self._n_gi_batch
        for i in xrange(0, len(starts), n):
            st = starts[i:i+n]
            paramset = query_params.copy()
            paramset['start_positions'] = st
            param_list.append(paramset)

        return param_list


    


    def ftp_files(self):
        """Connect to NCBI protein FTP. Retrieve protein sequence
        file list.

        Break the file list into batches.

        Args:
            None

        Returns:
            List with ftp download parameters

        """

        # Connect to ftp site
        ftp = FTP(self._ftp_host)
        ftp.login()
        ftp.cwd(self._ftp_dir)

        # Get file list
        files = ftp.nlst()
        
        if len(files) == 0:
            raise Exception("FTP directory %s contains no files" % self._ftp_dir)

        p = re.compile('^nr\..+tar.gz$')
        nrfiles = [ x for x in files if p.match(x) ]

        if len(nrfiles) == 0:
            raise Exception("FTP directory %s contains no files with pattern 'nr.XX.tar.gz" % self._ftp_dir)
        
        self.logger.info("%i gzipped nr blast db files found in %s." % (len(nrfiles), self._ftp_dir))

        ftp.quit()

        return self.ftp_file_block(nrfiles)


    def ftp_file_block(self, file_list):
        """Function produces List of ftp jobs i.e.
        blocks of ftp files for download.

        Args:
            file_list[list]: list of files in ftp directory

        Returns:
            List of dicts with key/values:
                ftp_batch[list]: contains lists.
                    Each element contains list with FTP address [0] and file basename [1]
                complete[bool] status

        """
        flist = []
        n = 5
        for i in xrange(0, len(file_list), n):
            f = [ ["ftp://{0}/{1}/{2}".format(self._ftp_host, self._ftp_dir, x), x] for x in file_list[i:i+n] ]
            query = {
                'ftp_batch': f,
                'complete': False
            }

            flist.append(query)

        return flist


    def download_ftp_batch(self, batch):
        """Download, extract ftp files

        Args:
            batch[dict]: ftp block where "ftp_batch" key
                contains list of lists with ftp address and file basename

        Returns:
            bool: True if successful

        """

        files = batch["ftp_batch"]
        self.logger.debug("Files in this batch %s " % ', '.join([x[1] for x in files]))

        i = 0
        n = len(files)
        for ftpset in files:
            url, f = ftpset

            # Download file
            dbfile = self._db_dir + '/' + f
            md5file = dbfile + '.md5'
            md5url = url + '.md5'

            try:
                urlretrieve(url, dbfile)
                urlretrieve(md5url, md5file)
            except:
                self.logger.warning("Download of %s failed" % url)
                return False
            self.logger.debug("Downloaded db file %s" % url)

            # Check checksum
            try:
                expected_sum = self.read_md5_file(md5file)
                observed_sum = self.md5sum(dbfile)

                if expected_sum != observed_sum:
                    self.logger.warning("Checksum values for download %s does not match.", url)
                    return False

            except Exception as e:
                self.logger.warning("Checksum check for %s failed: " % (url, str(e)))
                return False

            remove(md5file)

            # Extract file
            try:
                self.extract(dbfile)

            except Exception as e:
                self.logger.warning("Tar extraction for %s failed: " % (url, str(e)))
                return False

            i+=1
            remove(dbfile)
            self.logger.debug("%i of %i downloaded and extracted" % (i, n))
            

        return True

    def read_md5_file(self, filename):
        """Extract md5 checksum value from text file

        Args:
            filename[str]: file path

        Returns:
            str: MD5 checksum

        """

        with open(filename, "r") as f:
            chksum_line = f.read()
            chksum = chksum_line.split()

            return chksum[0]


    def md5sum(self, filename, blocksize=65536):
        """Compute MD5 checksum for binary file

        Args:
            filename[str]: file path
            blocksize[int]: bite size of read block

        Returns:
            str: MD5 checksum

        """

        hash = hashlib.md5()
        with open(filename, "rb") as f:
            for block in iter(lambda: f.read(blocksize), b""):
                hash.update(block)
        return hash.hexdigest()


    def extract(self, tgzfile):
        """
        Extract all files in tar.gz archive that start with nr*

        Args:
            tgzfile[str]: tar.gz file path
            location[str]: Destination directory
           
        Returns:
            bool: True on success

        """

        tar = tarfile.open(tgzfile)
        tar.extractall(path=self._db_dir, members=NCBIUtils._ok_files(tar))
        tar.close()


    @staticmethod
    def _ok_files(members):
        """Makes sure the files in the tar archive
        are not being extracted outside the intented
        directory

        """
        p = re.compile('^nr')
        for tarinfo in members:
            if p.match(tarinfo.name):
                yield tarinfo


    def transfer(self, dbi):
        """Transfer sequences from blastdb to psql protein DB
        using blastdbcmd.

        Args:


        """

        page = 0
        page_size = 10000
        while True:
            gi_batch = dbi.fetch_gi_page(page, page_size)
            self.logger.debug("Got of batch of GIs: %i" % len(gi_batch))

            if len(gi_batch) == 0:
                break

            self.get_sequences(gi_batch, dbi)
            self.logger.debug("Iterated through GIs from %i to %i" % (page,page+page_size))

            page+=page_size

        dbi.insert_sequences()
        


    def get_sequences(self, gi_list, dbi):
        """Extract sequences from blastdb
        using blastdbcmd.

        Args:


        """

        # Batch file
        batch_file = self._db_dir + '/batch.gi'
        with open(batch_file, "w") as w:
            for gi in gi_list:
                w.write('%i\n' % gi)

        # Run blastdbcmd
        fasta_file = self._db_dir + '/batch.ffp'
        db = self._db_dir + '/nr'
        call([self._bdbcmd, "-db", db, "-entry_batch", batch_file, "-out", fasta_file])

        # Parse fasta output
        gi_found = {}
        for seq in SeqIO.parse(fasta_file, 'fasta'):
            header = seq.description
            # split and use first entry in header
            ids = header.split('>')
            m = re.search(self._fasta_pattern, ids[0])
            if m:
                gi = m.group(1)
                dbsrc = m.group(2)
                acc = m.group(3)
                desc = m.group(4)
                proteinseq = str(seq.seq)
                dbi.add_sequences(gi, dbsrc, acc, desc.lstrip(), proteinseq)
                gi_found[gi] = gi_found.get(gi,0) + 1
           
            else:
                self.logger.warning("Invalid Fasta header %s. Skipping sequence" % seq.id)

        nseq = len(gi_found)
        ngi = sum(gi_found.values())
        self.logger.debug("Given %i GIs" % len(gi_list))
        self.logger.debug("Found and loaded %i sequences into DB" % nseq)
        self.logger.debug("%i sequences represent %i GIs (due to duplicate sequences)" % (nseq, ngi))

        
        if nseq == 0:
            self.logger.warning("No sequences loaded. No matching GIs were found in blast database for this batch.")

        remove(fasta_file)
        remove(batch_file)



    def build_diamond_db(self, fasta_file):
        """Build Diamond database from fasta input

        Args:
            fasta_file[str]: Fasta filepath containing protein sequences

        """
        db = self._db_dir + '/nr'
        check_call([self._dcmd, "makedb", "--db", db, "--in", fasta_file])




if __name__ == "__main__":
    """Run pipeline downloading sequences and loading into DB

    """

    # Parse command-line args
    parser = argparse.ArgumentParser(description='Download and store NCBI gamma-proteobacteria protein sequences in Postgresql DB')
    parser.add_argument('--dbname', action="store")
    parser.add_argument('--dbuser', action="store")
    parser.add_argument('--dbhost', action="store")
    parser.add_argument('--dbport', action="store")
    parser.add_argument('--dbpass', action="store")
    parser.add_argument('--blastdbcmd', action="store", default='blastdbcmd')
    parser.add_argument('--diamondcmd', action="store", default='diamond')
    parser.add_argument('--db_dir', action="store")
    parser.add_argument('--skip_gi', action="store_true")
    parser.add_argument('--skip_ftp', action="store_true")
    parser.add_argument('--update', action="store_true")
    
    args = parser.parse_args()

    # Strip out DB parameters
    dbargs = DBInterface.dbparams(dbname=args.dbname, dbuser=args.dbuser, dbhost=args.dbhost, dbport=args.dbport, dbpass=args.dbpass)

    # Run
    run(dbargs=dbargs, skip_gi=args.skip_gi, skip_ftp=args.skip_ftp, db_dir=args.db_dir, 
        blastdbcmd=args.blastdbcmd, diamondcmd=args.diamondcmd, update=args.update)
    

