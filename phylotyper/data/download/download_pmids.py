#!/usr/bin/env python

"""Download PMIDs for references related to pubmed journals


Example:
        $ python download_pmids.py outputfile

"""

import argparse
import logging
from Bio import Entrez
from collections import Counter

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mdwhitesi@gmail.com"


logging.basicConfig(
    filename='download_pmids.log',
    level=logging.DEBUG,
    format='%(asctime)s %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
    filemode='w')
logger = logging.getLogger(__name__)

def relevent_pubs(pmids, outputfile):
    """Find linked pubmed IDs using Entrez's ELink

    Results are written to <%outputfile%>


        Args:
            pmids (list): List of ints
            outputfile (str): Filepath to write results to


    """

    # set email
    Entrez.email = 'superphy.info@gmail.com'

    if not isinstance(pmids, list):
        pmids = [pmids]

    with open(outputfile, 'w') as f:
            
        pmid_string = ','.join(pmids)
        links = Entrez.elink(db='pubmed', dbfrom='pubmed', linkname='pubmed_pubmed', id=pmid_string, cmd='neighbor')
        link_results = Entrez.read(links)
        links.close()

        db_record = link_results[0]
        link_set = db_record['LinkSetDb'][0]
        if not link_set['LinkName'] == 'pubmed_pubmed':
            raise Exception('Elink result LinkSetDB is not expected target: pubmed_pubmed')

        for l in link_set['Link']:
            id = l['Id']
            f.write(id+"\n")


    return None


if __name__ == "__main__":
    """Run pipeline downloading sequences

    """

    # Parse command-line args
    parser = argparse.ArgumentParser(description='Retrieve related PMIDs')
    parser.add_argument('outputfile', action="store")
    parser.add_argument('pmidfile',action="store")
    
    args = parser.parse_args()

    with open(args.pmidfile, 'r') as f:
        pmids = f.read().splitlines()

    logger.debug("Searching linked pubmed articles to PMIDs: {}".format(','.join(pmids)))

    relevent_pubs(pmids, args.outputfile)




   