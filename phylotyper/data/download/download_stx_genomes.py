#!/usr/bin/env python

"""Download PMIDs for references related to pubmed journals


Example:
        $ python download_stx_genomes.py input output

"""

from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import argparse
import logging
import re
import os
import time

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mdwhitesi@gmail.com"


logging.basicConfig(
    #filename='download_.log',
    level=logging.DEBUG,
    format='%(asctime)s %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
    filemode='w')
logger = logging.getLogger(__name__)

Entrez.email = 'mdwhitesi@gmail.com'


def retrieve(idlist, output_directory):

    logger.debug("trying to retrieve {}".format(idlist))

    genome = None
    attempt = 1
    while attempt <= 3:
        try:
            print()
            fetch = Entrez.efetch(db='nuccore', id=idlist, retmode='text', rettype='fasta')

            header = fetch.readline()
            match = re.search(r'^>(\S\S\S\S+)', header)
            if not match:
                raise Exception('Invalid fasta header: {}'.format(header))

            genome = match.group(1)
            filename = os.path.join(output_directory, genome.replace('|','_') + '.fasta')

            with open(filename, 'w') as out:
                out.write(header)
                for line in fetch:
                    out.write(line)

            fetch.close()

            attempt = 4
            
        except Exception as e:
            logger.warning("Received error from server {}".format(e))
            logger.warning("Attempt %i of 3".format(attempt))
            attempt += 1
            time.sleep(15)

    logger.debug("{} retrieve returned {}".format(idlist, genome))

    return(genome)



if __name__ == "__main__":
    """Run

    """

    # Parse command-line args
    parser = argparse.ArgumentParser(description='Download and store stx genomes')
    parser.add_argument('output_directory', action="store")
    
    args = parser.parse_args()

    stx1fh = open('stx1.tsv', 'a')
    stx2fh = open('stx2.tsv', 'a')

    inputs = [
        # ["label",(),("2d"), "EF584538"],
        # ["label",(),("2g"), "AB048227"],
        ["label",("1a"),(),"Z36899"],
        ["label",(),("2c"), "AM982821"],
        ["label",(),("2e"), "AJ313016"],
        ["label",(),("2d"), "AF329817"],
        ["label",(),("2a"), "EF441620"],
        ["label",(),("2b"), "AB048238"],
        ["label",(),("2c"), "AY633473"],
        ["label",(),("2b"), "AB048222"],
        ["label",(),("2d"), "FN182287"],
        ["label",(),("2a"), "FM998856"],
        ["label",(),("2c"), "L11079"],
        ["label",(),("2e"), "AM939642"],
        ["label",(),("2b"), "AJ567997"],
        ["label",(),("2d"), "AY443047"],
        ["label",(),("2b"), "X65949"],
        ["label",(),("2a"), "FM998861"],
        ["label",(),("2d"), "AF500190"],
        ["label",(),("2e"), "X81415"],
        ["label",(),("2e"), "X81417"],
        ["label",(),("2d"), "EU816447"],
        ["label",("1c"),(),"AB071621"],
        ["label",(),("2a"), "AY443052"],
        ["label",(),("2e"), "FM998838"],
        ["label",(),("2b"), "AJ313015"],
        ["label",(),("2a"), "AY633459"],
        ["label",(),("2a"), "AY443054"],
        ["label",(),("2a"), "AF461173"],
        ["label",(),("2a"), "AF461171"],
        ["label",(),("2b"), "AB048225"],
        ["label",(),("2d"), "FM998848"],
        ["label",(),("2c"), "AF291819"],
        ["label",(),("2f"), "AB472687"],
        ["label",("1a"),(),"Z36900"],
        ["label",(),("2a"), "GQ429168"],
        ["label",(),("2a"), "FM998839"],
        ["label",(),("2a"), "Z37725"],
        ["label",(),("2a"), "AJ272135"],
        ["label",("1a"),(),"ECOSLTTI"],
        ["label",(),("2e"), "FM998846"],
        ["label",(),("2c"), "AY633469"],
        ["label",(),("2d"), "AF479829"],
        ["label",("1c"),(),"DQ449666"],
        ["label",(),("2a"), "EF441618"],
        ["label",(),("2a"), "AF524944"],
        ["label",(),("2d"), "EU816439"],
        ["label",("1c"),(),"AB048237"],
        ["label",(),("2g"), "AY286000"],
        ["label",("1d"),(),"AY170851"],
        ["label",(),("2e"), "FN182286"],
        ["label",(),("2b"), "AB048226"],
        ["label",(),("2c"), "AY739670"],
        ["label",(),("2e"), "AY332411"],
        ["label",("1d"),(),"AB050958"],
        ["label",(),("2c"), "AY633467"],
        ["label",(),("2c"), "FM998860"],
        ["label",(),("2b"), "AB012102"],
        ["label",(),("2a"), "FM998852"],
        ["label",(),("2d"), "X67514"],
        ["label",(),("2c"), "AY443044"],
        ["label",(),("2d"), "AY633457"],
        ["label",(),("2c"), "EU086525"],
        ["label",(),("2a"), "FM998854"],
        ["label",(),("2g"), "AB048236"],
        ["label",(),("2b"), "AB048223"],
        ["label",(),("2b"), "L11078"],
        ["label",("1c"),(),"AB071620"],
        ["label",(),("2d"), "FM998840"],
        ["label",("1a"),(),"AM230663"],
        ["label",(),("2a"), "EF441609"],
        ["label",("1c"),(),"AB071622"],
        ["label",(),("2a"), "GQ429170"],
        ["label",(),("2c"), "GQ429167"],
        ["label",("1c"),(),"AB071624"],
        ["label",(),("2a"), "FM998851"],
        ["label",(),("2c"), "AF461167"],
        ["label",("1d"),(),"AB050959"],
        ["label",(),("2e"), "U72191"],
        ["label",(),("2c"), "AY443049"],
        ["label",(),("2d"), "EF441621"],
        ["label",(),("2b"), "AB048224"],
        ["label",(),("2d"), "X67515"],
        ["label",(),("2d"), "AF500192"],
        ["label",(),("2b"), "AJ567995"],
        ["label",(),("2d"), "GQ429172"],
        ["label",("1a"),(),"AM230662"],
        ["label",(),("2e"), "AJ567998"],
        ["label",(),("2a"), "haemolyticus"],
        ["label",(),("2b"), "AF043627"],
        ["label",(),("2e"), "AM904726"],
        ["label",("1a"),(),"AB083044"],
        ["label",(),("2c"), "M59432"],
        ["label",(),("2a"), "GQ429163"],
        ["label",(),("2c"), "AY633470"],
        ["label",("1c"),(),"AB071619"],
        ["label",(),("2a"), "EF441599"],
        ["label",(),("2a"), "AB030484"],
        ["label",(),("2d"), "AF479828"],
        ["label",(),("2a"), "FM998853"],
        ["label",(),("2f"), "AJ010730"],
        ["label",(),("2g"), "AJ966783"],
        ["label",(),("2b"), "EF441616"],
        ["label",(),("2a"), "FM998842"],
        ["label",(),("2f"), "M29153"],
        ["label",(),("2c"), "AY739671"],
        ["label",(),("2e"), "M21534"],
        ["label",(),("2d"), "AY443048"],
        ["label",(),("2c"), "AY443043"],
        ["label",(),("2c"), "AY633464"],
        ["label",(),("2c"), "AB015057"],
        ["label",(),("2e"), "X81416"],
        ["label",(),("2a"), "EF441613"],
        ["label",(),("2c"), "AB071845"],
        ["label",(),("2a"), "X07865"],
        ["label",(),("2c"), "FM177471"],
        ["label",(),("2d"), "AF500191"],
        ["label",(),("2a"), "Z50754"],
        ["label",(),("2b"), "AB048229"],
        ["label",(),("2d"), "FM998855"],
        ["label",("1c"),(),"Z36901"],
        ["label",(),("2d"), "AY633458"],
        ["label",(),("2d"), "EF441605"],
        ["label",(),("2d"), "SNS"],
        ["label",(),("2a"), "FM998858"],
        ["label",(),("2c"), "AY443045"],
        ["label",(),("2a"), "AY633472"],
        ["label",(),("2b"), "AB048228"],
        ["label",(),("2e"), "AJ249351"],
        ["label",(),("2d"), "AF500189"],
        ["label",(),("2e"), "FM998844"],
        ["label",(),("2d"), "X61283"],
        ["label",("1c"),(),"AB071623"],
        ["label",(),("2a"), "EF441619"],
        ["label",(),("2a"), "GQ429162"],
        ["label",(),("2c"), "EF441604"],
        ["label",(),("2d"), "DQ235775"],
    ]

    for acc in inputs:

        logger.debug('Trying {}'.format(acc[3]))
        result_handle = NCBIWWW.qblast('blastn', 'nt', acc[3],
            megablast=True,
            expect=1e-10,
            perc_ident=100,
            alignments=10,
            descriptions=10,
            hitlist_size=10,
            filter='none',
            entrez_query='{}[Organism]'.format(
            "\""+'Escherichia coli'+"\"")
        )

        blast_records = NCBIXML.parse(result_handle)

        for blast_record in blast_records:
            
            for alignment in blast_record.alignments:

                subj = alignment.title
                match = re.search(r'^(\S\S\S\S+)', subj)
                if not match:
                    raise Exception('Invalid blast title: {}'.format(subj))
                subjid = match.group(1)
                hsp = alignment.hsps[0]
                slen = alignment.length

                # Only want large regions
                if slen > 200000:
                    # And not at the ends of contigs
                    if hsp.sbjct_start > 5000 and hsp.sbjct_end < (slen-5000):
                        # Is entire region aligned 
                        if hsp.query_start == 1 and hsp.query_end == blast_record.query_length:
                            # and 100% identical
                            if hsp.identities == hsp.align_length:
                                logger.debug('GOT ONE -- {}'.format(subj))
                                genome = retrieve(subjid, args.output_directory)

                                if len(acc[1]) > 0:
                                    if len(acc[1]) > 1:
                                        subtype = ','.join(acc[1])
                                    else:
                                         subtype = acc[1]

                                    stx1fh.write("{}\t{}\n".format(genome, subtype))

                                if len(acc[2]) > 0:
                                    if len(acc[2]) > 1:
                                        subtype = ','.join(acc[2])
                                    else:
                                         subtype = acc[2]

                                    stx2fh.write("{}\t{}\n".format(genome, subtype))

                            else:
                                logger.debug('Not identical -- {}'.format(subj))
                        else:
                            logger.debug('Not full coverage -- {}'.format(subj))
                    else:
                        logger.debug('Not in middle -- {}'.format(subj))
                else:
                    logger.debug('Too short -- {}'.format(subj))
    

    stx1fh.close()
    stx2fh.close()          


