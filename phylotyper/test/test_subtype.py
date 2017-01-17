import csv
import os
import re
import shutil
import tempfile
import unittest

from phylotyper.config import PhylotyperOptions
from phylotyper.run import subtype_pipeline
from phylotyper.subtypes_index import SubtypeConfig

class SubtypeTests(unittest.TestCase):

    def setUp(self):
        # Create temporary directory
        sandbox = os.path.abspath(os.path.join(os.path.dirname(__file__),'sandbox'))
        self.test_dir = os.path.abspath(tempfile.mkdtemp(dir=sandbox))

        # Package yaml config file
        self.yamlfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'subtypes_index.yaml')

        # Inputs
        self.data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'data'))


    def tearDown(self):
        # Remove previous directories created
        #shutil.rmtree(self.test_dir)
        pass


    def init(self, scheme):

        if not os.environ.get('PHYLOTYPER_CONFIG'):
            msg = 'Missing config file enviroment variable.\nMust provide Phylotyper config file using' \
                ' enviroment variable PHYLOTYPER_CONFIG for testing'
            raise Exception(msg)

        # Parse .ini config file       
        config = PhylotyperOptions(os.environ.get('PHYLOTYPER_CONFIG'))

        # Load subtype options
        stConfig = SubtypeConfig(self.yamlfile)
        subtype_options = stConfig.get_subtype_config(scheme)

        self.configObj = config
        self.subtype_options = subtype_options
        self.scheme = scheme
        self.subtype_options['fast'] = False
        self.subtype_options['noplots'] = True
        self.subtype_options['output_directory'] = self.test_dir


    def testcase1(self):
        test_genomes = ['test_stx2_genome1.fasta', 'test_stx2_genome2.fasta']
        self.subtype_options['genomes'] = [os.path.join(self.data_dir, f) for f in test_genomes]
        self.subtype_options['ngenomes'] = 2

        
    def testcase2(self):
        test_genomes = ['minigenome1.fasta', 'minigenome2.fasta', 'minigenome3.fasta', 'minigenome4.fasta']
        self.subtype_options['genomes'] = [os.path.join(self.data_dir, f) for f in test_genomes]
        self.subtype_options['ngenomes'] = 4


    def testStx2Unknowns(self):

        self.init('stx2')
        self.testcase1()

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        with open(os.path.join(self.test_dir, 'subtype_predictions.csv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            csvreader.next() # Header
            row1 = csvreader.next()
            row2 = csvreader.next()

        in_range1 = True if float(row1[3]) > .8 and float(row1[3]) < .9 else False
        above_value2 = True if float(row2[3]) > .95 else False
     
        self.assertTrue(all([row1[4] == 'non-significant/undetermined', row2[4] == 'c', in_range1, above_value2]))


    def testStx2Knowns(self):

        self.init('stx2')
        self.testcase2()

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        assignments = {}
        with open(os.path.join(self.test_dir, 'subtype_predictions.csv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            csvreader.next() # Header
            for row in csvreader:
                assignments[row[1]] = row[4]
            
        print assignments

        self.assertTrue(all([
            len(assignments) == 4,
            assignments['minigenome1-allele1'] == 'a',
            assignments['minigenome1-allele2'] == 'a',
            assignments['minigenome2'] == 'g',
            assignments['minigenome3'] == 'g'
        ]))


   
def main():
    unittest.main()

if __name__ == '__main__':

    main()