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


    def genomes1(self):
        test_genomes = ['test_stx2_genome1.fasta', 'test_stx2_genome2.fasta']
        self.subtype_options['genomes'] = [os.path.join(self.data_dir, f) for f in test_genomes]
        self.subtype_options['ngenomes'] = 2

        
    def genomes2(self):
        test_genomes = ['minigenome1.fasta', 'minigenome2.fasta', 'minigenome3.fasta', 'minigenome4.fasta']
        self.subtype_options['genomes'] = [os.path.join(self.data_dir, f) for f in test_genomes]
        self.subtype_options['ngenomes'] = 4


    def thegenome(self):
        test_genomes = ['test_stx2_genome2.fasta']
        self.subtype_options['genomes'] = [os.path.join(self.data_dir, f) for f in test_genomes]
        self.subtype_options['ngenomes'] = 1


    def thatgenome(self):
        test_genomes = ['ecoli_mgy_modified.fasta']
        self.subtype_options['genomes'] = [os.path.join(self.data_dir, f) for f in test_genomes]
        self.subtype_options['ngenomes'] = 1


    def theothergenome(self):
        test_genomes = ['test_eae_genome1.fasta']
        self.subtype_options['genomes'] = [os.path.join(self.data_dir, f) for f in test_genomes]
        self.subtype_options['ngenomes'] = 1


    def theotherothergenome(self):
        test_genomes = ['test_eae_genome2.fasta']
        self.subtype_options['genomes'] = [os.path.join(self.data_dir, f) for f in test_genomes]
        self.subtype_options['ngenomes'] = 1


    def theflicgenome(self):
        test_genomes = ['test_flic_genome1.fasta']
        self.subtype_options['genomes'] = [os.path.join(self.data_dir, f) for f in test_genomes]
        self.subtype_options['ngenomes'] = 1


    def testStx2Unknowns(self):

        self.init('stx2')
        self.genomes1()

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        with open(os.path.join(self.test_dir, 'subtype_predictions.tsv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            csvreader.next() # Header
            row1 = csvreader.next()
            row2 = csvreader.next()

        in_range1 = True if float(row1[3]) > .8 and float(row1[3]) < .9 else False
        above_value2 = True if float(row2[3]) > .95 else False
     
        self.assertTrue(all([row1[4] == 'non-significant/undetermined', row2[4] == 'c', in_range1, above_value2]))


    def testStx2(self):

        self.init('stx2')
        self.thegenome()

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        with open(os.path.join(self.test_dir, 'subtype_predictions.tsv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            csvreader.next() # Header
            row = csvreader.next()

        above_value = True if float(row[3]) > .95 else False
        self.assertTrue(all([row[4] == 'c', above_value]))


    def testStx1(self):

        self.init('stx1')
        self.thegenome()

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        with open(os.path.join(self.test_dir, 'subtype_predictions.tsv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            csvreader.next() # Header
            row = csvreader.next()

        above_value = True if float(row[3]) > .85 else False
        self.assertTrue(all([row[2] == 'a', above_value]))


    def testEaeMissing(self):

        self.init('eae')
        self.thegenome()

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        with open(os.path.join(self.test_dir, 'subtype_predictions.tsv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            csvreader.next() # Header
            row = csvreader.next()

        self.assertTrue(row[4] == 'Subtype loci not found in genome')


    def testEaeIdentical(self):

        self.init('eae')
        self.theothergenome()

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        with open(os.path.join(self.test_dir, 'subtype_predictions.tsv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            csvreader.next() # Header
            row = csvreader.next()

        self.assertTrue(all([bool(re.search(r'^identical', row[3])), row[2] == 'epsilon-1']))


    def testEae(self):

        self.init('eae')
        self.theotherothergenome()

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        with open(os.path.join(self.test_dir, 'subtype_predictions.tsv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            csvreader.next() # Header
            row = csvreader.next()

        above_value = True if float(row[3]) > .95 else False
        self.assertTrue(all([row[2] == 'epsilon-1', above_value]))


    def testFlicIdentical(self):

        self.init('flic')
        self.thegenome()

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        with open(os.path.join(self.test_dir, 'subtype_predictions.tsv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            csvreader.next() # Header
            row = csvreader.next()
            print row

        self.assertTrue(all([bool(re.search(r'^identical', row[3])), row[2] == 'h8']))


    def testFlic(self):

        self.init('flic')
        self.theflicgenome()

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        with open(os.path.join(self.test_dir, 'subtype_predictions.tsv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            csvreader.next() # Header
            row = csvreader.next()

        above_value = True if float(row[3]) > .95 else False
        self.assertTrue(all([row[2] == 'h12', above_value]))


    def testWz(self):

        self.init('wz')
        self.thatgenome()

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        with open(os.path.join(self.test_dir, 'subtype_predictions.tsv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            csvreader.next() # Header

            row = csvreader.next()
            
            above_value = True if float(row[3]) > .95 else False
            self.assertTrue(all([row[2] == 'o16', above_value]))


    def testStx2Knowns(self):

        self.init('stx2')
        self.genomes2()

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        assignments = {}
        with open(os.path.join(self.test_dir, 'subtype_predictions.tsv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            csvreader.next() # Header
            for row in csvreader:
                if row[1] == 'not applicable':
                    ID = row[0]
                    assignments[ID] = row[3]
                else:
                    ID = row[1]
                    assignments[ID] = row[4]
            
        print assignments

        self.assertTrue(all([
            len(assignments) == 4,
            assignments['minigenome1-allele1'] == 'a',
            assignments['minigenome1-allele2'] == 'a',
            assignments['minigenome3'] == 'g',
            bool(re.search(r'^identical', assignments['minigenome2']))
        ]))


   
def main():
    unittest.main()

if __name__ == '__main__':

    main()