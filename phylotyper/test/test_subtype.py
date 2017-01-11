import csv
import os
import re
import shutil
import tempfile
import unittest

from phylotyper.config import PhylotyperOptions
from phylotyper.main import build_pipeline, subtype_pipeline
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

        test_genomes = ['test_stx2_genome1.fasta', 'test_stx2_genome2.fasta']
        subtype_options['genomes'] = [os.path.join(self.data_dir, f) for f in test_genomes]
        subtype_options['output_directory'] = self.test_dir
        subtype_options['ngenomes'] = 2
        subtype_options['fast'] = False
        subtype_options['noplots'] = True

        self.configObj = config
        self.subtype_options = subtype_options
        self.scheme = scheme
    

    def testStx2Unknowns(self):

        self.init('stx2')

        subtype_pipeline(self.subtype_options, self.configObj)

        # Check predictions
        with open(os.path.join(self.test_dir, 'subtype_predictions.csv'), 'rb') as csvfile:
            csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
            row1 = csvreader.next()
            row2 = csvreader.next()

        found_by_identity = True if re.match(r'identical', row1[2]) else False
        above_value = True if float(row2[2]) > .95 else False
     
        self.assertTrue(all([row1[1] == 'a', row2[1] == 'a', above_value, found_by_identity]))


   

# class BuildTests(unittest.TestCase):

#     def setUp(self):
#         # Create temporary directory
#         sandbox = os.path.abspath(os.path.join(os.path.dirname(__file__),'sandbox'))
#         self.test_dir = os.path.abspath(tempfile.mkdtemp(dir=sandbox))
        
#         # Create temporary yaml file for testing
#         self.yamlfile = os.path.join(self.test_dir,'test_index.yaml')
#         output = '# Phylotyper filepaths and options for pre-built subtyping schemes\n\nroot_dir: {}\n\nsubtypes: {{}}\n'.format(
#             self.test_dir)
#         with open(self.yamlfile, 'w') as w:
#             w.write(output)

#         # Inputs
#         self.data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'data'))


#     def tearDown(self):
#         # Remove previous directories created
#         #shutil.rmtree(self.test_dir)
#         pass


#     def setUpBuild(self, aa=False):
#         # Initialize phylotyper arguments for new builds & subtype runs

#         # Test inputs
#         suffix = 'a' if aa else 'n'
#         ref = os.path.join(self.data_dir, 'stx2.ff'+suffix)
#         subt = os.path.join(self.data_dir, 'stx2_subtypes.csv')
#         inpu = os.path.join(self.data_dir, 'test_stx2.ff'+suffix)
#         scheme = 'test_gene'

#         if not os.environ.get('PHYLOTYPER_CONFIG'):
#             msg = 'Missing config file enviroment variable.\nMust provide Phylotyper config file using' \
#                 ' enviroment variable PHYLOTYPER_CONFIG for testing'
#             raise Exception(msg)

#         # Parse .ini config file       
#         config = PhylotyperOptions(os.environ.get('PHYLOTYPER_CONFIG'))

#         # Load subtype options
#         stConfig = SubtypeConfig(self.yamlfile)
    
#         # Check input file exists
#         if not os.path.isfile(ref):
#             msg = 'Invalid/missing reference file argument.'
#             raise Exception(msg)

#         # Check input file exists
#         if not os.path.isfile(inpu):
#             msg = 'Invalid/missing input file argument.'
#             raise Exception(msg)

#         # Check subtype file exists
#         if not os.path.isfile(subt):
#             msg = 'Invalid/missing subtype file argument.'
#             raise Exception(msg)

#         # Create subtype directory & file names
#         build_options = stConfig.create_subtype(scheme, aa)

#         # Save additional build options
#         build_options['input'] = os.path.abspath(ref)
#         build_options['subtype_orig'] = os.path.abspath(subt)
#         build_options['output_directory'] = self.test_dir
#         build_options['fast'] = False

#         self.configObj = config
#         self.build_options = build_options
#         self.subtypeOptionsObj = stConfig
#         self.scheme = scheme
#         self.input = inpu


#     def testSubtypeDNA(self):

#         self.setUpBuild(aa=False)

#         # Set up subtype files
#         build_pipeline(self.build_options, self.configObj)

#         # Save setup
#         self.subtypeOptionsObj.save()

#         # Initialize subtype options
#         subtype_options = self.subtypeOptionsObj.get_subtype_config(self.scheme)

#         # Add user-defined command-line arguments
#         subtype_options['input'] = self.input
#         subtype_options['output_directory'] = self.test_dir
#         subtype_options['fast'] = False
#         subtype_options['noplots'] = True

#         subtype_pipeline(subtype_options, self.configObj)

#         # Check predictions
#         with open(os.path.join(self.test_dir, 'subtype_predictions.csv'), 'rb') as csvfile:
#             csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
#             row1 = csvreader.next()
#             row2 = csvreader.next()

#         found_by_identity = True if re.match(r'identical', row1[2]) else False
#         above_value = True if float(row2[2]) > .95 else False
     
#         self.assertTrue(all([row1[1] == 'a', row2[1] == 'a', above_value, found_by_identity]))


#     def testSubtypeAminoAcid(self):

#         self.setUpBuild(aa=True)

#         # Set up subtype files
#         build_pipeline(self.build_options, self.configObj)

#         # Save setup
#         self.subtypeOptionsObj.save()

#         # Initialize subtype options
#         subtype_options = self.subtypeOptionsObj.get_subtype_config(self.scheme)

#         # Add user-defined command-line arguments
#         subtype_options['input'] = self.input
#         subtype_options['output_directory'] = self.test_dir
#         subtype_options['fast'] = False
#         subtype_options['noplots'] = True

#         subtype_pipeline(subtype_options, self.configObj)

#         # Check predictions
#         with open(os.path.join(self.test_dir, 'subtype_predictions.csv'), 'rb') as csvfile:
#             csvreader = csv.reader(filter(lambda row: row[0]!='#', csvfile), delimiter='\t')
#             row1 = csvreader.next()
#             row2 = csvreader.next()

#         found_by_identity = True if re.match(r'identical', row1[2]) else False
#         above_value = True if float(row2[2]) > .95 else False
     
#         self.assertTrue(all([row1[1] == 'a', row2[1] == 'a', above_value, found_by_identity]))


def main():
    unittest.main()

if __name__ == '__main__':

    main()