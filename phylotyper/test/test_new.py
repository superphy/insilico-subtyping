import os
import shutil
import tempfile
import unittest
from Bio import SeqIO

from phylotyper.config import PhylotyperOptions
from phylotyper.main import build_pipeline
from phylotyper.subtypes_index import SubtypeConfig
from phylotyper.tree.seq import SeqDict
from phylotyper.tree.seqaligner import SeqAligner

class NewTests(unittest.TestCase):

    def setUp(self):
        # Create temporary directory
        sandbox = os.path.abspath(os.path.join(os.path.dirname(__file__),'sandbox'))
        self.test_dir = os.path.abspath(tempfile.mkdtemp(dir=sandbox))
        
        # Create temporary yaml file for testing
        self.yamlfile = os.path.join(self.test_dir,'test_index.yaml')
        output = '# Phylotyper filepaths and options for pre-built subtyping schemes\n\nroot_dir: {}\n\nsubtypes: {{}}\n'.format(
            self.test_dir)
        with open(self.yamlfile, 'w') as w:
            w.write(output)

        # Inputs
        self.data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'data'))


    def tearDown(self):
        # Remove previous directories created
        shutil.rmtree(self.test_dir)


    def setUpPhylotyper(self, aa=False):
        # Initialize phylotyper arguments for new builds

        # Test inputs
        suffix = 'a' if aa else 'n'
        inpu = os.path.join(self.data_dir, 'stx2.ff'+suffix)
        subt = os.path.join(self.data_dir, 'stx2_subtypes.csv')
        scheme = 'test_gene'

        if not os.environ.get('PHYLOTYPER_CONFIG'):
            msg = 'Missing config file enviroment variable.\nMust provide Phylotyper config file using' \
                ' enviroment variable PHYLOTYPER_CONFIG for testing'
            raise Exception(msg)

        
        # Parse .ini config file       
        config = PhylotyperOptions(os.environ.get('PHYLOTYPER_CONFIG'))

        # Load subtype options
        stConfig = SubtypeConfig(self.yamlfile)
    
        # Check input file exists
        if not os.path.isfile(inpu):
            msg = 'Invalid/missing input file argument.'
            raise Exception(msg)

        # Check subtype file exists
        if not os.path.isfile(subt):
            msg = 'Invalid/missing subtype file argument.'
            raise Exception(msg)

        # Create subtype directory & file names
        subtype_options = stConfig.create_subtype(scheme, aa)

        # Save additional build options
        subtype_options['input'] = [os.path.abspath(inpu)]
        subtype_options['subtype_orig'] = os.path.abspath(subt)
        subtype_options['output_directory'] = self.test_dir
        subtype_options['fast'] = False

        self.configObj = config
        self.subtype_options = subtype_options
        self.subtypeOptionsObj = stConfig
        self.scheme = scheme


    def test_SeqDict_oneloci(self):

        self.setUpPhylotyper(aa=False)

        tmpfile = os.path.join(self.test_dir, 'tmp.fasta')
        fasta_file = os.path.join(self.data_dir, 'test_stx2.ffn')
        
        presd = SeqDict()
        presd.build(self.subtype_options['input'], self.subtype_options['subtype_orig'])
        # Output unique set
        presd.write(tmpfile, self.subtype_options['subtype'])
        # Save lookup object
        presd.store(self.subtype_options['lookup'])

        # Load output
        postsd = SeqDict()
        postsd.load(self.subtype_options['lookup'])

        # Test against identical/non-identcal sequence
        fasta_sequences = SeqIO.parse(open(fasta_file, 'r'),'fasta')
        results = []
        expected = [True, False]
        for fasta in fasta_sequences:
            rs = True if postsd.find(str(fasta.seq)) else False
            results.append(rs)

        self.assertEqual(results, expected)

    def testNewAminoAcid(self):

        self.setUpPhylotyper(aa=True)

        # Set up subtype files
        build_pipeline(self.subtype_options, self.configObj)

        # Save setup
        self.subtypeOptionsObj.save()

        # Check output files
        filepaths = self.subtypeOptionsObj.get_subtype_config(self.scheme)
        fasta = SeqIO.index(filepaths['alignment'], 'fasta')

        with open(filepaths['subtype']) as f:
            for i, l in enumerate(f):
                pass

        sd = SeqDict()
        sd.load(filepaths['lookup'])

        n = 91
        self.assertTrue(all([len(fasta) == n, i+1 == n, len(sd.seqs) == n]))

    def testNewDNA(self):

        self.setUpPhylotyper(aa=False)

        # Set up subtype files
        build_pipeline(self.subtype_options, self.configObj)

        # Save setup
        self.subtypeOptionsObj.save()

        # Check output files
        filepaths = self.subtypeOptionsObj.get_subtype_config(self.scheme)
        fasta = SeqIO.index(filepaths['alignment'], 'fasta')

        with open(filepaths['subtype']) as f:
            for i, l in enumerate(f):
                pass

        sd = SeqDict()
        sd.load(filepaths['lookup'])

        n = 120
        self.assertTrue(all([len(fasta) == n, i+1 == n, len(sd.seqs) == n]))


class NewTestsMultiLoci(unittest.TestCase):

    def setUp(self):
        # Create temporary directory
        sandbox = os.path.abspath(os.path.join(os.path.dirname(__file__),'sandbox'))
        self.test_dir = os.path.abspath(tempfile.mkdtemp(dir=sandbox))
        
        # Create temporary yaml file for testing
        self.yamlfile = os.path.join(self.test_dir,'test_index.yaml')
        output = '# Phylotyper filepaths and options for pre-built subtyping schemes\n\nroot_dir: {}\n\nsubtypes: {{}}\n'.format(
            self.test_dir)
        with open(self.yamlfile, 'w') as w:
            w.write(output)

        # Inputs
        self.data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'data'))


    def tearDown(self):
        # Remove previous directories created
        shutil.rmtree(self.test_dir)


    def setUpPhylotyper(self, aa=False):
        # Initialize phylotyper arguments for new builds

        # Test inputs
        suffix = 'a' if aa else 'n'
        inpu = [ os.path.join(self.data_dir, 'stx2a.ff'+suffix),
            os.path.join(self.data_dir, 'stx2b.ff'+suffix)
        ]
        subt = os.path.join(self.data_dir, 'stx2_subtypes.csv')
        scheme = 'test_gene'

        if not os.environ.get('PHYLOTYPER_CONFIG'):
            msg = 'Missing config file enviroment variable.\nMust provide Phylotyper config file using' \
                ' enviroment variable PHYLOTYPER_CONFIG for testing'
            raise Exception(msg)

        # Parse .ini config file       
        config = PhylotyperOptions(os.environ.get('PHYLOTYPER_CONFIG'))

        # Load subtype options
        stConfig = SubtypeConfig(self.yamlfile)
    
        # Check input file exists
        if not os.path.isfile(inpu[1]):
            msg = 'Invalid/missing input file argument.'
            raise Exception(msg)

        # Check subtype file exists
        if not os.path.isfile(subt):
            msg = 'Invalid/missing subtype file argument.'
            raise Exception(msg)

        # Create subtype directory & file names
        subtype_options = stConfig.create_subtype(scheme, aa)

        # Save additional build options
        subtype_options['input'] = [ os.path.abspath(i) for i in inpu ]
        subtype_options['subtype_orig'] = os.path.abspath(subt)
        subtype_options['output_directory'] = self.test_dir
        subtype_options['fast'] = False

        self.configObj = config
        self.subtype_options = subtype_options
        self.subtypeOptionsObj = stConfig
        self.scheme = scheme


    def test_SeqDict_twoloci(self):

        self.setUpPhylotyper(aa=False)

        tmpfile = [ os.path.join(self.test_dir, 'tmp1.fasta'),
            os.path.join(self.test_dir, 'tmp2.fasta')
        ]
        fasta_file = os.path.join(self.data_dir, 'test_stx2.ffn')
        
        presd = SeqDict(2)
        presd.build(self.subtype_options['input'], self.subtype_options['subtype_orig'])
        # Output unique set
        presd.write(tmpfile, self.subtype_options['subtype'])
        # Save lookup object
        presd.store(self.subtype_options['lookup'])

        # Load output
        postsd = SeqDict(2)
        postsd.load(self.subtype_options['lookup'])

        # Test against identical/non-identcal sequence
        fasta_sequences = SeqIO.parse(open(fasta_file, 'r'),'fasta')
        results = []
        expected = [True, False]
        for fasta in fasta_sequences:
            rs = True if postsd.find(str(fasta.seq)) else False
            results.append(rs)

        self.assertEqual(results, expected)


    def test_malign(self):

        self.setUpPhylotyper(aa=False)

        tmpfiles = [ os.path.join(self.test_dir, 'tmp1.fasta'),
            os.path.join(self.test_dir, 'tmp2.fasta')
        ]

        tmpfile = os.path.join(self.test_dir, 'tmp_saln.fasta')

        aln = SeqAligner(self.configObj)
        aln.malign(self.subtype_options['input'], tmpfiles, tmpfile)

        lengths = []
        tmpfiles.append(tmpfile)
        for f in tmpfiles:
            seqs = SeqIO.parse(open(f, 'r'),'fasta')
            lengths.append(len(str(seqs.next().seq)))

        self.assertEqual(lengths[0]+lengths[1], lengths[2])



class SubtypeIndexTests(unittest.TestCase):

    def setUp(self):
        # Create temporary directory
        sandbox = os.path.abspath(os.path.join(os.path.dirname(__file__),'sandbox'))
        self.root_dir = os.path.abspath(tempfile.mkdtemp(dir=sandbox))
        
        # Create temporary yaml file for testing
        self.yamlfile = os.path.join(self.root_dir,'test_index.yaml')
        output = '# Phylotyper filepaths and options for pre-built subtyping schemes\n\nroot_dir: {}\n\nsubtypes: {{}}\n'.format(
            self.root_dir)
        with open(self.yamlfile, 'w') as w:
            w.write(output)


    def tearDown(self):
        # Remove previous directories created
        shutil.rmtree(self.root_dir)
        

    def testCreate1(self):
        # Test creation of options for new subtype scheme
        sc = SubtypeConfig(self.yamlfile)
        options = sc.create_subtype('test_gene',False)
        keys = ['alignment','subtype','lookup','seq']
        self.assertTrue(all(k in options for k in keys))


    def testCreate2(self):

        def touch(path):
            with open(path, 'a'):
                os.utime(path, None)

        # Test creation of directory for new subtype scheme
        scheme = 'test_gene'
        pre = SubtypeConfig(self.yamlfile)
        pre_options = pre.create_subtype(scheme,False)
        pre.save()
        
        # Create files
        paths = ['alignment','subtype','lookup']
        for p in paths:
            touch(os.path.join(self.root_dir, pre_options[p]))

        post = SubtypeConfig(self.yamlfile)
        post_options = post.get_subtype_config(scheme)

        self.assertEqual(pre_options, post_options)

   

def main():
    unittest.main()

if __name__ == '__main__':

    main()