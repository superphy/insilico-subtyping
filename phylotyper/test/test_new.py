import os
import re
import shutil
import tempfile
import unittest
#import yaml

from phylotyper import subtypes_index

class NewTests(unittest.TestCase):

    def testNew(self):
        self.failUnless(False)


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
        #shutil.rmtree(self.root_dir)
        pass


    def testCreate1(self):
        # Test creation of options for new subtype scheme
        sc = subtypes_index.SubtypeConfig(self.yamlfile)
        options = sc.create_subtype('test_gene',False)
        keys = ['alignment','subtype','lookup','seq']
        self.assertTrue(all(k in options for k in keys))


    def testCreate2(self):

        def touch(path):
            with open(path, 'a'):
                os.utime(path, None)

        # Test creation of directory for new subtype scheme
        scheme = 'test_gene'
        pre = subtypes_index.SubtypeConfig(self.yamlfile)
        pre_options = pre.create_subtype(scheme,False)
        pre.save()
        
        # Create files
        paths = ['alignment','subtype','lookup']
        for p in paths:
            touch(os.path.join(self.root_dir, pre_options[p]))

        post = subtypes_index.SubtypeConfig(self.yamlfile)
        post_options = post.get_subtype_config(scheme)

        self.assertEqual(pre_options, post_options)

   

def main():
    unittest.main()

if __name__ == '__main__':
    main()