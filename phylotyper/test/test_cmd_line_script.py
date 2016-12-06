import contextlib
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

def spawn_python(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT):
    """Run a Python subprocess with the given arguments.

        Returns a Popen object.
    """
    cmd_line = [sys.executable]
    cmd_line.extend(args)
    
    return subprocess.Popen(cmd_line, stdin=subprocess.PIPE,
                        stdout=stdout, stderr=stderr)

def kill_python(p):
    """Run the given Popen process until completion and return stdout."""
    p.stdin.close()
    data = p.stdout.read()
    p.stdout.close()
   
    return data


class ScriptTest(unittest.TestCase):

    def setUp(self):
        # Create temporary directory
        self.abspath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

        sandbox = os.path.abspath(os.path.join(os.path.dirname(__file__),'sandbox'))
        self.test_dir = os.path.abspath(tempfile.mkdtemp(dir=sandbox))
        
        # Inputs
        self.data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'data'))
        aa = True
        suffix = 'a' if aa else 'n'
        ref = os.path.join(self.data_dir, 'stx2.ff'+suffix)
        subt = os.path.join(self.data_dir, 'stx2_subtypes.csv')
        inpu = os.path.join(self.data_dir, 'test_stx2.ff'+suffix)
        scheme = 'test_gene'

        if not os.environ.get('PHYLOTYPER_CONFIG'):
            msg = 'Missing config file enviroment variable.\nMust provide Phylotyper config file using' \
                ' enviroment variable PHYLOTYPER_CONFIG for testing'
            raise Exception(msg)


        


    def tearDown(self):
        # Remove previous directories created
        shutil.rmtree(self.test_dir)

    @contextlib.contextmanager
    def run_script(self, args, separate_stderr=False):
        if separate_stderr:
            p = spawn_python(args, stderr=subprocess.PIPE)
            stderr = p.stderr
        else:
            p = spawn_python(args, stderr=subprocess.STDOUT)
            stderr = p.stdout

        try:
            yield p
        finally:
            kill_python(p)
            stderr.close()

    
    def test_missing_arguments(self):

        cmd = [os.path.join(self.abspath, 'search.py')]
        with self.run_script(cmd, separate_stderr=True) as p:
            response = p.stderr.readline().strip()
            self.assertRegexpMatches(response, r'usage: search.py')


