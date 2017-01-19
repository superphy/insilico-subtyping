import os
import unittest

from phylotyper.tree.seq import LociConcat

class SequenceManipulationTests(unittest.TestCase):

    def setUp(self):

        # Inputs
        self.data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'data'))

    def testLociConcat(self):

        concat = LociConcat()
        inputs = []
        for i in ('locus1.fasta','locus2.fasta','locus3.fasta'):
            inputs.append(os.path.join(self.data_dir, i))
        sequences = concat.load(inputs)[0]

	print sequences

        correctlen = [ len(ts) == 3 for ts in sequences['genome1'] ]
       
        self.assertTrue(all([
            len(sequences.keys()) == 4, 
            sequences['genome2'].nloci == 3,
            all(correctlen),
            len(sequences['genome3']) == 2,
            len(sequences['genome1']) == 4
        ]))

    def testLociConcat2(self):

        concat = LociConcat()
        inputs = []
        for i in ('locus1.fasta','locus2.fasta','locus3.fasta'):
            inputs.append(os.path.join(self.data_dir, i))
        sequences = concat.load(inputs)[0]

        it = sequences['genome1'].iteralleles()
        testseq1 = it.next()
        testseq2 = it.next()
        testseq3 = it.next()
        testseq4 = it.next()

        self.assertTrue(all([
            testseq1.seq() == 'ATGCTACGAA',
            testseq2.seq() == 'ATGGTACGAA',
            testseq3.seq() == 'ATGCTACGAT',
            testseq4.seq() == 'ATGGTACGAT',
        ]))

   


def main():
    unittest.main()

if __name__ == '__main__':

    main()
