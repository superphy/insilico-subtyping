#!/usr/bin/env python

"""Bio Sequences Utilities

Tools for manipulating / validating input DNA and protein sequences


"""

import csv
import hashlib
import json
import re
import warnings
from Bio import SeqIO
from collections import Counter, defaultdict, namedtuple

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mwhiteside@canada.ca"


class SeqDict(object):
    """Library of non-redundant set of sequences for quick lookup



    """

    def __init__(self, n_loci=1, format_name=True):
        """Constructor

        Args:
            None

        """

        self._seqdict = {}
        self._genenum = 0
        self._hash_algorithm = 'md5'
        self._nloci = n_loci
        self._unique_names = set()
        self._format_names = format_name


    @property
    def seqs(self):
        """dict: ."""

        return self._seqdict

    @property
    def num(self):
        """Number sequences"""

        return self._genenum

   
    def build(self, fasta_files, subtype_file):
        """Populate SeqDict with entries in fasta and subtype files

        Args:
            fasta_file (str|list): Filepath to input sequence file or list of such files
            subtype_file (str): Filepath to input subtype

        Returns:
            None

        """

        # Check subtype file
        subtypes = {}

        for row in csv.reader(open(subtype_file,'r'),delimiter='\t'):
            name = row[0]
            subt = row[1]

            subtypes[name] = subt.lower()
            
        if isinstance(fasta_files, str):
            # Create list
            fasta_files = [fasta_files]

        if len(fasta_files) != self._nloci:
            raise Exception("Missing fasta file. {} fasta files provided for {} number of loci.".format(len(fasta_files), self._nloci))

        concat = LociConcat()
        sequences = concat.collect(fasta_files)
      
        for name,seqslist in sequences.iteritems():
            this_subt = subtypes[name]

            for seqs in seqslist:
                self.add(seqs, name, this_subt)


    def find(self, seq):
        """Search for identical sequence in dictionary
        
        Args:
            seq (str|list): sequence or list of sequences

        Returns:
            None if no identical sequence found -or- SeqDict entry
            
        """
        
        if not isinstance(seq, str):
            seqstr = ''
            # Concatenate list
            for s in seq:
                seqstr += self.prepseq(s.upper())
        else:
            seqstr = self.prepseq(seq.upper())
        
        searchstr = self.digest(seqstr)

        if searchstr in self.seqs:
            hits = self.seqs[searchstr]

            if seqstr in hits:
                return hits[seqstr]
            else:
                return None
            
        else:
            return None
            

    def prepseq(self, seq):
        """Sequence manipulations

        Ignore trailing stop codons, if any

        """

        wtf = re.sub(r'\*$', '', seq)
        return wtf

   
    def add(self, seq, name, subt):
        """Add new entry to SecDict
        
        Args:
            seq (list|str): list of sequences or single sequence if one loci
            name (str): accession
            subt (str): subtype

        Returns:
            None
            
        """
        
        if self._nloci > 1:
            
            if isinstance(seq, str):
                raise Exception('Invalid seq parameter. Multiple loci need to be passed as list')
           
            keystr = ''
            for s in xrange(len(seq)):
                thisseq = seq[s]
                keystr += self.prepseq(thisseq.upper())

        else:
            if not isinstance(seq, str):
                keystr = seq[0].upper()
            keystr = self.prepseq(keystr)

        keystr = keystr.replace('-','')
        matched = self.find(keystr)

        if matched:
            # Existing identical sequence
            # Append accession
            if not matched['subtype'] == subt:
                raise Exception('Identical sequence {} has different subtype assignment than {} (expected: {}, assigned: {})'.format(
                    name, ','.join(matched['accessions']), matched['subtype'], subt))

            # Record gene/genome
            matched['accessions'].append(name)

        else:
            # Add new sequence
            searchstr = self.digest(keystr)
            if not searchstr in self.seqs:
                self.seqs[searchstr] = {}

            self._genenum += 1

            # Set up sequence name
            this_name = self.format_name(name, subt)
            if not self._format_names:
                if name in self._unique_names:
                    raise Exception("Non-unique name for distict sequences {}".format(name))
                else:
                    this_name = name
                    self._unique_names.add(name)

            self.seqs[searchstr][keystr] = {
                'subtype': subt,
                'name': this_name,
                'accessions': [name]
            }

            if self._nloci > 1:
                self.seqs[searchstr][keystr]['loci'] = seq


    def accession_map(self):
        """Return dict of SeqDict names mapped to lists of original accessions"""

        acc_map = {}
        for seqkey,seqs in self._seqdict.iteritems():
            for seq,seqentry in seqs.iteritems():
                acc_map[seqentry['name']] = seqentry['accessions']

        return acc_map


    def format_name(self, name, subtype):
        """Standard gene naming"""

        # Try to isolate NCBI accession
        nm = re.sub(r'^(?:\d+\-)?([A-Z0-9\.]+)_.+$', r'\1', name.upper())
        genename = '{}-{}__{}'.format(self._genenum, nm, subtype.lower())

        return genename


    def store(self, reffile):
        """Save SeqDict to file"""
        with open(reffile, 'w') as rfh:
            json.dump(self.seqs,rfh)

        return None


    def load(self, reffile):
        """Load SeqDict from file"""
        with open(reffile, 'r') as rfh:
            self._seqdict = json.load(rfh)

        # Check format
        keys =  ['name','subtype','accessions']
        if self._nloci > 1:
            keys.append('loci')
        for seqkey,seqs in self._seqdict.iteritems():
            for seq,seqentry in seqs.iteritems():
                self._genenum += 1

                for k in keys:
                    if not k in seqentry:
                        raise Exception('Improperly formated SeqDict object')

                if self._nloci > 1:
                    if len(seqentry['loci']) != self._nloci:
                        raise Exception('Improperly formated SeqDict object')

        return None


    def digest(self, seq):
        """Compute quick lookup hash value from sequence"""

        h = hashlib.new(self._hash_algorithm)
        h.update(seq)
        dig = h.hexdigest()

        return dig


    def write(self, fasta_filepaths, subtype_filepath):
        """Write unique sequences in fasta format
        Args:
            fasta_filepaths (list): List of filepaths to output fasta files
            subtype_filepath (str): Filepath to output subtype file
            
        Returns:
            None
        """

        if self._nloci == 1 and isinstance(fasta_filepaths, str):
            fasta_filepaths = [fasta_filepaths]

        if len(fasta_filepaths) != self._nloci:
            raise Exception("Missing fasta file. {} fasta files provided for {} loci.".format(len(fasta_filepaths), self._nloci))


        # Open all needed files
        try:
            fhs = [ open(filename, 'w') for filename in fasta_filepaths ]

            with open(subtype_filepath, 'w') as s:
                for seqkey,seqs in self.seqs.iteritems():
                    for seq,seqentry in seqs.iteritems():

                        if self._nloci == 1:
                            f = fhs[0]
                            f.write(">{}\n{}\n".format(seqentry['name'], seq))
                            
                        else:
                            for f,loci_seq in zip(fhs, seqentry['loci']):
                                f.write(">{}\n{}\n".format(seqentry['name'], loci_seq))
                               
                    s.write("{}\t{}\n".format(seqentry['name'], seqentry['subtype']))


        except IOError as e:
            for f in fhs:
                f.close()

            raise e


        return None


    def subtype_occurences(self):
        """Count occurences of each subtype in SeqDict

        Returns:
            Counter object

        """

        subtype_counts = Counter()

        for seqkey,seqs in self.seqs.iteritems():
            for seq,seqentry in seqs.iteritems():

                subtype_counts[seqentry['subtype']] += 1

        return subtype_counts



class LociConcat(object):
    """Concatenate loci from multiple fasta input files into supersequences.

    Occaisonly multiple loci are used to predict traits.  These can be concatenated
    to form supersequences.

    """
    
    def safe_collapse(self, inputs, fasta_filepath=None):
        """Concatenate loci from mulitple fasta files

        Returns dict indexed by genome ID containing supersequences.
        Optionally, writes to fasta file.

        Assumes that there is one supersequence per genome.  Checks
        this assertion. This ususually guranteed by manipulations to the
        fasta headers/files

        Args:
            inputs (list): Filepaths to fasta sequences.
            fasta_filepath (str)[Optional]: If provided, write supersequences to this file
        
        Returns:
            dict

        Raise:
            Exception when missing loci in a file
            Duplicate alleles per genome

        """

        seqdict = {}
        superseqs, incomplete = self.load(inputs, parseheaders=False)
        for genome in superseqs:
            seqlist = [ a.seq() for a in superseqs[genome].iteralleles() ]
            if len(seqlist) > 1:
                raise Exception("Multiple alleles for genome {}".format(genome))
            seqdict[genome] = seqlist[0]
           
        if fasta_filepath:
            with open(fasta_filepath, 'w') as f:
                for name,seq in seqdict.iteritems():
                    f.write(">{}\n{}\n".format(name, seq))
                
        return seqdict


    def collect(self, inputs):
        """Concatenate loci from mulitple fasta files

        Returns dict indexed by fasta ID containg list of loci sequences
       
        Args:
            inputs (list): Filepaths to fasta sequences.
        
        Returns:
            dict

        Raises:
            Exception when missing loci in a file

        """

        sequencesonly = {}
        superseqs, incomplete = self.load(inputs)
        for genome in superseqs:
            sequencesonly[genome] = [ a.seqlist() for a in superseqs[genome].iteralleles() ]
        
        return sequencesonly


    def load(self, inputs, parseheaders=True, missing='raise'):
        """Read loci from mulitple fasta files

        Returns dict indexed by fasta ID containg list of loci sequences.
        Tracks all individual allele fasta headers and sequences.
       
        Args:
            inputs (list): Filepaths to fasta sequences.
            parseheaders (bool): Split header into genome|allele parts
            missing (str): Behavior when genome is missing one loci. Options are 'raise' Exception or
              'warn'.
        
        Returns:
            dict of typing sequences, set of genomes with missing loci

        Raise:
            Exception when missing loci in a file

        """

        i = 0
        
        sequences = defaultdict(TypingSequence)
        uniq = Counter()
        genomes = set()
        incomplete = set()

        # Load gene allele sequences into memory
        # for each loci
        for ff in inputs:
            
            fasta = SeqIO.parse(ff, 'fasta')
            for record in fasta:
                if parseheaders:
                    name = LociConcat.parse_fasta_id(record.id)
                    genome = name.genome
                else:
                    name = genome = record.id

               
                uniq[genome] = i

                if i == 0:
                    genomes.add(genome)
                else:
                    if not genome in genomes:
                        if missing == 'raise':
                            raise Exception('Unknown fasta record for {} in file {}'.format(name, ff))
                        elif missing == 'warn':
                            warnings.warn('Loci missing for {}.'.format(name))
                        incomplete.add(genome)

                # Record allele for this genome
                sequences[genome].add(i, str(record.seq), str(record.description))
                print('Added: {} in position {}'.format(str(record.description), i))

            for name in uniq:
                if uniq[name] != i:
                    if missing == 'raise':
                        raise Exception('Missing fasta record for {} in file {}'.format(name, ff))
                    elif missing == 'warn':
                        warnings.warn('Loci missing for {}.'.format(name))
                    incomplete.add(name)

            i += 1

        return (sequences, incomplete)
        

    # def combinations(self, loci_sets, loci):
    #     """Generate all combinations of alleles for each of the loci sets

    #     Returns dict indexed by genome name containg list of tuples. Each
    #     tuple represents one allele combination. Tuples contain two lists:
    #     1. Names of all concatenated loci in this combination
    #     2. Sequences of all concatenated loci in this combination

    #     Args:
    #         loci_sets (list): Names of loci sets
    #         loci (dict): dict of genome alleles 
        
    #     Returns:
    #         dict

    #     Raise:
    #         Exception when missing loci in a file

    #     """
    #     # 

    #     ls = loci_sets.pop()

    #     if not loci_sets:
    #         return loci[ls].iteritems()

    #     else:
    #         seqs = self.combinations(loci_sets, loci)
    #         newseqs = defaultdict(list)
    #         for genome, hits in loci[ls].iteritems():

    #             # Genomes should have at least one copy in each loci set
    #             if not genome in seqs:
    #                 raise Exception('Genome {} missing loci {}'.format(genome, ls))

    #             # Append all copies for this loci
    #             for identifier, seq in hits.iteritems():
    #                 for s in seqs[genome]:
    #                     newseqs[genome].append((identifier+s[0],seq+s[1]))

    #         return newseqs


    @classmethod
    def parse_fasta_id(cls, id):
        """Extract Genome and gene identifier fasta header

        From the ID component of a fasta header, determine
        for a particular genome copy, whether there is multiple 
        alleles of a gene in the genome.  Format is:

            >[lcl|]genome|allele

            Args:
                id (str): Fasta ID

            Returns:
                AlleleID namedtuple

        """

        # Remove database id
        id = re.sub(r'^(?:lcl\|)|(?:gi\|)', '', id)

        # Extract genome and gene Ids if they exist
        parts = re.split(r'\|', id, maxsplit=1)
        
        return AlleleID(*parts)



# Factory for AlleleID class
# Stores genome and allele ID components
class AlleleID(namedtuple('AlleleID', ['genome', 'allele'])):
    __slots__ = ()
    def __new__(cls, genome, allele=None):
        # add default values
        return super(AlleleID, cls).__new__(cls, genome, allele)
    def __str__(self):
        # String
        allele = '|'+str(self.allele) if self.allele else ''
        return '{}{}'.format(self.genome, allele)



class TypingSequence(object):
    """Stores supersequences for single genomes 

    In phylotyper, multiple loci can be concatenated to form single typing sequence.
    Genomes can have multiple copy of any one loci.  TypingSequence is a storage class
    for managing this data for a genome.

    """

    def __init__(self):
        """Constructor

        Args:
            None

        """

        self.alleles = []
        self.positions = {}
        self.nloci = 0
        self.typing_sequences = [[]]
        self.combinations = []
       

    def __iter__(self):
        return iter(self.typing_sequences)


    def __len__(self):
        return len(self.typing_sequences)


    def add(self, loci, sequence, header):
        """Add allele for given loci

        Either appends new loci to typing sequence,
        or if this loci is already in sequence, splits
        typing sequence into multiple instances for each
        allele.

        Args:
            loci(str): Loci identifier
            sequence(str): bio sequence string
            header(str): Fasta header associated with this loci copy

        Returns:
            None

        """

        alleleid = len(self.alleles)
        newallele = Allele(alleleid, sequence, header)
        self.alleles.append(newallele)

        print(self.typing_sequences)
        print(self.combinations)

        if not isinstance(loci, str):
            loci = str(loci)

        if not loci in self.positions:
            # Adding new loci

            self.positions[loci] = self.nloci

            for ts in self.typing_sequences:
                ts.append(alleleid)

            if self.nloci < len(self.combinations):
                self.combinations[self.nloci].append(alleleid)
            else:
                self.combinations.append([alleleid])

            self.nloci += 1
            print('extending ts')

        else:
            # Add allele to existing loci
            pos = self.positions[loci]
            newseqs = []

            # Generate all combinations
            for c in xrange(len(self.combinations)):
                print(c)
                print(newseqs)
                if c == pos:
                    if newseqs:
                        for s in newseqs:
                            copy = self._copy(s)
                            print(copy)
                            print(c)
                            copy[c] = alleleid
                            newseqs.append(copy)
                    else:
                        newseqs.append([alleleid])

                else:
                    for a in self.combinations[c]:
                        if newseqs:
                            for s in newseqs:
                                copy = self._copy(s)
                                copy[c] = a
                                newseqs.append(copy)
                        else:
                            newseqs.append([a])
                            
            print('Adding variations:')
            print(newseqs)
            # for ts in self.typing_sequences:
            #     # make copy and change allele
            #     copy = self._copy(ts)
                
            #     newseqs.append(copy)

            # Save new typing sequences with alternate
            # allele
            self.typing_sequences.extend(newseqs)
            print(self.typing_sequences)


        print(len(self.typing_sequences))

    def _copy(self, ts):
        # Make a copy of a typing sequence

        newts = []
        for allele in ts:
            newts.append(allele)

        return newts


    def iteralleles(self):
        """generator iterating over each typing sequence

        Each iteration returns an AlleleList object

        """

        for ts in self.typing_sequences:
            yield AlleleList([ self.alleles[a] for a in ts ])


# Helper class for working with lists of Alleles
class AlleleList(object):

    def __init__(self, allele_list):
        self.allele_list = allele_list

    def seqlist(self):
        return [ a.sequence for a in self.allele_list ]

    def seq(self):
        return ''.join(self.seqlist())

    def idlist(self):
        return [ a.allele for a in self.allele_list ]

    def indexlist(self):
        return [ a.alleleid for a in self.allele_list ]

    def iddump(self):
        return json.dumps(self.idlist())

    def alleles(self):
        return self.allele_list


# Factory for Allele class
# Stores sequence and allele ID components
class Allele(namedtuple('Allele', ['alleleid', 'sequence', 'allele'])):
    __slots__ = ()
    def __str__(self):
        # String
        return '{}'.format(self.sequence)
