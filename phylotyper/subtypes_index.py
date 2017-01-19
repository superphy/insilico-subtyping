#!/usr/bin/env python

"""Phylotyper Pre-built Subtypes Files

Provide paths to alignments and subtype files

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

.. _Google Python Style Guide:
   http://google.github.io/styleguide/pyguide.html
"""

import os
import yaml
from datetime import datetime

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2016, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


class SubtypeConfig(object):
    """Store paths for builtin subtype schemes.

    Paths are loaded from yaml file

    Raises:
        Exception if required option is missing/invalid.

    """

    def __init__(self, yamlfile):
        """Constructor

        Parse config options from YAML file
        
        Args:
            yamlfile (str): Location of YAML config file
            
        """
        
        if not yamlfile:
            raise Exception("Missing config file")

        self._yamlfile = yamlfile
        self._config = {}

        f = open(self._yamlfile)
        config_map = yaml.safe_load(f)
        f.close()

        self._rel_paths = config_map

        if not os.path.isabs( self._rel_paths['root_dir'] ):
            self._root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), config_map['root_dir']))
        else:
            self._root_dir = self._rel_paths['root_dir']

        if not os.path.isdir(self._root_dir):
            raise Exception("Invalid/missing <root_dir> config parameter in YAML file %s" % (yamlfile))

        if not 'subtypes' in config_map:
            raise Exception("Invalid/missing <subtypes> config parameter in YAML file %s" % (yamlfile))

        for subt in config_map['subtypes']:
            e = self._is_invalid_subtype_config(config_map['subtypes'][subt])
            if e:
                raise Exception("Invalid/missing config in YAML file %s for subtype %s: %s" % (yamlfile, subt, e))
            else:
                self._config[subt] = config_map['subtypes'][subt]
                


    def _is_invalid_subtype_config(self, config):
        """Check for required config and valid filepaths.
        self._root_dir will be prepended to all filepaths.


        Args:
            config (dict): Config dictionary

        Returns:
            None
                -or-
            str when error encountered

        """
        root_dir = self._root_dir

        for filepath_parameter in ['subtype', 'lookup', 'rate_matrix']:

            if filepath_parameter in config:
                config[filepath_parameter] = os.path.join(root_dir, config[filepath_parameter])

                if not os.path.isfile(config[filepath_parameter]):
                    return "%s file %s not found" % (filepath_parameter, config[filepath_parameter])

            else:
                return "no <%s> parameter" % (filepath_parameter)

        for filepaths_parameter in ['alignment']:

            if filepaths_parameter in config:
                fp_list = []
                for fp in config[filepaths_parameter]:

                    full_fp = os.path.join(root_dir, fp)
                    fp_list.append(full_fp)

                    if not os.path.isfile(full_fp):
                        return "%s file %s not found" % (filepaths_parameter, full_fp)

                config[filepaths_parameter] = fp_list

            else:
                return "no <%s> parameter" % (filepaths_parameter)

        if (not 'seq' in config) or not (config['seq'] != 'nt' or config['seq'] != 'aa'):
            return "missing/invalid <seq> parameter"

        if not 'nloci' in config or config['nloci'] < 1:
            return "missing/invalid <nloci> parameter"

        if not 'desc' in config:
            return "missing/invalid <desc> parameter"

        for blastdb_parameter in ['search_database']:

            if blastdb_parameter in config:
                config[blastdb_parameter] = os.path.join(root_dir, config[blastdb_parameter])

                blast_file = config[blastdb_parameter]+'.psq'
                if not os.path.isfile(blast_file):
                    return "%s blast database %s not found" % (blastdb_parameter, config[blastdb_parameter])

            else:
                return "no <%s> parameter" % (blastdb_parameter)



    def get_subtype_config(self, subtype_name):
        """Config getter

        Args:
            subtype_name (str): YAML subtype key

        Returns:
            dictionary with keys:
                alignment: list of filepath
                subtype: Filepath
                lookup: Filepath
                search_database:
                seq: str (aa|nt)
                rate_matrix: Filepath

        """

        if not self._config[subtype_name]:
            raise Exception("Unknown subtype: %s" % (subtype_name))

        return self._config[subtype_name]


    def create_subtype(self, scheme, num_loci, is_aa=False, description='No description available'):
        """Create subtype index file list

        Creates directory under root_dir named by <scheme>
        argument.  Generates filepaths for Phylotyper inputs.

        Note: subtype is not saved to the YAML file unless
        save() is called.

        Args:
            scheme (str): YAML subtype key
            num_loci (int): Number if loci/alignments in this scheme
            is_aa (bool): True = amino acid sequences
            description (str): Description of subtype scheme for help info

        Returns:
            dictionary with keys:
                alignment: list of filepaths
                subtype: Filepath
                lookup: Filepath
                search_database: Filepath
                seq: str (aa|nt)
                rate_matrix: Filepath
                nloci: int

        """

        # Create directory
        newdir = os.path.join(self._root_dir, scheme)
        if not os.path.exists(newdir):
            os.makedirs(newdir)

        # Create index options
        alignment_file_list = []
        for i in xrange(num_loci):
            alignment_file_list.append(os.path.join(scheme, '{}_reference{}.affn'.format(scheme, i)))

        print description

        rel_paths = {
            'alignment': alignment_file_list,
            'subtype': os.path.join(scheme, '{}_subtypes.csv'.format(scheme)),
            'lookup': os.path.join(scheme, '{}_dictionary.json'.format(scheme)),
            'search_database': os.path.join(scheme, '{}_search_database'.format(scheme)),
            'rate_matrix': os.path.join(scheme, '{}_rate_matrix.rds'.format(scheme)),
            'nloci': num_loci,
            'desc': description
        }

        subtype_options = {
            'alignment': [ os.path.join(self._root_dir, f) for f in rel_paths['alignment'] ],
            'subtype': os.path.join(self._root_dir, rel_paths['subtype']),
            'lookup': os.path.join(self._root_dir, rel_paths['lookup']),
            'search_database': os.path.join(self._root_dir, rel_paths['search_database']),
            'rate_matrix': os.path.join(self._root_dir, rel_paths['rate_matrix']),
            'nloci': num_loci,
            'desc': description
        }

        # Sequence type
        if is_aa:
            subtype_options['seq'] = 'aa'
            rel_paths['seq'] = 'aa'
        else:
            subtype_options['seq'] = 'nt'
            rel_paths['seq'] = 'nt'

        # Store in memory
        self._config[scheme] = subtype_options
        self._rel_paths['subtypes'][scheme] = rel_paths

        return subtype_options


    def save(self):
        """Save current options to YAML file

        Creates backup before overwritting file

        Args:
            None

        """

        # Backup existing file
        yamlfile = self._yamlfile
        modified_time = os.path.getmtime(yamlfile) 

        time_stamp =  datetime.fromtimestamp(modified_time).strftime("%b-%d-%y-%H:%M:%S")
        os.rename(yamlfile, yamlfile+"_"+time_stamp)

        output = '{}\n\n{}\n'.format('# Phylotyper filepaths and options for pre-built subtyping schemes',
            yaml.dump(self._rel_paths))

        with open(self._yamlfile, 'w') as w:
            w.write(output)
            w.flush()


    def list(self):
        """List available subtypes

        Args:
            None

        Returns:
            tuple containing:
                0: list of subtype scheme names
                1: nt|aa status
                2: number of loci
                3: description

        """

        names = self._config.keys()
        types = []
        nloci = []
        desc = []
        for n in names:
            opt = self._config[n]
            types.append(opt['seq'])
            nloci.append(opt['nloci'])
            desc.append(opt['desc'])

        return (names, types, nloci, desc)



        



        


   