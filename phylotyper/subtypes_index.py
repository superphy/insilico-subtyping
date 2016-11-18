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

        for filepath_parameter in ['alignment', 'subtype', 'lookup']:

            if filepath_parameter in config:
                config[filepath_parameter] = os.path.join(root_dir, config[filepath_parameter])

                if not os.path.isfile(config[filepath_parameter]):
                    return "%s file %s not found" % (filepath_parameter, config[filepath_parameter])

            else:
                return "no <%s> parameter" % (filepath_parameter)

        if not config['seq'] and (config['seq'] != 'nt' or config['seq'] != 'aa'):
            return "missing/invalid <seq> parameter"


    def get_subtype_config(self, subtype_name):
        """Config getter

        Args:
            subtype_name (str): YAML subtype key

        Returns:
            dictionary with keys:
                alignment: Filepath
                subtype: Filepath
                lookup: Filepath
                seq: str (aa|nt)

        """

        if not self._config[subtype_name]:
            raise Exception("Unknown subtype: %s" % (subtype_name))

        return self._config[subtype_name]


    def create_subtype(self, scheme, is_aa=False):
        """Create subtype index file list

        Creates directory under root_dir named by <scheme>
        argument.  Generates filepaths for Phylotyper inputs.

        Note: subtype is not saved to the YAML file unless
        save() is called.

        Args:
            scheme (str): YAML subtype key
            is_aa (bool): True = amino acid sequences

        Returns:
            dictionary with keys:
                alignment: Filepath
                subtype: Filepath
                lookup: Filepath
                seq: str (aa|nt)

        Raises:
            Exception if <scheme> directory already exists

        """

        # Create directory
        newdir = os.path.join(self._root_dir, scheme)
        if os.path.exists(newdir):
            raise Exception('Directory {} exists. Cannot create new phylotyper subtype scheme {}.'.format(newdir, scheme))
        else:
            os.makedirs(newdir)

        # Create index options
        rel_paths = {
            'alignment': os.path.join(scheme, '{}_reference.affn'.format(scheme)),
            'subtype': os.path.join(scheme, '{}_subtypes.csv'.format(scheme)),
            'lookup': os.path.join(scheme, '{}_dictionary.json'.format(scheme)),
        }

        subtype_options = {
            'alignment': os.path.join(self._root_dir, rel_paths['alignment']),
            'subtype': os.path.join(self._root_dir, rel_paths['subtype']),
            'lookup': os.path.join(self._root_dir, rel_paths['lookup']),
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



        



        


   