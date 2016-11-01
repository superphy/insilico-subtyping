#!/usr/bin/env python

"""Phylotyper Built-in Subtypes

For subtypes packaged in phylotyper, provide paths to alignments and subtype
files

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


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
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


        self._root_dir = os.path.abspath(config_map['root_dir'])
        if not os.path.isdir(self._root_dir):
            raise Exception("Invalid/missing <root_dir> config parameter in YAML file %s" % (yamlfile))

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

        """

        if not self._config[subtype_name]:
            raise Exception("Unknown subtype: %s" % (subtype_name))

        return self._config[subtype_name]


    def add_subtype(self, scheme, options):
        pass

   