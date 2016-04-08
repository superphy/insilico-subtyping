#!/usr/bin/env python

"""Phylotyper Config

Parse and validate Phylotyper config options

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

import ConfigParser
import os

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mwhiteside@canada.ca"


def which(program):
    """Checks if program is executable.

    Args:
        program (str): program path
        
    Returns:
        bool: True if successful, False otherwise

    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None



class PhylotyperOptions(object):
    """Store phylotyper options.

    Parse and validate INI config. Config INI file
    format is as follows:

    # External programs required in Phylotyper
    [external]
    fasttree=fastree_exe # (Required) Path to FastTree double-precision executable 
    
    Raises:
        Exception if required option is missing/invalid.

    """

    def __init__(self, inifile):
        """Constructor

        Parse config optionas from INI file
        
        Args:
            inifile (str): Location of INI config file
            
        """

        # if not inifile:
        #     raise Exception("Missing config file")

        self._inifile = inifile
        self._options = {}
        self._valid_option_names = []

        required = {
            'external': [
                'fasttree'
            ]
        }

        # optional = {
        #     'external': [
        #     ]
        # }

        config = ConfigParser.ConfigParser()
        config.read(inifile)

        for section in required.keys():
            # Define dict for section
            self._options[section] = {}

            # Load options for section
            for option in required[section]:
                option_name = "%s.%s" % (section,option)

                if not config.has_option(section, option):
                    raise Exception("Missing option in config file: %s." % option_name)

                value = config.get(section, option)

                if section == 'external':
                    # Make sure path is executable
                    if not which(value):
                        raise Exception("External program option not executable: %s" % option_name)

                self._options[section][option] = value
                self._valid_option_names.append(option_name)


    @property
    def ini_file(self):
        """str: Location of INI config file"""

        return self._inifile


    def _valid_option(self, section, option):

        option_name = "%s.%s" % (section,option)
        if option_name in self._valid_option_names:
            return True
        else:
            return False


    def get(self, section, option):
        """Option getter


        Args:
            section (str): INI section name
            option (str): option name

        Returns:
            option value

        """

        if not self._valid_option(section, option):
            raise Exception("Unknown config option: %s.%s" % (section.option))

        return self._options[section][option]

   