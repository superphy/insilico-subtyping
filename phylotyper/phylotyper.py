#!/usr/bin/env python

"""Phylotyper module

Phylotyper functions


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

import logging


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


class Phylotyper(object):
    """Phylotyper functions.

    
    Raises:
        Exception if required option is missing/invalid.

    """

    def __init__(self):
        """Constructor

        Description
        
        Args:
            something (str): Description
            
        """

        self.logger = logging.getLogger('phylotyper.phylotyper.Phylotyper')

  
  	def subtype(self, tree, subtypes, unknowns):
  		"""Constructor

        Description
        
        Args:
            something (str): Description
            
        """
  		pass




  
	