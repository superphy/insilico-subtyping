#!/usr/bin/env python

"""Phylotyper module

Phylotyper functions

"""

import logging
import os
import pkg_resources
import rpy2.robjects as robjects
from subprocess import check_output, CalledProcessError, STDOUT


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


class Phylotyper(object):
    """Phylotyper functions.


    """

    def __init__(self, config):
        """Constructor

        Initializes R environment and script settings. Initializes logger.
        
        Args:
            config (obj): PhylotyperConfig with .ini file settings
            
        """

        self.logger = logging.getLogger('phylotyper.phylotyper.Phylotyper')
        self.rwd = pkg_resources.resource_filename(__name__,  '/'.join(('R')))
        self.rlib = config.get('R', 'lib')
        self.rrepo = config.get('R', 'repo')
        self.rscript = config.get('R', 'rscript')
        self.rfiles = {
            'phylotyper': os.path.join(self.rwd, 'phylotyper.R'),
            'tests': os.path.join(self.rwd, 'test_functions.R'),
            'performance_tests': os.path.join(self.rwd, 'performance_tests.R')
        }


    def subtype(self, tree, subtypes, output_dir, plot_name='posterior_probability_tree.png'):
    	"""Run phylotyper subtyping method

        Wrapper around the phylotyper.R functions. Calls
        R functions using rpy2 interface. Produces png showing 
        posterior probability as pie charts on tree.
        
        Args:
            tree (str): Filepath to newick tree containing both typed & untyped genes
            subtypes (str): Filepath to tab-delimited subtype assignments for typed genes
            output_dir (str): Filepath to output directory
            noplots (bool): Generate png images of tree

        Returns:
            assignments (dict): Each untyped gene key will contain a tuple with:
                1. One or more predicted subtypes
                2. Posterior probability for that assignment
            
        """

        # Set work dir
        robjects.r('setwd("%s")' % self.rwd)

        # Load Phylotyper R functions
        rcode = 'source("%s", chdir=TRUE)' % (self.rfiles['phylotyper']) 
        robjects.r(rcode)

        # Load libraries
        rcode = 'suppressMessages(phylotyper$loadInstallLibraries(libloc="%s",repo="%s"))' % (self.rlib, self.rrepo) 
        robjects.r(rcode)

        # Load data files
        rcode = 'rs = phylotyper$loadSubtype("%s","%s"); tree = rs$tree; subtypes = rs$subtypes; untyped = rs$untyped' % (tree, subtypes)
        robjects.r(rcode)

        untyped = robjects.r('untyped')
        self.logger.debug('phylotyper.R detected the following genes have no subtype: %s' % (','.join(untyped)))

        # Make prior matrix
        rcode = 'priorR = phylotyper$makePriors(tree, subtypes); priorM = priorR$prior.matrix'
        robjects.r(rcode)
       
        # Run subtype procedure
        rcode = 'est.scheme = 5; result = phylotyper$runSubtypeProcedure(tree, priorM, est.scheme, tips=untyped)'
        robjects.r(rcode)

        # Write subtype predictions
        rcode = '''
        cn = colnames(result$tip.pp);
        matrix(round(result$tip.pp[untyped,], digits=7),ncol=length(cn),nrow=length(untyped),dimnames=list(untyped,cn), byrow=TRUE)
        '''
        predictions = robjects.r(rcode)
        subtype_states = predictions.colnames
        assignments = {}

        # Find largest pp
        for genome in untyped:
            pprow = predictions.rx(genome, True)
            
            maxpp = -1
            maxst = []
            i = 0
            for pp in pprow:
                if pp > maxpp:
                    maxst = [subtype_states[i]]
                    maxpp = pp
                elif pp == maxpp:
                    maxst.append(subtype_states[i])

                i += 1

            assignments[genome] = (maxst, maxpp)

        # Plot pp on the tree
        if plot_name:
            rcode = '''
            dim = phylotyper$plotDim(tree)
            graphics.off()
            png(filename="%s", width=dim[['x']],height=dim[['y']],res=dim[['res']])
            do.call(result$plot.function, list(tree=tree, fit=result$result, subtypes=subtypes))
            graphics.off()
            ''' % (plot_name)
            robjects.r(rcode)

        # Return assignments
        return assignments
  		

    def evaluate(self, tree, subtypes, output_dir):
        """Run phylotyper performance tests
        
        Calls external R script
        
        Args:
            tree (str): Filepath to newick tree containing both typed & untyped genes
            subtypes (str): Filepath to tab-delimited subtype assignments for typed genes
            output_dir (str): Filepath to output directory

        Returns:
            results are output to output_dir
            
        """

        cmd = "{} {} {} {} --out {}".format(self.rscript, self.rfiles['performance_tests'],
            tree, subtypes, output_dir)
        cwd = os.getcwd()
        os.chdir(self.rwd)

        try:
            check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
        except CalledProcessError as e:
            os.chdir(cwd)
            msg = "R script {} failed: {} (return code: {}).".format(cmd, e.output, e.returncode)                                                                                                   
            raise Exception(msg)

        os.chdir(cwd)
        return None





        


       


  
	