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
        self.rfiles = {
            'phylotyper': os.path.join(self.rwd, 'phylotyper.R'),
            'tests': os.path.join(self.rwd, 'test_functions.R')
        }


    def subtype(self, tree, subtypes, output_dir):
    	"""Run phylotyper subtyping method

        Wrapper around the phylotyper.R functions. Calls
        R functions using rpy2 interface. Produces png showing 
        posterior probability as pie charts on tree.
        
        Args:
            tree (str): Filepath to newick tree containing both typed & untyped genes
            subtypes (str): Filepath to tab-delimited subtype assignments for typed genes
            output_dir (str): Filepath to output directory

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
        rcode = 'priorR = phylotyper$makePriors(tree, subtypes); priorM = priorM = priorR$prior.matrix'
        robjects.r(rcode)
       
        # Run subtype procedure
        rcode = 'est.scheme = 1; result = phylotyper$runSubtypeProcedure(tree, priorM, est.scheme)'
        robjects.r(rcode)

        # Write subtype predictions
        predictions = robjects.r('round(result$tip.pp[untyped,], digits=7)')
        subtype_states = predictions.colnames
        assignments = {}

        # Find largest pp over thresold
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
        rcode = '''
        png(filename=file.path("%s", "posterior_probability_tree.png"), width=1024, height=1024)
        do.call(result$plot.function, list(tree=tree, fit=result$result, subtypes=subtypes))
        dev.off()
        ''' % (output_dir)
        robjects.r(rcode)

        # Return assignments
        return assignments
  		

    def evaluate(self, tree, subtypes, output_dir):
        """Run phylotyper performance tests

        Calls functions in test_functions.R. Produces png showing 
        posterior probability as pie charts on tree, as well as accuracy
        metrics.
        
        Args:
            tree (str): Filepath to newick tree containing both typed & untyped genes
            subtypes (str): Filepath to tab-delimited subtype assignments for typed genes
            output_dir (str): Filepath to output directory

        Returns:
            results are output to output_dir
            
        """

        # Set work dir
        robjects.r('setwd("%s")' % self.rwd)
        
        # Load libraries and source
        rcode = 'source("{}")'.format(self.rfiles['phylotyper']) 
        robjects.r(rcode)

        rcode = 'suppressMessages(phylotyper$loadInstallLibraries(libloc="%s",repo="%s"))' % (self.rlib, self.rrepo) 
        robjects.r(rcode)

        rcode = 'source("{}")'.format(self.rfiles['tests']) 
        robjects.r(rcode)

        # Set output directory
        rcode = 'output_dir = "{}"'.format(output_dir)
        robjects.r(rcode)

        # Load data files
        rcode = 'rs = phylotyper$loadSubtype("{}","{}"); tree = rs$tree; subtypes = rs$subtypes; untyped = rs$untyped'.format(tree, subtypes)
        robjects.r(rcode)

        rcode = '''
        # Run Mk model evaluation (this markov model is used in simmap and rerootingMethod)
        aic = evaluateModels(tree,subtypes)
        file = 'model_aic'
        write.table(aic, file=file.path(output_dir, paste(file, '.csv', sep='')),
            sep="\t",
            quote=FALSE)

        # Iterate through esimtation procedures
        estimation.methods = list(rerooting=1)
        for(i in 1:length(estimation.methods)) {
            est.scheme = estimation.methods[[i]]
            est.name = names(estimation.methods)[i]

            # Overlay posterior probabilities in tree plot
            print("Running estimation procedure: ", est.name)

            priorR = phylotyper$makePriors(tree, subtypes)
            priorM = priorR$prior.matrix
            result = phylotyper$runSubtypeProcedure(tree, priorM, est.scheme)
            file = 'posterior_probability_tree'
            dim = phylotyper$plotDim(tree)
            png(filename=file.path(output_dir, paste(est.name, '_', file, '.png', sep='')),
                width=dim[['x']],height=dim[['y']],res=dim[['res']]
            )
            do.call(result$plot.function, list(tree=tree, fit=result$result, subtypes=subtypes))
            dev.off()

            # Iterate through validation procedures
            # Leave-One-Out CV
            # 5-fold CV
            for(validation in c('loocv', 'kfcv')) {
                print("Running validation: ", validation)

                pp = do.call(validation, c(tree=tree, subtypes=subtypes, scheme=est.scheme))

                # Summarize performance
                results = simulationSummary(subtypes, pp)

                # Write performance metrics to file
                file = 'performance_metrics'
                write.table(results$metrics, file=file.path(output_dir, paste(est.name, '_', validation, '_', file, '.csv', sep='')),
                    sep="\t",
                    quote=FALSE)


                # Plots

                # Plot confusion matrix
                file = 'confusion_matrix'
                png(filename=file.path(output_dir, paste(est.name, '_', validation, '_', file, '.png', sep='')))
                plotConfusionMatrix(results$confusion.matrix)
                dev.off()

                # Plot posterior probability histogram
                file = 'posterior_probability_histogram'
                png(filename=file.path(output_dir, paste(est.name, '_', validation, '_', file, '.png', sep='')))
                plotPPHistogram(results$test.results, subtypes)
                dev.off()

                print(validation, "complete")
            }

            print(est.name, "complete")
        }
        '''

        robjects.r(rcode)

        return None


       


  
	