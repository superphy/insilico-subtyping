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
        self.cwd = os.getcwd()
        self.rlib = config.get('R', 'lib')
        self.rrepo = config.get('R', 'repo')
        self.rscript = config.get('R', 'rscript')
        self.rfiles = {
            'phylotyper': os.path.join(self.rwd, 'phylotyper.R'),
            'tests': os.path.join(self.rwd, 'test_functions.R'),
            'performance_tests': os.path.join(self.rwd, 'performance_tests.R'),
            'subtype_properties': os.path.join(self.rwd, 'subtype_properties.R')
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
        robjects.r('setwd("%s")' % self.cwd)
        return assignments
  		

    def run(self, tree, subtypes, output_dir):
        """Run phylotyper R script
        
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


    def evaluate(self, tree, subtypes, rate_matrix, output_dir, accession_map):
        """Evaluate new subtype scheme
        
        Runs performance tests and reports accuracy. Identifies anomylous subtypes.
        Builds transition matrix to use in subtyping.
        
        Args:
            tree (str): Filepath to newick tree containing both typed & untyped genes
            subtypes (str): Filepath to tab-delimited subtype assignments for typed genes
            rate_matrix (str): Filepath for R serialized file containing transition rate matrix from winning model
            output_dir (str): Filepath to output directory
            accession_map (dict): Value from SeqDict.accession_map()
            

        Returns:
            results are output to output_dir
            
        """

        # Set work dir
        robjects.r('setwd("{}")'.format(self.rwd))

        # Load Phylotyper R functions
        rcode = 'source("{}")'.format(self.rfiles['phylotyper']) 
        robjects.r(rcode)

        # Load libraries
        rcode = 'suppressMessages(phylotyper$loadInstallLibraries(libloc="{}",repo="{}"))'.format(self.rlib, self.rrepo) 
        robjects.r(rcode)

        # Load test functions
        rcode = 'source("{}"); source("{}");'.format(self.rfiles['tests'], self.rfiles['subtype_properties']) 
        robjects.r(rcode)

        # Load data files
        rcode = 'rs = phylotyper$loadSubtype("{}","{}"); tree = rs$tree; subtypes = rs$subtypes'.format(tree, subtypes)
        robjects.r(rcode)

        # Check for weird subtypes
        rcode = '''
        fixes = reassign.subtypes(tree, subtypes);
        row.names(fixes)[as.character(fixes$new) != as.character(fixes$prev)]
        '''
        suspect = robjects.r(rcode)

        if suspect:
            newvalues = robjects.r('as.character(fixes$new)')
            genomes = robjects.r('row.names(fixes)')
           
            # Save to file
            newsubtype_file = os.path.join(output_dir, 'phylotyper_proposed_subtypes.csv')
            with open(newsubtype_file, 'w') as outfh:
                for i in xrange(len(newvalues)):
                    original_genome_labels = accession_map[str(genomes[i])]
                    for g in original_genome_labels:
                        outfh.write('{}\t{}\n'.format(g, newvalues[i]))

            self.logger.warn(
                '''\nSuspicious Subtypes! Detected anomylous subtype assignments for the following sequences
                that do not follow phylogenetic groupings:
                  {}
                A proposed subtype assignment has been generated and is available in:
                  {}
                Review the subtypes and if you want to use this subtype assignment, replace your subtype file 
                with this one and re-run the build script.
                '''.format('\n  '.join(tuple(suspect)), newsubtype_file, newsubtype_file))


        # Examine different transition matrix models
        rcode = '''
        models <- list('equal_rates'='ER')
        if(length(levels(subtypes)) < 10) {
            models <- c(models, 'symmetric_rates'='SYM')
        }
        
        rs = transition.rate.parameters(tree, subtypes)
        models[['iterative_binning']] = rs$model

        rs = transition.rate.parameters2(tree, subtypes)
        models[['small_distance_binning']] = rs$model
        '''
        robjects.r(rcode)

        # Assess performance of models via leave-one-out cross-validation
        rcode = '''
        showdown <- list()
        for(modname in names(models)) {
            mod <- models[[modname]]
            est.scheme = 5 # only computes pp for tips
            pp = loocv(tree=tree, subtypes=subtypes, scheme=est.scheme, model=mod)

            # Summarize performance
            results = simulationSummary(subtypes, pp)

            # Save a copy of the fit
            fit <- fitMk(tree, subtypes, model=mod)

            showdown[[modname]] = list(pp=pp, results=results, anc=fit)
        }
        '''
        robjects.r(rcode)

        # Retrieve F-scores
        rcode = '''
        performance <- data.frame(model=names(models), fscore=0, stringsAsFactors=FALSE)

        for(modname in names(showdown)) {
            rs = showdown[[modname]]$results$metrics
            fs = rs['ave / total','F-score']
        
            row = performance[,'model'] == modname
            performance[row,'fscore'] = fs
        }

        # find best model, with the following tie-break precedence: ER, SYM, Custom, Custom2
        modrow = min(which.max(performance$fscore))
        bestmodel = performance$model[modrow]
        bestfscore = performance$fscore[modrow]
        '''
        robjects.r(rcode)
        bestmodel = tuple(robjects.r('bestmodel'))[0]
        bestfscore = tuple(robjects.r('bestfscore'))[0]

        self.logger.info("The rate matrix model with highest accuracy is {}. Its associated F1-score: {}".format(bestmodel, bestfscore))

        if bestfscore < 0.9:
            # Not an ideal phylogenetic tree for predicting subtypes
            self.logger.warn(
                '''\nLow accuracy! The calculated F1-score from a leave-one-out cross-validation is {},
                suggesting that the ability of this phylogenetic tree to predict subtypes is questionable (highest scoring rate matrix model: {}).
                '''.format(bestfscore, bestmodel))

        # Save rate matrix
        rcode = 'Q = phylotyper$makeQ(showdown[[bestmodel]]$anc); saveRDS(Q, "{}")'.format(rate_matrix)
        robjects.r(rcode)

        # Return to the main directory
        robjects.r('setwd("%s")' % self.cwd)




        


       


  
	