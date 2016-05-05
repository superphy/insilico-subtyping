#!/usr/bin/env Rscript

##########################################################
## File: pp_tests.R
##
##  Compare model-estimated posterior probability with 
##  empirical p-values obtained from simulated data.
##
##  The tree that is loaded will have its subtype labels shuffled
##  to create a larger test space.
##
##  Leave-One-Out cross-validation will be performed for each
##  label in shuffled tree and pp values and correct labels 
##  recorded in each iteration.
##
##
## author: Matt Whiteside
## copyright: Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative
## license: ASL
## version: 2.0
## maintainer: Matt Whiteside
## email: mdwhitesi@gmail.com
##
##########################################################


# Load libraries and source
source('phylotyper.R')
suppressMessages(phylotyper$loadInstallLibraries())



## MAIN 

# Process command-line args
# option_list = list(
# 	make_option(c("-o", "--out"), type="character", default="out.txt", 
#     	help="output file name [default: %default]", metavar="character_string")
# ); 
 
argparser = OptionParser(usage="%prog [options] tree_file subtype_file");
arguments = parse_args(argparser, positional_arguments=2);

args = arguments$args
treefile = args[1]
subtypefile = args[2]
output_dir = '../test/'

# Check arguments
if( file.access(treefile) == -1) {
	stop(sprintf("Specified tree file ( %s ) does not exist", treefile))
}

if( file.access(subtypefile) == -1) {
	stop(sprintf("Specified subtype file ( %s ) does not exist", treefile))
}


# Run tests

# Load tree
rs = loadSubtype(treefile,subtypefile)
tree = rs$tree; subtypes = rs$subtypes


# Iterate through esimtation procedures
estimation.methods = list(rerooting=1, simmap=4)
for(i in 1:length(estimation.methods)) {
	est.scheme = estimation.methods[[i]]
	est.name = names(estimation.methods)[i]
	
	# Shuffle tree subtypes and compute pp through Leave-One-Out CV
	for(j in 1:100) {


		pp = loocv(tree=tree, subtypes=subtypes, scheme=est.scheme))

		
	}

	# Summarize performance
	results = simulationSummary(subtypes, pp)

	# Plot posterior probability histogram
	file = 'posterior_probability_histogram_for_shuffled_tip_labels'
	png(filename=file.path(output_dir, paste(est.name, '_', validation, '_', file, '.png', sep='')))
	plotPPHistogram(results$test.results, subtypes)
	dev.off()

	print(est.name, "complete")

}











