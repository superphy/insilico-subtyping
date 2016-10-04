#!/usr/bin/env Rscript

##########################################################
## File: performance_testing.R
##
##  Evaluate methods/models for predicting subtype from
##  phylogenic trees
##
##  Needs to be run from phylotyper/R directory
## 
##
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

## MAIN

# Check work directory
wd = normalizePath(getwd())
if(!grepl("phylotyper/R/?$", wd, perl=TRUE)) stop(c("Incorrect work directory: ", wd))

# Load libraries and source
source('phylotyper.R')
suppressMessages(phylotyper$loadInstallLibraries())
source('test_functions.R')

# Process command-line args
option_list = list(
	make_option(c("-o", "--out"), type="character", default="../output/", 
    	help="output directory [default: %default]", metavar="character_string")
); 
 
argparser = OptionParser(usage="%prog [options] tree_file subtype_file", option_list=option_list);
arguments = parse_args(argparser, positional_arguments=2);

args = arguments$args
opts = arguments$options
treefile = args[1]
subtypefile = args[2]
output_dir = opts$out


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

# Run Mk model evaluation (this markov model is used in simmap and rerootingMethod)
# For large datasets, ARD and SYM computation takes a long time
# aic = evaluateModels(tree,subtypes)
# file = 'model_aic'
# write.table(aic, file=file.path(output_dir, paste(file, '.csv', sep='')),
# 	sep="\t",
# 	quote=FALSE)


	
# Overlay posterior probabilities in tree plot
cat("Computing posterior probabilites for all internal nodes: ", est.name, "\n")

priorR = phylotyper$makePriors(tree, subtypes)
priorM = priorR$prior.matrix
fit <- rerootingMethod(tree, priorM, model='ER', tips=FALSE)

file = 'posterior_probability_tree'
dim = phylotyper$plotDim(tree)
graphics.off()
png(filename=file.path(output_dir, paste(est.name, '_', file, '.png', sep='')),
    width=dim[['x']],height=dim[['y']],res=dim[['res']]
)
do.call(result$plot.function, list(tree=tree, fit=result$result, subtypes=subtypes))
graphics.off()

# Iterate through validation procedures
# Leave-One-Out CV
# 5-fold CV
est.scheme = 5 # only computes pp for tips
for(validation in c('loocv', 'kfcv')) {
	cat("Running validation: ", validation, "\n")

	pp = do.call(validation, c(tree=tree, subtypes=subtypes, scheme=est.scheme))

	# Summarize performance
	results = simulationSummary(subtypes, pp)

	# Write performance metrics to file
	file = 'performance_metrics'
	write.table(results$metrics, file=file.path(output_dir, paste(validation, '_', file, '.csv', sep='')),
		sep="\t",
		quote=FALSE)

	# Plots

	# Plot confusion matrix
	file <- 'confusion_matrix'
	fn = file.path(output_dir, paste(validation, '_', file, '.png', sep=''))
	plotConfusionMatrix(results$confusion.matrix, fn)
	

	# Plot posterior probability histogram
	file <- 'posterior_probability_histogram'
	fn = file.path(output_dir, paste(validation, '_', file, '.png', sep=''))
	plotPPHistogram(results$test.results, subtypes, fn)
	
	cat(validation, "complete", "\n")
}

	













