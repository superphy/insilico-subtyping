#!/usr/bin/env Rscript

##########################################################
## File: performance_testing.R
##
##  Evaluate methods/models for predicting subtype from
##  phylogenic trees
##
## 
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


# Load libraries and source
source('phylotyper.R')
suppressMessages(phylotyper$loadInstallLibraries())
source('test_functions.R')


## MAIN 

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
aic = evaluateModels(tree,subtypes)
file = 'model_aic'
write.table(aic, file=file.path(output_dir, paste(file, '.csv', sep='')),
	sep="\t",
	quote=FALSE)

# Iterate through esimtation procedures
estimation.methods = list(rerooting=1, simmap=4)
for(i in 1:length(estimation.methods)) {
	est.scheme = estimation.methods[[i]]
	est.name = names(estimation.methods)[i]

	# Overlay posterior probabilities in tree plot
	print("Running estimation procedure: ", est.name)

	priorR = phylotyper$makePriors(tree, subtypes)
	priorM = priorR$prior.matrix
	result = phylotyper$runSubtypeProcedure(tree, priorM, est.scheme)
	file = 'posterior_probability_tree'
	png(filename=file.path(output_dir, paste(est.name, '_', file, '.png', sep='')))
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











