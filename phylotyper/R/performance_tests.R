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
	
# Overlay posterior probabilities in tree plot
cat("Computing posterior probabilites for all internal nodes\n")

priorR = phylotyper$makePriors(tree, subtypes)
priorM = priorR$prior.matrix
fit <- rerootingMethod(tree, priorM, model='SYM', tips=FALSE)



source('phylotyper.R')
file = 'posterior_probability_tree_small'
#file = 'posterior_probability_tree_large'
dim = phylotyper$plotDim(tree)
r = 300
w = 3.38 * r
graphics.off()
# png(filename=file.path(output_dir, paste(file, '.png', sep='')),
#     width=dim[['x']],height=dim[['y']],res=dim[['res']]
# )
setEPS()
postscript(file.path(output_dir, paste(file, '.eps', sep='')),
	width=3.38,height=4,pointsize=6)
# tiff(filename=file.path(output_dir, paste(file, '.tiff', sep='')),
# 	width=w, height=w*1.25, res=r, pointsize=6, units='px', compression='zip')
phylotyper$plot.anc(tree,fit,subtypes)
graphics.off()



# Iterate through validation procedures
# Leave-One-Out CV
# 5-fold CV
est.scheme = 5 # only computes pp for tips
cvs = c('loocv')
for(validation in cvs) {
	cat("Running validation: ", validation, "\n")

	pp = do.call(validation, list(tree=tree, subtypes=subtypes, scheme=est.scheme))

	# Summarize performance
	results = simulationSummary(subtypes, pp, threshold=.85)

	# Write performance metrics to file
	file = 'performance_metrics'
	write.table(results$metrics, file=file.path(output_dir, paste(validation, '_', file, '.csv', sep='')),
		sep="\t",
		quote=FALSE)

	# Plots

	# Plot confusion matrix
	# file <- 'confusion_matrix'
	# fn = file.path(output_dir, paste(validation, '_', file, '.png', sep=''))
	# plotConfusionMatrix(results$confusion.matrix, fn)
	

	# Plot posterior probability histogram
	# file <- 'posterior_probability_histogram'
	# fn = file.path(output_dir, paste(validation, '_', file, '.png', sep=''))
	# plotPPHistogram(results$test.results, subtypes, fn)
	
	cat(validation, "complete", "\n")
}

	













