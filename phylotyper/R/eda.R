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
source('subtype_properties.R')

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

# Load tree
rs = loadSubtype(treefile,subtypefile)
tree = rs$tree; subtypes = rs$subtypes


# Subtype fixes
subtype.treatments <- list('original'=subtypes)
rs = reassign.subtypes(tree, subtypes)
subtype.treatments[['reassign']] = rs$new
names(subtype.treatments[['reassign']]) <- rownames(rs)
rs = monophyletic.subtypes(tree, subtypes)
subtype.treatments[['monophyletic']] = rs$new
names(subtype.treatments[['monophyletic']]) <- rownames(rs)


showdown <- list()

for(stname in names(subtype.treatments)) {
	subt <- subtype.treatments[[stname]]
	
	# Models
	models <- list('er'='ER')
	if(length(levels(subt)) < 10) {
		models <- c(models, 'sym'='SYM')
	}
	cat("Computing subtype-specific rate model parameters\n")
	rs = transition.rate.parameters(tree, subt)
	models[['custom']] = rs$model

	rs = transition.rate.parameters2(tree, subt)
	models[['custom2']] = rs$model

	# Running simulation
	for(modname in names(models)) {
		mod <- models[[modname]]
		est.scheme = 5 # only computes pp for tips
		pp = loocv(tree=tree, subtypes=subt, scheme=est.scheme, model=mod)

		# Summarize performance
		results = simulationSummary(subt, pp)

		# Save a copy of the fit
		fit <- fitMk(tree, subt, model=mod)

		showdown[[stname]][[modname]] = list(pp=pp, results=results, anc=fit)


		cat(modname, "complete", "\n")
	}

	cat(stname, "complete", "\n")
}

# Retrieve F-scores
performance <- data.frame(do.call('rbind', sapply(names(showdown), function(n) cbind(model=names(showdown[[n]]), treatment=n, fscore=0))), stringsAsFactors=FALSE)

for(stname in names(showdown)) {
	for(modname in names(showdown[[stname]])) {
		rs = showdown[[stname]][[modname]]$results$metrics

		fs = rs['ave / total','F-score']
		
		row = performance[,'treatment'] == stname & performance[,'model'] == modname
		performance[row,'fscore'] = fs
	}
}

performance$fscore = as.numeric(performance$fscore)

ggplot(performance, aes(x=model, y=fscore, color=treatment)) +
 	geom_point(stat='identity', size=3) + 
	ggtitle("F-score for different transistion model\nparameteritizations and subtype treatments")

	


	













