#!/usr/bin/env Rscript

##########################################################
## File: performance_tests_roc.R
##
##  Remove entire subtype from training data, count number of 
##  test set that are incorrectly classified as other subtype
##
##  Needs to be run from phylotyper/R directory
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
 
argparser = OptionParser(usage="%prog [options] tree_file subtype_file modelname", option_list=option_list);
arguments = parse_args(argparser, positional_arguments=3);

args = arguments$args
opts = arguments$options
treefile = args[1]
subtypefile = args[2]
modelname = args[3]
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

# Run cross-validation
res = kfcv(tree, subtypes, scheme=5, model=modelname)

# Compute FPR / TPR vs cutoff
pred = prediction(res$all$prediction, res$all$label)
fpr = performance(pred, 'fpr')
tpr = performance(pred, 'tpr')

plot(fpr, avg='vertical', lwd=1, col='red', spread.estimate="stderror", 
	ylab='Average rate across subtypes', xlab='Phylotyper posterior probability cutoff')
plot(tpr, avg='vertical', lwd=1, col='blue', spread.estimate="stderror", add=TRUE)
legend(0.6,0.6,c('FPR','TPR'),col=c('red','blue'),lwd=1)

# Compute PPV / TPR for max scheme
mx.pred = prediction(res$max$prediction, res$max$label)
plot(performance(mx.pred, 'ppv'), col='red', ylab='Precision', xlab='Phylotyper posterior probability cutoff')
plot(performance(mx.pred, 'tpr'), col='red', ylab='Recall', xlab='Phylotyper posterior probability cutoff')
plot(performance(mx.pred, 'f'), col='red', ylab='F1-Statistic', xlab='Phylotyper posterior probability cutoff')
	
# Compute the averate FPR when entire subtype is removed
fp = lsocv(tree, subtypes, scheme=5, model=modelname, threshold=0.85)

lso.fpr = mean(fp[,1]/fp[,2])






	













