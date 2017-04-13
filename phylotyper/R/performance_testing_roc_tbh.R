#!/usr/bin/env Rscript

##########################################################
## File: performance_tests_roc_tbh.R
##
##  Summarize Top-BLAST hit approach performance results
##
##  R script uses outputs from:
##    1. python performance/cv.py kfcv ...
##    2. python performance/cv.py lco ...
##  
##  Needs to be run from phylotyper/R directory
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


# Process command-line args
option_list = list(
	make_option(c("-o", "--out"), type="character", default="../output/", 
    	help="output directory [default: %default]", metavar="character_string")
); 
 
argparser = OptionParser(usage="%prog [options] lco_results.tsv kfcv_results.tsv", option_list=option_list);
arguments = parse_args(argparser, positional_arguments=2);

args = arguments$args
opts = arguments$options
lcofile = args[1]
kfcvfile = args[2]
output_dir = opts$out


# Check arguments
if( file.access(treefile) == -1) {
	stop(sprintf("Specified tree file ( %s ) does not exist", treefile))
}

if( file.access(subtypefile) == -1) {
	stop(sprintf("Specified subtype file ( %s ) does not exist", treefile))
}

# Load data

# Read the max classifier results
fp = read.table(lcofile, sep='\t', as.is=TRUE, row.names=1, col.names=c('subtype','fp','n'), header=FALSE)

# Read % identity vs postives
con  <- file(kfcvfile, open = "r")
line <- readLines(con, n = 1, warn = FALSE)
mxpredictions <- as.numeric(strsplit(line, "\t")[[1]])
line <- readLines(con, n = 1, warn = FALSE)
mxlabels <- as.numeric(strsplit(line, "\t")[[1]])

results <- list(predictions=NULL, labels=NULL)

while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {

	subtype <- line
	line2 <- readLines(con, n = 1, warn = FALSE)
	predictions <- (strsplit(line2, "\t"))
    predictions <- c(as.numeric(predictions[[1]]), 0) # Add a superfluous TN to prevent NaN in FPR calculations
                                                      # i.e. something with 0% percent identity should get a negative assignment
    line3 <- readLines(con, n = 1, warn = FALSE)
	labels <- (strsplit(line3, "\t"))
    labels <- c(as.numeric(labels[[1]]),0)

    results$predictions[[subtype]] = predictions
    results$labels[[subtype]] = labels

} 

close(con)


# Compute FPR / TPR vs cutoff
pred = prediction(results$predictions, results$labels)
fpr = performance(pred, 'fpr')
tpr = performance(pred, 'tpr')

plot(fpr, avg='vertical', lwd=1, col='red', spread.estimate="stderror", 
	ylab='Average rate across subtypes', xlab='Phylotyper posterior probability cutoff')
plot(tpr, avg='vertical', lwd=1, col='blue', spread.estimate="stderror", add=TRUE)
legend(0.6,0.6,c('FPR','TPR'),col=c('red','blue'),lwd=1)

# Compute PPV / TPR for max scheme
mx.pred = prediction(mxpredictions, mxlabels)
plot(performance(mx.pred, 'ppv'), col='red', ylab='Precision', xlab='Percent Identify in BLAST hit')
plot(performance(mx.pred, 'tpr'), col='red', ylab='Recall', xlab='Percent Identify in BLAST hit')
plot(performance(mx.pred, 'f'), col='red', ylab='F1-Statistic', xlab='Percent Identify in BLAST hit')
	
# Compute the averate FPR when entire subtype is removed
lso.fpr = mean(fp[,1]/fp[,2])






	













