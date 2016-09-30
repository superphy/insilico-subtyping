#!/usr/bin/env Rscript

##########################################################
## File: denovo.R
##
##  Generate subtype assignment based on tree
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
source('subtype_properties.R')

# Process command-line args
option_list = list(
	make_option(c("-o", "--out"), type="character", default="../output/", 
    	help="output directory [default: %default]", metavar="character_string")
); 
 
argparser = OptionParser(usage="%prog [options] tree_file", option_list=option_list);
arguments = parse_args(argparser, positional_arguments=1);

args = arguments$args
opts = arguments$options
treefile = args[1]
output_dir = opts$out

# Load inputs
tree = loadSubtype(treefile)$tree

# Make new subtype

# Find distribution of branch lengths for genomes in same
# subtype
file <- 'tree_patristic_length_histogram'
fn = file.path(output_dir, paste(file, '.png', sep=''))
bldist = branchDistances(tree, plot.name=fn)
file <- 'tree_patristic_distribution_fit'
fn = file.path(output_dir, paste(file, '.png', sep=''))
fit = fitSubtypeDistribution(bldist, plot.name=fn)

# Assign tips based on subtype branch distribution
file <- 'denovo_subtype'
fn = file.path(output_dir, paste(file, '.png', sep=''))
subt = assignSubtypes(fit$same, tree, plot.name=fn)

file <- 'denovo_subtype'
fn = file.path(output_dir, paste(file, '.txt', sep=''))
write.table(data.frame(names(subt), as.character(subt)), file=fn, quote=FALSE, sep="\t", 
	row.names=FALSE, col.names=FALSE)

