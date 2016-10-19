#!/usr/bin/env Rscript

##########################################################
## File: reassign.R
##
##  Change subtype assignments that significantly deviate 
##  from main phylogenic distance distribution
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
 
argparser = OptionParser(usage="%prog [options] tree_file subtype_file revised_subtype_file", option_list=option_list);
arguments = parse_args(argparser, positional_arguments=1);

args = arguments$args
opts = arguments$options
treefile = args[1]
subtypefile = args[2]
newsubtypefile = args[3]
output_dir = opts$out

# Load inputs
rs = loadSubtype(treefile, subtypefile); tree <- rs$tree; subtypes <- rs$subtypes

# Compute branch lengths
file <- 'subtype_patristic_length_histogram'
fn = file.path(output_dir, paste(file, '.png', sep=''))
bdist = branchDistances(tree, subtypes, plot.name=fn)

# Isolate distances for genomes in same subtype
same = bdist[bdist$subtype == 'same',]
x = same$distance

# Outlier detection - remove the worst offenders
x <- x[!x %in% adjboxStats(x)$out]

# Try fitting the following models:
models <- c('exp','norm')
aics <- array(NA, dim=length(models))
names(aics) = models

for(m in models) {
	aics[m] = fitdist(x, distr=m)$aic
}

# Best fitting distribution
distr = names(aics)[which.min(aics)]
fit = fitdist(x, distr=distr)
file <- 'subtype_patristic_distribution_fit'
fn = file.path(output_dir, paste(file, '.png', sep=''))
graphics.off()
png(filename=fn)
plot(fit)
graphics.off()

# Observations that belong to this distribution
# Uses conservative cutoff - prefer to merge over splitting subtypes
pcutoff = 0.01
qmethod = paste('q',distr,sep='')
param = c(p=pcutoff, as.list(fit$estimate), lower.tail=FALSE)
dcutoff = do.call(qmethod, param)

# Find subtype groups that fall under this distance threshold
revised_same = bdist[bdist$distance < dcutoff,]
subt = assignSubtypes(revised_same, tree, plot.name=NULL)

# Make naming scheme
subsubtypes <- list()
for(s in levels(subt)) {
	grp <- names(subt)[subt == s]
	old_subtypes = unique(as.character(subtypes[grp]))
	
	grp_subtypes = array(NA, length(old_subtypes))
	names(grp_subtypes) = old_subtypes

	for(os in old_subtypes) {
		# Encountered old subtype second time
		# Indicates that it was split into different groups
		# Need to append sub-subtype designation
		
		if(os %in% names(subsubtypes)) {
			st = subsubtypes[[os]]

			if(st > 26) {
				stop("Run out of sub-subtype designations a-z")
			}

			grp_subtypes[os] = paste(os, letters[st], sep='')
			subsubtypes[[os]] = st + 1

		} else {
			# First time subtype was encountered
			# Re-use old subtype name
			grp_subtypes[os] = os
			subsubtypes[[os]] = 1

		}
	}

	# Collapse multiple subtypes (if any) into single name
	subtnm = paste(grp_subtypes, collapse='/')
	
	# Rename
	levels(subt)[levels(subt)==s] <- subtnm
}


# Save to file
genomes = names(subtypes)
ord = sort(genomes)
reassigned = data.frame(cbind(as.character(subt[ord]),as.character(subtypes[ord])), row.names=ord)
write.table(reassigned, file=newsubtypefile, quote=FALSE, row.names=TRUE, col.names=FALSE, sep='\t')


# Plot
file <- 'reassigned_subtype'
fn = file.path(output_dir, paste(file, '.png', sep=''))
dim = phylotyper$plotDim(tree)
graphics.off()
png(filename=file.path(output_dir, paste(file, '.png', sep='')),
    width=dim[['x']],height=dim[['y']],res=dim[['res']]
)
s1 <- reassigned[,1]
names(s1) <- rownames(reassigned)
s2 <- reassigned[,2]
names(s2) <- rownames(reassigned)
phylotyper$plot.subtype(tree, s1, s2)
graphics.off()







