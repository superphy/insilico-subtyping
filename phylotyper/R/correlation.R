#!/usr/bin/env Rscript

##########################################################
## File: correlation.R
##
##  Examine how the leave-one-out F-score changes when
##  the subtype labels in tree are shuffled.
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
x='dplyr'
if (!require(x,character.only = TRUE)) {
	install.packages(x,dep=TRUE)
	library(x)
	if(!require(x,character.only = TRUE)) stop(paste("Package not loaded: ",x))
}


# Process command-line args
option_list = list(
	make_option(c("-o", "--out"), type="character", default="../output/", 
    	help="output directory [default: %default]", metavar="character_string")
	make_option(c("-m", "--model"), type="character", default="ER", 
    	help="transition rate model [default: %default]", metavar="character_string")
); 
 
argparser = OptionParser(usage="%prog [options] tree_file subtype_file modelname", option_list=option_list);
arguments = parse_args(argparser, positional_arguments=3);

args = arguments$args
opts = arguments$options
treefile = args[1]
subtypefile = args[2]
modelname = opts$model
output_dir = opts$out


# Check arguments
if( file.access(treefile) == -1) {
	stop(sprintf("Specified tree file ( %s ) does not exist", treefile))
}
if( file.access(subtypefile) == -1) {
	stop(sprintf("Specified subtype file ( %s ) does not exist", treefile))
}
if(! modelname %in% c('ER','SYM','GROUP','ITER')) {
	stop(sprintf("Invalid model ( %s )", modelname))
}

# Load tree
rs = loadSubtype(treefile,subtypefile)
tree = rs$tree; subtypes = rs$subtypes

model <- NULL
if(modelname == 'GROUP') {
	rs = transition.rate.parameters(tree, subtypes)
	model = rs$model
} else if(modelname == 'ITER') {
	rs = transition.rate.parameters2(tree, subtypes)
	model = rs$model
} else {
	model = modelname
}

print(model)

shuffling <- function(subtdf, perc) {
	shuffledf = subtdf

	c = levels(subtypes)
	l = length(subtypes)
	n = floor(perc*l)
	
	sampled = sample(l, n)

	for(i in sampled) {
		orig = as.character(subtypes[i])
		new = sample(c[orig != c], 1)

		shuffledf[i] = new
	}

	return(shuffledf)
}

percents = seq(0,0.5,by=0.05)

iters = 30
nr = iters * 2 * (length(percents)-1) + 2 

rates = data.frame('metric'=character(nr), 'percent.shuffled'=numeric(nr), 'rate'=numeric(nr), stringsAsFactors = FALSE)

j = 1
for(p in percents) {
	last = iters
	if(p == 0) {
		last = 1
	}
	for(i in 1:last) {

		shuffled = shuffling(subtypes, p)

		# Run validation
		est.scheme = 5 # only computes pp for tips
		pp = loocv(tree=tree, subtypes=shuffled, scheme=5, model=model)

		# Summarize performance
		results = simulationSummary(subtypes, pp, threshold=0.85)

		rs = results$metrics
        fs = rs['ave / total','F-score']
        prec = rs['ave / total','precision']
        tpr = rs['ave / total','recall']

        # Save
        rates[j,1] <- 'fscore'
        rates[j,2] <- p
        rates[j,3] <- fs
        rates[j+1,1] <- 'recall'
        rates[j+1,2] <- p
        rates[j+1,3] <- tpr
       	
       	j <- j+2
	}
	
    print(paste("percent:",p,"complete", sep=' '))
}

# Summerize and plot
rate_summary <- rates %>%
	group_by(metric, percent.shuffled) %>%   # the grouping variable
    summarise(  # calculates the mean & SD of each group
    	mean.rate = mean(rate),  
    	sd.rate = sd(rate),
        n.rate = n()) %>%
    mutate(  # calculate the SE and CI for each group
    	se.rate = sd.rate / sqrt(n.rate),
        lower.ci.rate = mean.rate - qt(1 - (0.05 / 2), n.rate - 1) * se.rate,
        upper.ci.rate = mean.rate + qt(1 - (0.05 / 2), n.rate - 1) * se.rate )

# Plot 
pd <- position_dodge(0.005) # move them to the left and right
ggplot(rate_summary, aes(x=percent.shuffled, y=mean.rate, colour=metric, group=metric)) +
	geom_errorbar(aes(ymin=lower.ci.rate, ymax=upper.ci.rate), colour="black", width=.01, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=2, shape=21, fill="white") + # 21 is filled circle
    xlab("Percent of subtype labels randomized") +
    ylab("Performance metric value") +
    scale_colour_hue(name="Metric",    # Legend label, use darker colors
                     breaks=c("fscore", "recall"),
                     labels=c("F-score", "Recall"),
                     l=40) +                    # Use darker colors, lightness=40
    expand_limits(y=0) +                        # Expand y range
    scale_y_continuous(breaks=seq(0,1,0.1)) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(0,0),
          legend.position=c(0.01,0.01))

ggsave(file.path(output_dir, 'performance_for_randomized_labels.png')
