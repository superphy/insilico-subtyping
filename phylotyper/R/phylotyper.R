##########################################################
## File: phylotyper.R
##
##  Predict subtype state for a small number of unknowns in
##  the subtype DNA sequence tree
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


# Create enviroment to keep phylotyper functions
phylotyper = new.env()


phylotyper$runSubtypeProcedure = function(tree, priorM, scheme=1) {
	# A Wrapper that runs one of several possible subtype 
	# estimation procedures
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  priorM: a matrix of prior values indicating tip node subtype state
	#  scheme: an integer between 1:? indicating procedure to use.
	#    Possible procedures (phytool method, model):
	#      1=rerootingMethod, ER
	#      2=rerootingMethod, SYM
	#      3=rerootingMethod, ARD
	#      4=make.simmap, ER
	#      5=make.simmap, SYM
	#      6=make.simmap, ARD
	#
	# Returns list with:
	#   ?
	#

	if(scheme == 1) {
		# Rerooting method, ER model
		fit = rerootingMethod(tree, priorM, model='ER', tips=TRUE)
		return(list(result=fit, tip.pp=fit$marginal.anc[tree$tip.label,], plot.function='plotRR'))

	} else if(scheme == 2) {
		# Rerooting method, SYM model
		fit = rerootingMethod(tree, priorM, model='SYM', tips=TRUE)
		return(list(result=fit, tip.pp=fit$marginal.anc[tree$tip.label,], plot.function='plotRR'))

	} else if(scheme == 3) {
		# Rerooting method, SYM model
		fit = rerootingMethod(tree, priorM, model='SYM', tips=TRUE)
		return(list(result=fit, tip.pp=fit$marginal.anc[tree$tip.label,], plot.function='plotRR'))
	} else if(scheme == 4) {
		# simmap method, ER model
		trees = make.simmap(tree, priorM, model='ER', nsim=1000)
		fit = describe.simmap(trees, plot=FALSE)
		return(list(result=fit, tip.pp=fit$tips, plot.function='plotSM'))
	} 

}

phylotyper$makePriors = function(tree, subtypes) {
	# Generate prior matrix used as input for ACE
	# 
	# Tree tip nodes with no subtype in x, will be assigned
	# flat priors (equal likelihood of having any subtype)
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  subtypes: factor list of subtype assignments. List names
	#    must match tip names in tree, except for tips
	#    with unknown subtypes (which subtypes
	#    will be predicted for)
	#
	# Returns list with:
	#   prior.matrix: matrix containing prior values
	#   untyped: vector of tree tip names with no subtype 
	#
	
	# Find leaves with no subtype
	undefined = setdiff(tree$tip.label, names(subtypes))

	# Build prior matrix
	priors <- to.matrix(subtypes,seq=levels(subtypes))

	# Add flat priors for leaves with no subtype
	numstates = length(levels(subtypes))
	if(length(undefined) > 0)
		priors <- rbind(priors, matrix(1/numstates,
			length(undefined), numstates, dimnames=list(undefined)))


	return(list(prior.matrix=priors, untyped=undefined))
}

phylotyper$loadInstallLibraries = function(libloc="~/R/", repo="http://cran.stat.sfu.ca/") {
	# Attempt to load required libraries
	# If library not installed, attempt install
	#
	# Args:
	#  libloc: If defined, add path string to R's default internal Library paths variable
	#  repo: R repository to use in install.packages
	#
	# Returns:
	#   nothing
	#
	# Function will throw error if load/install fails
	# 

	# Location to install libraries to
	if(! is.null(libloc)) {
		.libPaths(c(libloc,.libPaths()))
	}
	
	# Set repo
	r = getOption("repos")
	r["CRAN"] = repo
	options(repos = r)

	# Install libraries from Bioconductor
	bioc.libs = c("Biostrings")
	if(!exists('biocLite')) {
		source("http://bioconductor.org/biocLite.R")
	}
	for(x in bioc.libs) {
		if (!require(x,character.only = TRUE)) {
	  		biocLite(x)
	    	if(!require(x,character.only = TRUE)) stop(paste("Package not loaded: ",x))
		}
	}
	
	# Install libraries from CRAN
	cran.libs = c("devtools", "ape", "phangorn", "RColorBrewer", "ggplot2")
	for(x in cran.libs) {
		if (!require(x,character.only = TRUE)) {
	  		install.packages(x,dep=TRUE)
	    	if(!require(x,character.only = TRUE)) stop(paste("Package not loaded: ",x))
		}
	}
	
	# Install libraries from Github
	git.libs = data.frame(lib=c("phytools"), git.path=c("liamrevell/phytools"))
	for(i in 1:nrow(git.libs)) {
    	lib = git.libs[i,1]
   		path = git.libs[i,2]
    	if (!require(lib, character.only = TRUE)) {
  			install_github(path)
    		if(!require(lib, character.only = TRUE)) stop(paste("Package not loaded:",lib))
		}
	}

	# Install older version of root, unroot, and is.root functions from APE packages
	# New version have bug
	source('root.R')
	source('rerootingMethod.R')
	
	TRUE
}

phylotyper$loadSubtype = function(treefile, stfile) {
	# Load tree and subtype assignments associated with subtype scheme
	#
	# Args:
	#  treefile: Path to tree in newick format
	#  stfile: Path to subtype assignments in format: tree_tip_label\tsubtype
	#
	# Returns:
	#   list:
	#     tree = phylo object
    #     untyped = list of tree_tip_labels not found in stfile
    #     subtypes = factor of subtype assignments. Names() match tree_tip_labels
	#
	# 
	# 

	tree = read.tree(treefile)
	st = read.table(stfile, sep="\t", row.names=1)

	# Convert to factor list
	subtypes = st[,1]
	names(subtypes) = rownames(st)

	# Make note of labels to find subtypes for
	undefined = setdiff(tree$tip.label, names(subtypes))

	return(list(tree=tree, untyped=undefined, subtypes=subtypes))
}

phylotyper$palette = function(subtypes) {
	# Generate a set of colors representing subtypes
	#
	# Args:
	#   subtypes: factor list of subtype assignments. List names
	#    must match tip names in tree, except for tips
	#    with unknown subtypes (which subtypes
	#    will be predicted for)
	#
	# Returns:
	#   list of color assignments. List names match subtype names
	# 
	# 

	colors = brewer.pal(12,"Set3")
	states = levels(subtypes)
	n = length(states)

	if(n > length(colors)) {
		stop("Number of subtypes exceeds available colors in palette. Please defined your own color palette.")
	}

	pal = colors[1:n]
	names(pal) = states

	return(pal)
}

phylotyper$plotRR = function(tree, fit, subtypes, plot.nodes=TRUE) {
	# plot tree with posterior probabilities from 
	# rerootingMethod function displayed as pies
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  fit: output object from rerootingMethod()
	#  subtypes: factor list of subtype assignments
	#  plot.nodes: display interior tree node ancestral state
	#    posterior probability
	#
	# Returns:
	#   nothing
	#

	cols = phylotyper$palette(subtypes)

	plot(tree)
	tiplabels(pie=fit$marginal.anc[tree$tip.label,], 
		piecol=cols,
		cex=0.3)
	if(plot.nodes) {
		nodelabels(pie=fit$marginal.anc[as.character(1:tree$Nnode+Ntip(tree)),],
			piecol=cols,
			cex=0.6)
	}
	add.simmap.legend(colors=cols,x=0.9*par()$usr[2],
		y=0.9*par()$usr[4],prompt=FALSE)
}

while("phylotyper" %in% search())
  detach("phylotyper")
attach(phylotyper)

