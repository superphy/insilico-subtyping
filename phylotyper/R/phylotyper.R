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

phylotyper$plotRR = function(tree, result, plot.nodes=TRUE) {
	# plot tree with posterior probabilities from 
	# rerootingMethod function displayed as pies
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  result: output object from rerootingMethod()
	#  plot.nodes: display interior tree node ancestral state
	#    posterior probability
	#
	# Returns:
	#   nothing
	#

	plot(tree)
	tiplabels(pie=result$marginal.anc[tree$tip.label,], cex=0.5)
	if(plot.nodes) {
		nodelabels(pie=result$marginal.anc[as.character(1:tree$Nnode+Ntip(tree)),],cex=0.6)
	}
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
	source("http://bioconductor.org/biocLite.R")
	for(x in bioc.libs) {
		if (!require(x,character.only = TRUE)) {
	  		biocLite(x)
	    	if(!require(x,character.only = TRUE)) stop(paste("Package not loaded: ",x))
		}
	}
	
	# Install libraries from CRAN
	cran.libs = c("devtools", "ape", "phangorn")
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
	
	TRUE
}



while("phylotyper" %in% search())
  detach("phylotyper")
attach(phylotyper)

