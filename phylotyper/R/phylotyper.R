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
	undefined = setdiff(tree$tip.label, names(x))

	# Build prior matrix
	priors <- to.matrix(x,seq=levels(x))

	# Add flat priors for leaves with no subtype
	numstates = length(levels(x))
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

phylotyper$testTree = function() {

	# Simulate stochastic pure-birth tree
	tree<-pbtree(n=26,tip.label=LETTERS,scale=1)

	# Generate character transition matrix
	Q<-matrix(c(-1,0.8,0.2,0.8,-1.2,0.4,0.2,0.4,-0.6),3,3, dimnames=list(c("red","green","blue"),c("red","green","blue")))
	
	# Use transition matrix to simulate states on tree
	y<-as.factor(sim.history(tree,Q)$states)

	# Set some states as unknown
	y.sampled<-sample(y,20)
	y.sampled


	return(list(tree=tree, subtypes=y.sampled))
}

phylotyper$testSimulation = function(tree, subtypes, scheme=1) {
	# Estimate tip subtype posterior probability for each 
	# by setting tip prior to flat/unassigned and running
	# esimtation procedure
	# 
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  subtypes: factor list of subtype assignments. List names
	#    must match tip names in tree
	#
	# Returns list with:
	#   prior.matrix: matrix containing prior values
	#   untyped: vector of tree tip names with no subtype 
	#

	n = length(subtypes)

	# Make prior matrix
	priorR = phylotyper$makePriors(tree, subtypes)
	priorM = priorR$prior.matrix
	nstates = ncol(priorM)

	# Make result matrix
	pp = matrix(0,n,nstates)
	rownames(pp) = names(subtypes)
	colnames(pp) = colnames(priorM)

	# Flat prior
	flat = matrix(1/nstates,1,nstates)

	for(i in 1:n) {

		tip = names(subtypes)[i]
		testprior = priorM
		testprior[i,] = flat

		testfit = phylotyper$runSubtypeProcedure(tree, testprior, scheme)

		pp[i, ] = testfit$tip.pp[tip, ]

		cat("Completed test iteration",i,"...\n")
	}

	return(pp)
}

phylotyper$simulationSummary = function(subtypes, pp) {
	# Summarize performance of simulation
	# 
	# 
	#
	# Args:
	#  subtypes: factor list of subtype assignments. List names
	#    must match tip names in tree
	#  pp: posterior probability for each subtype obtained for 
	#    tip when true subtype masked and prediction procedure was run
	#
	# Returns list with:
	#   prior.matrix: matrix containing prior values
	#   untyped: vector of tree tip names with no subtype 
	#

	# Record frequency of correctly predicted tip subtype
	# when selecting subtype with:
	# 'max': highest posterior probability
	# 'over50': highest posterior probability, when pp is > 50
	# 'over70': highest posterior probability, when pp is > 70
	# 'over90': highest posterior probability, when pp is > 90

	
	numCorrect = matrix(0, 3, 4)
	rownames(numCorrect) = c('frequency', 'n', 'correct')
	colnames(numCorrect) = c('max', 'over50', 'over70', 'over90')

	subtype.states = colnames(pp)

	for(i in names(subtypes)) {
		st = subtypes[i]
		mxi = which.max(pp[i,])
		mx = pp[i,mxi]
		
		# Total
		numCorrect['n', 'max'] = numCorrect['n', 'max'] + 1

		# Correct
		gotit = FALSE
		if(st == subtype.states[mxi]) {
			numCorrect['correct', 'max'] = numCorrect['correct', 'max'] + 1
			gotit = TRUE
		}

		if(mx > 0.5) {
			numCorrect['n', 'over50'] = numCorrect['n', 'over50'] + 1
			if(gotit)
				numCorrect['correct', 'over50'] = numCorrect['correct', 'over50'] + 1

			if(mx > 0.7) {
				numCorrect['n', 'over70'] = numCorrect['n', 'over70'] + 1
				if(gotit)
					numCorrect['correct', 'over70'] = numCorrect['correct', 'over70'] + 1

				if(mx > 0.9) {
					numCorrect['n', 'over90'] = numCorrect['n', 'over90'] + 1
					if(gotit)
						numCorrect['correct', 'over90'] = numCorrect['correct', 'over90'] + 1

				}
			}
		}

		# Frequencies
		numCorrect['frequency',] = apply(numCorrect, 2, function(x) if(x['n'] == 0){ 0 } else { x['correct']/x['n'] } )

	}

	# Check accuracy of posterior probabilities
	w.size = 2 
	for(j in seq(from=w.size, to=(100-w.size), by=2) ) {
		l = j-w.size
		u = j+w.size

		apply(
	}


	return(numCorrect)
}


while("phylotyper" %in% search())
  detach("phylotyper")
attach(phylotyper)

