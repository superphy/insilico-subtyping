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

wd <- setwd(".")
setwd(wd)

# Create enviroment to keep phylotyper functions
phylotyper = new.env()


phylotyper$runSubtypeProcedure <- function(tree, priorM, scheme=1, ...) {
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
	#      5=tip.posterior.probability, ER
	#
	# Returns list with:
	#   ?
	#

	if(scheme == 1) {
		# Rerooting method, ER model
		capture.output(fit <- rerootingMethod(tree, priorM, model='ER', tips=TRUE))
		return(list(result=fit, tip.pp=fit$marginal.anc[tree$tip.label,], plot.function='plotRR'))

	} else if(scheme == 2) {
		# Rerooting method, SYM model
		fit = rerootingMethod(tree, priorM, model='SYM', tips=TRUE)
		return(list(result=fit, tip.pp=fit$marginal.anc[tree$tip.label,], plot.function='plotRR'))

	} else if(scheme == 3) {
		# Rerooting method, ARD model
		fit = rerootingMethod(tree, priorM, model='ARD', tips=TRUE)
		return(list(result=fit, tip.pp=fit$marginal.anc[tree$tip.label,], plot.function='plotRR'))

	} else if(scheme == 4) {
		# simmap method, ER model
		trees = make.simmap(tree, priorM, model='ER', nsim=1000)
		fit = describe.simmap(trees, plot=FALSE)
		return(list(result=fit, tip.pp=fit$tips, plot.function='plotSM'))

	} else if(scheme == 5) {
		# tip.posterior.probability method, ER model or specified using model= parameter
		if(hasArg(tips)) tips<-list(...)$tips
		else stop('Missing argument: tips')

		model='ER'
		if(hasArg(model)) model<-list(...)$model

		fixedQ=NULL
		if(hasArg(fixedQ)) fixedQ<-list(...)$fixedQ
		
		fit = tip.posterior.probability(tree, priorM, tips, model=model, fixedQ=fixedQ)
		return(list(result=fit, tip.pp=fit$marginal.anc, plot.function='plotTPP'))

	}

}

phylotyper$makePriors <- function(tree, subtypes) {
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

	# Reorder based on tip order
	priors <- priors[tree$tip.label,]

	return(list(prior.matrix=priors, untyped=undefined))
}

phylotyper$loadInstallLibraries <- function(libloc="~/R/", repo="http://cran.stat.sfu.ca/") {
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
	cran.libs = c("devtools", "ape", "phangorn", "RColorBrewer", "ggplot2", "optparse", "mclust",
		"robustbase", "fitdistrplus", "igraph")
	for(x in cran.libs) {
		if (!require(x,character.only = TRUE)) {
	  		install.packages(x,dep=TRUE)
	  		library(x)
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
	# Reload methods that rely on root.R methods
	source('rerootingMethod.R')
	source('utilities.R')
}

phylotyper$loadSubtype <- function(treefile, stfile=NULL, do.root=TRUE, resolve.polytomies=TRUE) {
	# Load tree and subtype assignments associated with subtype scheme
	#
	# Root tree at midpoint.
	#
	# Args:
	#  treefile: Path to tree in newick format
	#  stfile: Path to subtype assignments in format: tree_tip_label\tsubtype
	#  do.root: Boolean indicating if midpoint.root should be called on tree
	#  resolve.polytomies: Boolean indicating if multifurcating nodes should be expanded to bifurcating
	#
	# Returns:
	#   list:
	#     tree = phylo object
    #     untyped = list of tree_tip_labels not found in stfile
    #     subtypes = factor of subtype assignments. Names() match tree_tip_labels
	#
	# 
	# 

	otree = read.tree(treefile)
	tree = otree
	if(do.root) tree <- midpoint.root(tree)
	if(resolve.polytomies) {
		tree <- multi2di(tree,random=TRUE)
	}

	res=list(tree=tree, multif=otree)

	if(!is.null(stfile)) {
		# Load subtypes
		st = read.table(stfile, sep="\t", row.names=1)

		# Convert to factor list
		subtypes = st[,1]
		names(subtypes) = rownames(st)
		res[['subtypes']] = subtypes

		# Make note of labels to find subtypes for
		undefined = setdiff(tree$tip.label, names(subtypes))
		res[['untyped']] = undefined
	}
	
	return(res)
}

phylotyper$mypalette <- function(subtypes) {
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

	
	states = levels(subtypes)
	n = length(states)
	cols = rainbow(n, s=0.6,v=1)

	if(n > length(cols)) {
		stop("Number of subtypes exceeds available colors in palette. Please defined your own color palette.")
	}

	pal = cols[sample(1:n,n)]
	names(pal) = states

	pal = c(pal,'other'="#D3D3D3")

	return(pal)
}

phylotyper$piecolors <- function(marginals) {

	n = ncol(marginals)
	cutoff = 1/n
	reduced = cbind(marginals,other=0)

	reduced[,'other'] = apply(marginals, 1, function(x) { sum(x[which(x <= cutoff)]) })
	reduced[reduced <= cutoff] = 0

	return(reduced)
}

phylotyper$plot.subtype <- function(tree, subtypes, tip.subtypes=NULL) {
	# plot tree with subtype assignments overlayed
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  subtypes: factor list of subtype assignments
	#
	# Returns:
	#   nothing
	#

	cols = cols2 = phylotyper$mypalette(subtypes)
	subtypes2 = subtypes
	diff.tips = FALSE
	if(!is.null(tip.subtypes)){
		subtypes2 = tip.subtypes
		cols2 = phylotyper$mypalette(subtypes2)
		diff.tips = TRUE
	}
	pies = matrix(0, ncol=length(cols),nrow=length(subtypes))
	colnames(pies) = names(cols)
	rownames(pies) = tree$tip.label
	for(i in 1:nrow(pies)) {
		st = subtypes[rownames(pies)[i]]
		pies[i,st] = 1
	}

	plot(tree,label.offset=0.001,cex=0.7,type='fan',align.tip.label=TRUE,tip.col=cols2[subtypes2[tree$tip.label]])

	tiplabels(pie=pies, 
		piecol=cols,
		cex=0.2)
	
	add.simmap.legend(colors=cols,x=0.9*par()$usr[2],
		y=0.9*par()$usr[4],prompt=FALSE)

	# if(diff.tips) {
	# 	add.simmap.legend(colors=cols2,x=0.7*par()$usr[2],
	# 	y=0.9*par()$usr[4],prompt=FALSE)
	# }
}

phylotyper$plot.subtype2 <- function(tree, subtypes, tip.subtypes=NULL) {
	# plot tree with subtype assignments overlayed
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  subtypes: factor list of subtype assignments
	#
	# Returns:
	#   nothing
	#

	cols = cols2 = phylotyper$mypalette(subtypes)
	subtypes2 = subtypes
	diff.tips = FALSE
	if(!is.null(tip.subtypes)){
		subtypes2 = tip.subtypes
		cols2 = phylotyper$mypalette(subtypes2)
		diff.tips = TRUE
	}
	pies = matrix(0, ncol=length(cols),nrow=length(subtypes))
	colnames(pies) = names(cols)
	rownames(pies) = tree$tip.label
	for(i in 1:nrow(pies)) {
		st = subtypes[rownames(pies)[i]]
		pies[i,st] = 1
	}

	plot(tree,label.offset=0.001,cex=0.7,type='phylogram',align.tip.label=TRUE,tip.col=cols2[subtypes2[tree$tip.label]])

	tiplabels(pie=pies, 
		piecol=cols,
		cex=0.2)

	nodelabels(cex=0.2)
	
	add.simmap.legend(colors=cols,x=0.9*par()$usr[2],
		y=0.9*par()$usr[4],prompt=FALSE)

	# if(diff.tips) {
	# 	add.simmap.legend(colors=cols2,x=0.7*par()$usr[2],
	# 	y=0.9*par()$usr[4],prompt=FALSE)
	# }
}

phylotyper$plot.anc <- function(tree, fit, subtypes) {
	# plot tree with internal node marginal anc overlayed
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  subtypes: factor list of subtype assignments
	#
	# Returns:
	#   nothing
	#

	cols = phylotyper$mypalette(subtypes)

	plot(tree,label.offset=0.001,cex=0.8,type='fan',align.tip.label=TRUE,tip.col=cols[subtypes[tree$tip.label]],
		no.margin=TRUE)

	nodelabels(pie=phylotyper$piecolors(fit$marginal.anc),
		piecol=cols,
		cex=0.6)
	
	add.simmap.legend(colors=cols,x=0.01*par()$usr[2],
		y=0.95*par()$usr[4],prompt=FALSE,cex=0.8)
} 

phylotyper$plotRR <- function(tree, fit, subtypes) {
	# plot tree with posterior probabilities from 
	# rerootingMethod function displayed as pies
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  fit: output object from rerootingMethod()
	#  subtypes: factor list of subtype assignments
	#
	# Returns:
	#   nothing
	#

	cols = phylotyper$mypalette(subtypes)

	plot(tree,label.offset=0.001,cex=0.7,type='fan',align.tip.label=TRUE,tip.col=cols[subtypes[tree$tip.label]])

	tiplabels(pie=fit$marginal.anc[tree$tip.label,], 
		piecol=cols,
		cex=0.2)

	nodelabels(pie=phylotyper$piecolors(fit$marginal.anc[as.character(1:tree$Nnode+Ntip(tree)),]),
		piecol=cols,
		cex=0.4)
	
	add.simmap.legend(colors=cols,x=0.9*par()$usr[2],
		y=0.9*par()$usr[4],prompt=FALSE)
}

phylotyper$plotSM <- function(tree, fit, subtypes) {
	# plot tree with posterior probabilities from 
	# make.simmap function displayed as pies
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

	cols = phylotyper$mypalette(subtypes)

	plot(fit, colors=cols, fsize=0.6)

	add.simmap.legend(colors=cols,x=0.9*par()$usr[2],
		y=0.9*par()$usr[4],prompt=FALSE)
}

phylotyper$plotTPP <- function(fit, tree, subtypes) {
	# plot tree with tip posterior probability and conditional
	# likelihoods from tip.posterior.probability function
	#
	# Args:
	#  fit: output object from tip.posterior.probabilities()
	#  tree: not used, use fit$tree provided by the fit object
	#  subtypes: factor list of subtype assignments
	#
	# Returns:
	#   nothing
	#

	tree = fit$rerootedTree

	# Save original node assignments in the node.label field
	tree$node.label = as.character(1:tree$Nnode+Ntip(tree))

	# Reroot for better viewing
	tree = midpoint.root(tree)

	cols = phylotyper$mypalette(subtypes)
	probs = phylotyper$makePriors(tree, subtypes)$prior.matrix
	probs[rownames(fit$marginal.anc),] = fit$marginal.anc

	plot(tree,label.offset=0.001,cex=0.7,type='phylo',align.tip.label=TRUE,
		tip.col=cols[subtypes[tree$tip.label]], main=paste('Posterior probability and associated subtree conditional\nlikelihoods for tip',fit$tip))

	tiplabels(pie=probs[tree$tip.label,],
		piecol=cols,
		cex=0.2)

	original.nodes = tree$node.label[tree$node.label != "Root"]
	nstates = ncol(fit$conditional.likelihoods)
	root.prior = rep(1/nstates, nstates)
	probs2 = rbind(root.prior, fit$conditional.likelihoods[original.nodes,])
	nodelabels(pie=phylotyper$piecolors(probs2),
		piecol=cols,
		cex=0.4)
	
	add.simmap.legend(colors=cols,x=0.1*par()$usr[2],
		y=0.9*par()$usr[4],prompt=FALSE)

}


phylotyper$plotDim <- function(tree, type='fan') {
	# Compute X, Y dimensions for a plot to fit tree
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#
	# Returns:
	#   list with x, y, res names
	#

	x_scale = 3600
	y_scale = 16
	x_max = 1728
	y_max = 1728
	res = 72

	# X
	h = max(nodeHeights(tree))
	x = h * x_scale
	x = ifelse(x > x_max, x_max, x)

	if(type != 'fan') {
		# Y
		l = length(tree$tip.label)
		y = l * y_scale
		y = ifelse(y > y_max, y_max, y)
	} else {
		y <- x
	}

	return(list(x=x,y=y,res=res))
}

phylotyper$tip.posterior.probability <- function(tree,priorM,uncertain,model=c("ER","SYM"),fixedQ=NULL) {
	# Adapted from phytools rerootingMethod
	#
	# It computes the marginal posterior probabiliyt of a tip in tree by rooting at lowest point
	# on edge descending to tip. Unlike the original rerootingMethod that this function
	# is adapted from, it only iterates for each tip listed in the uncertain argument.
	#
	# Marginal ancestral state posterior probability is equivalent to the conditional 
	# scaled likelihoods at the root node of the tree (obtained by fitMk). The original 
	# rerootingMethod roots the tree at each internal node and tip edge to generate the 
	# marginals for each node tip in the tree.
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  priorM: a matrix of prior values indicating tip node subtype state
	#
	# Returns:
	#   list with x, y, res names
	#

	if(!inherits(tree,"phylo")) 
		stop("tree should be an object of class \"phylo\".")
	if(!is.matrix(model)) model<-model[1]
	n<-Ntip(tree)
	
	if(!is.matrix(priorM)){ 
		stop("priorM should be a matrix")
	}
	yy<-priorM
	yy<-yy[tree$tip.label,]
	yy<-yy/rowSums(yy)

	# Check the targets
	if(length(uncertain) < 1 || length(uncertain) > n) {
		stop('Missing / invalid argument: uncertain')
	}
	tl <- uncertain %in% tree$tip.label
	if(any(!tl)) {
		stop("Unrecognized tip label in uncertain argument")
	}
	nn = which(tree$tip.label %in% uncertain)
	if(any(nn > n)) {
		stop("Not a valid tip node in tree")
	}

	# Start with first target tip
	nn1 = nn[1]
	nn = nn[-1]

	# Compute marginal probability for first unknown tip and also transition matrix

	# Root tree at lowest point on edge to tip
	tt <- reroot(tree,nn1,tree$edge.length[which(tree$edge[,2]==nn1)])

	if(is.null(fixedQ)) {
		YY<-fitMk(tt,yy,model=model,output.liks=TRUE)
		Q <- phylotyper$makeQ(YY)
	} else {
		if(!is.matrix(fixedQ)){ 
			stop("fixedQ should be a matrix")
		}
		Q <- fixedQ
		YY<-fitMk(tt,yy,model=model,fixedQ=Q,output.liks=TRUE)
	}
	
	# Repeat for remaining tips
	ff<-function(nn){
		tt <- reroot(tree,nn,tree$edge.length[which(tree$edge[,2]==nn)])
		res = fitMk(tt,yy,model=model,fixedQ=Q,output.liks=TRUE)
		#print(paste('tip node:',nn))
		res$lik.anc[1,]
	}
	if(length(nn) > 0) {
		XX<-t(sapply(nn,ff))
		XX<-rbind(YY$lik.anc[1,],XX)
	}
	else {
		XX <- YY$lik.anc[1,,drop=FALSE]
	}
	rownames(XX)<-uncertain
	liks <- YY$lik.anc
	rownames(liks) <- 1:tt$Nnode+n

	return(list(loglik=YY$logLik,Q=Q,marginal.anc=XX,conditional.likelihoods=liks,rerootedTree=tt,tip=uncertain[1]))
}

phylotyper$makeQ <- function(YY) {
	# Compute rate matrix Q from a fitMk object
	Q<-matrix(c(0,YY$rates)[YY$index.matrix+1],length(YY$states),
			length(YY$states),dimnames=list(YY$states,YY$states))
		diag(Q)<--colSums(Q,na.rm=TRUE)

	return(Q)
}

while("phylotyper" %in% search())
  detach("phylotyper")
attach(phylotyper)

