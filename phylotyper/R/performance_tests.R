##########################################################
## File: performance_testing.R
##
##  Evaluate methods/models for predicting subtype from
##  phylogenic trees
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

# Load libraries and source
source('phylotyper.R')
phylotyper$loadInstallLibraries()

# Functions
testTree = function() {

	# Simulate stochastic pure-birth tree
	tree<-pbtree(n=26,tip.label=LETTERS,scale=1)

	# Generate character transition matrix
	Q<-matrix(c(-1,0.8,0.2,0.8,-1.2,0.4,0.2,0.4,-0.6),3,3, dimnames=list(c("red","green","blue"),c("red","green","blue")))
	
	# Use transition matrix to simulate states on tree
	y<-as.factor(sim.history(tree,Q)$states)

	return(list(tree=tree, subtypes=y))
}

testSimulation = function(tree, subtypes, scheme=1) {
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

simulationSummary = function(subtypes, pp) {
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

	# Record true state for each tip in each state column
	match.colsubtype = t(sapply(names(subtypes), function(x) subtypes[x] == colnames(pp)))
	colnames(match.colsubtype) = colnames(pp)

	# Record frequencies for different pp windows
	w.size = 0.1 
	boundaries = seq(from=w.size, to=(1-w.size), by=0.05)
	pp.frequencies = matrix(0, length(boundaries), 4)
	colnames(pp.frequencies) = c('frequency', 'n', 'correct', 'pp')
	pp.frequencies[,'pp'] = boundaries
	rownames(pp.frequencies) = as.character(boundaries)

	
	for(j in  boundaries) {
		l = j-w.size
		u = j+w.size

		# Find pp that fall in window
		for(k in 1:ncol(pp)) {
			in.win = which(pp[,k] >= l & pp[,k] <= u)

			n = length(in.win)
			if(n > 0) {
				ncorrect = sum(match.colsubtype[in.win,k])

				i = which(pp.frequencies[,'pp'] == j)
				pp.frequencies[i,'n'] = n
				pp.frequencies[i,'correct'] = ncorrect
				pp.frequencies[i, 'frequency'] = ncorrect / n
			}
		}
	}


	return('frequency'=numCorrect, 'pp.evalulation'=pp.frequencies)
}

evaluateModels = function(tree, subtypes, models=list(equal="ER",
	symmetrical="SYM", different="ARD")) {
	# Compute AIC values for models used for state transition 
	# matrix Q
	#
	# Models tested:
	# 1. ER = equal rates
	# 2. SYM = symmetrical rates
	# 3. ARD = All rates different
	# 4. PGN = A single subtype gives rise to all other types
	# 
	# Args:
	#  tree: phylo object containing subtype tree
	#  subtypes: factor list of subtype assignments. List names
	#    must match tip names in tree
	#
	# Returns matrix containing AIC values for each model
	#

	nmodels = length(models)
	aic = matrix(nrow=nmodels)
	rownames(aic) = names(models)
	colnames(aic) = 'AIC'

	for(i in 1:nmodels) {
		m = models[[i]]
		n = names(models)[i]
		fit = fitMk(tree,subtypes,model=m)
		aic[i,1] = AIC(fit)
	}

	return(aic)
}


# Run tests

# Get tree
rs = testTree()
tree = rs$tree; subtypes = rs$subtypes

# Run model evaluation
aic = evaluateModels(tree,subtypes)

# Run one method/model parameter set
#pp = testSimulation(tree, subtypes)

# Summarize performance
#results = simulationSummary(subtypes, pp)

