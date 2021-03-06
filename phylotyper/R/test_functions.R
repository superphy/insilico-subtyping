#!/usr/bin/env Rscript

##########################################################
## File: test_functions.R
##
##  Functions to evaluate performance and models
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


loocv = function(tree, subtypes, scheme=5, model='ER') {
	# Estimate prediction accuracy by predicting subtypes 
	# for each tip using leave-one-out validation. 
	# Method sets tip prior to flat/unassigned and runs
	# esimtation procedure to get posterior probability for
	# unassigned tip
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  subtypes: factor list of subtype assignments. List names
	#    must match tip names in tree
	#
	# Returns:
	# 	data.frame of posterior probabilities 
	#

	n = length(subtypes)

	# Make prior matrix
	priorR = phylotyper$makePriors(tree, subtypes)
	priorM = priorR$prior.matrix
	nstates = ncol(priorM)

	# Make result matrix
	pp = data.frame(matrix(0,n,nstates))
	colnames(pp) = colnames(priorM)
	tips = rownames(priorM)
	pp = cbind("tip"=tips, "true_value"=as.character(subtypes[tips]), pp)
	nc = ncol(pp)

	# Flat prior
	flat = matrix(1/nstates,1,nstates)
	Q <- NULL

	for(i in 1:n) {

		tip = as.character(pp$tip[i])
		testprior = priorM
		testprior[i,] = flat

		testfit = phylotyper$runSubtypeProcedure(tree, testprior, scheme, tips=tip, model=model, fixedQ=Q)
		#print(testfit$result$Q)

		if(is.null(Q)) {
			Q<-testfit$result$Q
		}

		pp[i,3:nc] = testfit$tip.pp[tip, ]

		#cat("Completed test iteration",i,"...\n")
	}

	return(pp)
}

lsocv = function(tree, subtypes, scheme=5, model='ER', threshold=.7) {
	# Estimate false positive rates that occur when entire
	# subtype is removed from training set and used as input 
	# Method sets tip prior to flat/unassigned and runs
	# esimtation procedure to get posterior probability for
	# unassigned tip
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  subtypes: factor list of subtype assignments. List names
	#    must match tip names in tree
	#  scheme: Integer. See phylotyper$runSubtypeProcedure for list
	#  model: Transition matrix model ('ER','SYM', or parameter matrix)
	#  threshold: posterior probability cutoff for assigning subtype class
	#
	# Returns:
	# 	data.frame of posterior probabilities 
	#

	# Make prior matrix
	priorR = phylotyper$makePriors(tree, subtypes)
	priorM = priorR$prior.matrix
	states = colnames(priorM)
	nstates = length(states)

	# Make result matrix
	pp = data.frame(matrix(0,nstates,2))
	rownames(pp) = states
	colnames(pp) = c('fp','n')
	
	# Flat prior
	flat = matrix(1/(nstates-1),1,(nstates-1))
	Q <- NULL

	for(i in 1:nstates) {

		state = as.character(states[i])
		testprior = priorM[,-i]
		testers = which(apply(testprior, 1, sum) == 0)
		testprior[testers,] = flat
		tipnames = names(testers)
		
		stopifnot(all(tipnames %in% tree$tip.label))

		# Hide columns/rows if model is matrix
		if(!is.character(model)) {
			testmodel = model[-which(rownames(model) %in% state),]
			testmodel = testmodel[,-which(colnames(model) %in% state)]
		} else{
			testmodel = model
		}

		testfit = phylotyper$runSubtypeProcedure(tree, testprior, scheme, tips=tipnames, model=testmodel, fixedQ=Q)
		pp[state,'fp'] = sum(testfit$tip.pp[tipnames, ] >= threshold)
		pp[state, 'n'] = length(testers)

	}
	
	return(pp)
}


kfcv = function(tree, subtypes, scheme=5, model='ER') {
	# Estimate prediction accuracy by predicting subtypes 
	# for each tip using k-fold cross validation. 
	# Method sets tip priors to flat/unassigned and runs
	# esimtation procedure to get posterior probability for
	# unassigned tips.
	#
	# Args:
	#  tree: phylo object containing subtype tree
	#  subtypes: factor list of subtype assignments. List names
	#    must match tip names in tree
	#  scheme: Integer. See phylotyper$runSubtypeProcedure for list
	#  model: Transition matrix model ('ER','SYM', or parameter matrix)
	#
	# Returns:
	# 	 
	#

	# Compute k-fold size
	n = length(subtypes)
	iter = 100
	k = 5
	ksize = floor(n/k)
	cat("Cross-validation test set size: ",ksize,"\n")

	if(ksize < 5)
		stop("Test set too small for 10-fold cross validation. Test set size is less than 5.")

	# Make prior matrix
	priorR = phylotyper$makePriors(tree, subtypes)
	priorM = priorR$prior.matrix
	nstates = ncol(priorM)

	# Make result matrix
	subtype_classes = levels(subtypes)
	results = list()
	results$predictions = vector("list", length(subtype_classes))
	results$labels = vector("list", length(subtype_classes))
	names(results$predictions) <- c(subtype_classes)
	names(results$labels) <- c(subtype_classes)
	pickmax.results = list(predictions=NULL, labels=NULL)

	# Flat prior
	flat = matrix(1/nstates,1,nstates)
	
	for(i in 1:iter) {

		testset = sample(1:n, ksize)

		tips = names(subtypes)[testset]
		true_values = as.character(subtypes[testset])
		testprior = priorM
		testprior[tips,] = flat
		
		testfit = phylotyper$runSubtypeProcedure(tree, testprior, scheme, tips=tips, model=model)

		scores = testfit$tip.pp[tips, ]
		print(scores)
		
		for(j in 1:length(true_values)) {
			actual = true_values[j]

			for(c in colnames(scores)) {
				results$predictions[[c]] <- c(results$predictions[[c]], scores[j,c])
				results$labels[[c]] <- c(results$labels[[c]], if(c == actual) 1 else 0)
			}

			mx = which.max(scores[j,])
			c = colnames(scores)[mx]
			score = scores[j,mx]
			
			pickmax.results$predictions = c(pickmax.results$predictions, score)
			pickmax.results$labels = c(pickmax.results$labels, if(c == actual) 1 else 0)
		}

		cat("Completed test iteration",i,"...\n")
	}

	return(list(all=results,max=pickmax.results))
}


updateClassificationSums = function(true_class, predicted_class, count_matrix) {

	correct = true_class == predicted_class

	for(i in 1:nrow(count_matrix)) {
		c = rownames(count_matrix)[i]
		if(c == true_class) {
			# condition positive
			if(correct) {
				# tp
				count_matrix[i,'tp'] = count_matrix[i,'tp'] + 1
			} else {
				# fn
				count_matrix[i,'fn'] = count_matrix[i,'fn'] + 1
			}

		} else {
			# condition negative
			if(correct) {
				# tn
				count_matrix[i,'tn'] = count_matrix[i,'tn'] + 1
			} else {
				if(c == predicted_class) {
					# fp
					count_matrix[i,'fp'] = count_matrix[i,'fp'] + 1
				} else {
					# tn
					count_matrix[i,'tn'] = count_matrix[i,'tn'] + 1
				}
			}
		}
	}

	return(count_matrix)
}


simulationSummary = function(subtypes, pp, threshold=.85) {
	# Summarize performance of validation
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
	# 'over95': highest posterior probability, when pp is > 95
	# 

	subtype.states = levels(subtypes)
	if(!all(subtype.states %in% colnames(pp)))
		stop("Posterior probability data.frame is missing subtype states")

	# Record true positive, false positive etc for each subtype class
	classification_results = matrix(0, nrow=length(subtype.states), ncol=4)
	rownames(classification_results) = subtype.states
	colnames(classification_results) = c('tp','fp','tn','fn')
	metrics = matrix(0, nrow=length(subtype.states)+1, 4)
	rownames(metrics) = c(subtype.states,'ave / total')
	colnames(metrics) = c('precision','recall','F-score','support')

	# Process result for each test sample case using no cutoff, and cutoffs of 50, 70, 90 & 95
	cases = cbind(pp, 'max'=factor(c('unassigned'),levels=c('unassigned','correct','wrong')),
		'over50'=factor(c('unassigned'),levels=c('unassigned','correct','wrong')),
		'over70'=factor(c('unassigned'),levels=c('unassigned','correct','wrong')),
		'over90'=factor(c('unassigned'),levels=c('unassigned','correct','wrong')),
		'over95'=factor(c('unassigned'),levels=c('unassigned','correct','wrong'))
	)
	ppvals = cases[,subtype.states]

	# Confusion matrix
	confusion = matrix(0, nrow=length(subtype.states), ncol=length(subtype.states))
	rownames(confusion) = paste('true_',subtype.states,sep="")
	colnames(confusion) = paste('predicted_',subtype.states,sep="")

	for(i in 1:nrow(pp)) {
		true_st = as.character(pp$true_value[i])
		mxi = which.max(ppvals[i,])
		mx = ppvals[i,mxi]
		pred_st = subtype.states[mxi]

		# Increment confusion matrix counts
		k = which(subtype.states == true_st)
		confusion[k,mxi] = confusion[k,mxi] + 1
		
		# Correct
		gotit = FALSE
		if(true_st == pred_st) {
			gotit = TRUE
			cases[i, 'max'] = 'correct'
		} else {
			cases[i, 'max'] = 'wrong'
		}

		if(mx > 0.5) {
			if(gotit) {
				cases[i, 'over50'] = 'correct'
			} else {
				cases[i, 'over50'] = 'wrong'
			}

			if(mx > 0.7) {
				if(gotit) {
					cases[i, 'over70'] = 'correct'
				} else {
					cases[i, 'over70'] = 'wrong'
				}

				if(mx > 0.9) {
					if(gotit) {
						cases[i, 'over90'] = 'correct'
					} else {
						cases[i, 'over90'] = 'wrong'
					}

					if(mx > 0.95) {
						if(gotit) {
							cases[i, 'over95'] = 'correct'
						} else {
							cases[i, 'over95'] = 'wrong'
						}
					}
				}
			}
		}

		if(mx > threshold) {
			classification_results = updateClassificationSums(true_st,pred_st,classification_results)
		} else {
			# Missed result = false negative... and  a bunch of true negatives
			classification_results[true_st, 'fn'] = classification_results[true_st, 'fn'] + 1
			classification_results[which(rownames(classification_results) != true_st), 'tn'] = 
				classification_results[which(rownames(classification_results) != true_st), 'tn'] + 1
		}
	}

	# Relabel FN when n = 1.  There is no possibility of a true positive because the subtype is not
	# in the training set, so its not really a FN
	ncls = nrow(classification_results)
	testset = classification_results[1:ncls,'tp'] + classification_results[1:ncls,'fn']
	classification_results[testset == 1,'fn'] = 0

	# Compute metrics from classication results
	tots = apply(classification_results, 2, sum)
	classification_results = rbind(classification_results, "totals"=tots)
	ncls = nrow(classification_results)

	# precision
	metrics[1:ncls,'precision'] = classification_results[1:ncls,'tp'] / (classification_results[1:ncls,'tp']+classification_results[1:ncls,'fp'])
	# recall
	metrics[1:ncls,'recall'] = classification_results[1:ncls,'tp'] / (classification_results[1:ncls,'tp']+classification_results[1:ncls,'fn'])
	# F-score
	beta = 1
	denom = (beta^2*metrics[1:ncls,'precision'])+metrics[1:ncls,'recall']
	metrics[1:ncls,'F-score'] = (1+beta^2)*(metrics[1:ncls,'precision']*metrics[1:ncls,'recall']/denom)
	# support
	metrics[1:ncls,'support'] = classification_results[1:ncls,'tp'] + classification_results[1:ncls,'fn']

	
	# Check accuracy of posterior probabilities

	# # Record true state for each tip in each state column
	# match.colsubtype = t(sapply(pp$true_value, function(x) x == colnames(ppvals)))
	# colnames(match.colsubtype) = colnames(ppvals)

	# # Record frequencies for different pp windows
	# w.size = 0.1
	# boundaries = seq(from=w.size, to=(1-w.size), by=0.05)
	# pp.frequencies = matrix(0, length(boundaries), 4)
	# colnames(pp.frequencies) = c('accuracy', 'n', 'correct', 'pp')
	# pp.frequencies[,'pp'] = boundaries
	# rownames(pp.frequencies) = as.character(boundaries)

	# for(j in  boundaries) {
	# 	l = j-w.size
	# 	u = j+w.size

	# 	# Find pp that fall in window
	# 	for(k in 1:ncol(ppvals)) {
	# 		# Find rows with pp in that window
	# 		in.win = which(ppvals[,k] >= l & ppvals[,k] <= u)
	# 		n = length(in.win)

	# 		if(n > 0) {
	# 			ncorrect = sum(match.colsubtype[in.win,k])

	# 			i = which(pp.frequencies[,'pp'] == j)
	# 			pp.frequencies[i,'n'] = n + pp.frequencies[i,'n']
	# 			pp.frequencies[i,'correct'] = ncorrect + pp.frequencies[i,'correct']
	# 		}
	# 	}
	# }

	# pp.frequencies[,'accuracy'] = apply(pp.frequencies, 1, function(x) if(x['n'] == 0){ 0 } else { x['correct']/x['n'] } )


	return(list('metrics'=metrics, 'test.results'=cases, 'confusion.matrix'=confusion))
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

plotConfusionMatrix = function(counts, fn) {

	graphics.off()

	input.matrix.normalized <- sweep(counts, 1, rowSums(counts), FUN="/")

	colnames(input.matrix.normalized) = gsub("^predicted_","",colnames(results$confusion.matrix))
	rownames(input.matrix.normalized) = colnames(input.matrix.normalized)

	confusion <- as.data.frame(as.table(input.matrix.normalized))

	plot <- ggplot(confusion)
	plot + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_x_discrete(name="Actual Subtype") + scale_y_discrete(name="Predicted Subtype") + scale_fill_gradient(breaks=seq(from=0, to=1, by=.2)) + labs(fill="Normalized\nFrequency")
	ggsave(file=fn)

}

plotPPHistogram = function(test.results, subtypes, fn) {

	graphics.off()

	subtype.states = levels(subtypes)
	if(!all(subtype.states %in% colnames(test.results)))
		stop("Posterior probability result data.frame is missing subtype states")

	# posterior probability
	pp = test.results[,subtype.states]

	# The subset that can actually be predicted correct (i.e. has > 2 copies, in loocv)
	# kfcv is more complicated
	counts = table(test.results$true_value)
	states.subset = names(counts)[which(counts > 1)]

	# Create correct/incorrect distributions
	correct = unlist(sapply(1:length(subtype.states), function(i) { ss = subtype.states[i]; if(ss %in% states.subset) as.numeric(pp[test.results$true_value == ss,i]) } ))
	wrong = unlist(sapply(1:length(subtype.states), function(i) as.numeric(pp[test.results$true_value != subtype.states[i],i])))

	worst = unlist(sapply(1:length(subtype.states), function(i) max(as.numeric(pp[test.results$true_value != subtype.states[i],i]))))
	names(worst) = subtype.states

	dist = data.frame("result"=c(rep('positive',length(correct)), rep('negative',length(wrong))), 
		"posterior probability"=c(correct, wrong))

	# Plot with density and rug
	ggplot(dist, aes(x=posterior.probability, colour = result)) +
		geom_freqpoly(binwidth = .05, aes(y=..density..)) +
		geom_rug()

	ggsave(file=fn)
}












