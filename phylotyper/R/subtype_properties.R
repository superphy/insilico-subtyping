#!/usr/bin/env Rscript

##########################################################
## File: subtype_properties.R
##
##  Functions to 
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
branchDistances <- function(tree, subtypes=NULL, plot.name=NULL) {
# Compute patristic distance between tips in same vs different subtype groups
#
# Args:
#  tree: phylo object containing subtype tree
#  subtypes[OPTIONAL]: factor list of subtype assignments. If provided, row with 
#      same/different classifcation will be provided
#  plot.name[OPTIONAL]: Plot output filename
#
# Returns:
#     data.frame of same/different subtype classification and corresponding distance 
#

	do.subtypes = !is.null(subtypes)

	d = cophenetic(tree)
	tips = tree$tip.label
	pairs = combn(tips,2)

	nr = 3
	rn = c('t1','t2','distance')
	if(do.subtypes) {
		nr <- 4
		rn <- c(rn, 'subtype')
	}
	df = matrix(ncol=ncol(pairs),nrow=nr)
	rownames(df) = rn
	

	for(c in 1:ncol(pairs)) {
		pair = pairs[,c]

		a = pair[1]
		b = pair[2]

		df[1:3,c] = c(a,b,d[a,b])

		if(do.subtypes) {
			if(subtypes[a] == subtypes[b]) {
				df[4,c] = 'same'
			}
			else {
				df[4,c] = 'diff'
			}
		}
	}

	# Plot with density and rug
	df = data.frame(t(df), stringsAsFactors=FALSE)
	df$distance = as.numeric(df$distance)
	if(!is.null(plot.name)) {
		if(do.subtypes) {
			ggplot(df, aes(x = distance, colour = subtype)) +
				geom_freqpoly(binwidth = .01)
		} else {
			ggplot(df, aes(x = distance)) +
				geom_histogram(binwidth = .01)
		}
		
		ggsave(file=plot.name)
	}

	return(df)
}


fitSubtypeDistribution <- function(sdist, plot.name=NULL, ...) {
# Identify the branch length distribution that correspond to 
# tips in the same subtype
#
# Args:
#  sdist: data.frame with branch distance (named $distance)
#  plot.name[OPTIONAL]: Plot output filename
#
# Returns:
#   list:
#     fit = Mclust object
#     same = subset of sdist for tips that are in same subtype
#

	fit = Mclust(sdist$distance, ...)

	if(!is.null(plot.name)) {
		plotMclustFit(fit, plot.name)
	}

	# Find the smallest distribution
	comp = as.character(names(which.min(fit$parameters$mean)))

	# Check the properties of the distribution

	# Retrieve subset of tip pairs in same subtype
	same = sdist[fit$classification %in% comp,]

	return(list(fit=fit, same=same))
}

plotMclustFit <- function(fit, plot.name) {
# Identify the branch length distribution that correspond to 
# tips in the same subtype
#
# Args:
#  fit: Mclust object
#  plot.name: Plot output filename
#
# Returns:
#   Nothing
#

	if(!is.null(plot.name)) {
		graphics.off()
		png(filename=plot.name,
	        width=12*72,height=12*72
	    )
		layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
		plot(fit, what="density")
		plot(fit, what="BIC")
		plot(fit, what="uncertainty")
		plot(fit, what="classification")
		graphics.off()
	}
}

mykeyFunc = function(trow) {
	# Internal function

	return(paste(as.character(sort(trow)),collapse='__'))
}


assignSubtrees <- function(pdist, tree, patristic.distance, plot.name=NULL, verbose=FALSE) {
# Identify subtrees where all patristic distances between leaves are below 
# specified threshold
#
# Args:
#  pdist: data.frame with tip.labels pairs and patristic distances
#  tree: phylo object containing subtype tree
#  patristic.distance: patristic distance threshold
#  plot.name[OPTIONAL]: Plot output filename
#
# Returns:
#   factor of tip.label names assigned to subtype
#

	ntips<-Ntip(tree)
	nnodes <- Nnode(tree)

	# Subtype group assignment
	cnum <- 0 # cluster number
	assign <- rep(0,ntips)
	names(assign) <- tree$tip.label

	# Lookup key for patristic data.frame
	row.keys <- as.character(apply(pdist[,c('t1','t2')], 1, mykeyFunc))

	# Subtree leaves 
	subtree.tips <- list()

	# Record nodes that are break points
	breaks <- c(rep(FALSE, ntips), rep(FALSE, nnodes))

	# DFS tree
	igraph.tree <- graph.edgelist(tree$edge) # tree in igraph form
	dfs <- graph.dfs(igraph.tree,root=ntips+1,neimode='out', order=TRUE,dist=TRUE)

	# Travese the tree in reverse depth first order
	ord <- rev(dfs$order)
	for(i in 1:length(ord)) {
		node <- ord[i]
    	
    	if(node <= ntips) {
    		# Skip leaves
    		next
    	}
    	
    	# Save new leaves encountered in this internal node
    	this.children <- Children(tree, node)
    	is.leaf <- sapply(this.children, function(x) x <= ntips)
    	inodes <- this.children[!is.leaf]

    	leaves <- list()
    	if(length(inodes) > 0) {
    		if(!all(as.character(inodes) %in% names(subtree.tips))) {
    			stop('Error in dfs traversal of tree')
    		}
    		leaves <- lapply(inodes, function(j) { subtree.tips[[as.character(j)]] } )
    	}
    	leaves <- as.list(c(unlist(leaves), this.children[is.leaf]))
    	if(verbose)
    		print(paste("At internal node",node,"with tips:",paste(tree$tip.label[unlist(leaves)], collapse=',')))

    	subtree.tips[[as.character(node)]] <- leaves

    	# If one of the descendants of this node was a break point
    	# Assign all subtrees to separate subtypes if they haven't been assigned
    	broken = any(sapply(inodes, function(x) breaks[x]))

    	if(broken) {
    		# Descendant was break point

    		if(verbose) {
	    		print(paste("Descendant in node",node,"is a break point"))
	    	}

    		# Assign any unassigned subtrees to subtypes
    		for(c in this.children) {
    			if(c <= ntips) {
    				# Leaf
    				cnum = cnum + 1
    				assign[c] = cnum
    			}
    			else {
    				# Subtree
    				if(!breaks[c]) {
    					# Subtree was not assiged
    					
	    				descendants = Descendants(tree, c, 'tips')[[1]]
	    				cnum = cnum + 1
	    				assign[descendants] = cnum
    				}
    			}
    		}

    		breaks[node] = TRUE

    	} else {
    		# Check if this node is a break point

    		# Compute inter patristic distances between tips in subtree
	    	# Assumes lower subtrees have already passed condition, so
	    	# only checks new tip pairs added by this subtree
	    	r = length(leaves)
	    	pairwise.dists <- unlist(sapply(1:r, function(x) if(x < r) sapply((x+1):r, function(y) {

				l1 = leaves[[x]]
				l2 = leaves[[y]]

				for(s in l1) {
					for(t in l2) {
						k = mykeyFunc(c(tree$tip.label[s], tree$tip.label[t]))
						return(pdist$distance[row.keys == k])
					}
				}
				
	    	})))

	    	if(verbose) {
	    		mx <- max(pairwise.dists)
	    		print(paste("Node",node,"largest pairwise patristic distance between tips:", mx))
	    	}
	    	
	    	if(all(pairwise.dists < patristic.distance)) {
	    		# This subtree passes condition
	    		if(verbose) {
		    		print(paste("Node",node,"is valid subtype cluster"))
		    	}
		    	next

	    	} else {
	    		# This node is the break point, subtrees must be separate subtypes

	    		if(verbose) {
		    		print(paste("Node",node,"is a break point"))
		    	}

	    		for(c in this.children) {
	    			if(c <= ntips) {
	    				# Leaf
	    				cnum = cnum + 1
	    				assign[c] = cnum
	    			}
	    			else {
	    				# Subtree
	    				descendants = Descendants(tree, c, 'tips')[[1]]
	    				cnum = cnum + 1
	    				assign[descendants] = cnum
	    			}
	    		}

	    		breaks[node] = TRUE

	    	}
    	}
    }

    subtype.subtrees = factor(assign)
	if(!is.null(plot.name)) {
		graphics.off()
		png(filename=plot.name,
	        width=12*72,height=12*72
	    )
		phylotyper$plot.subtype2(tree, subtype.subtrees)
		graphics.off()
	}

    	
    return(subtype.subtrees)

}

lineup <- function(pdist, tree, subtypes, Q=NULL) {
# Compare transition 
#
# Args:
#  pdist: data.frame with tip.labels pairs and patristic distances
#  tree: phylo object containing subtype tree
#  patristic.distance: patristic distance threshold
#  plot.name[OPTIONAL]: Plot output filename
#
# Returns:
#   factor of tip.label names assigned to subtype
#

	if(is.null(Q)) {
		fit <- fitMk(tree, subtypes, model='SYM', tips=FALSE)
		Q <- matrix(c(0,fit$rates)[fit$index.matrix+1],length(fit$states),
			length(fit$states),dimnames=list(fit$states,fit$states))
		diag(Q)<--colSums(Q,na.rm=TRUE)
	}

	xx <- mppd(pdist, subtypes)

	states <- sort(levels(subtypes))
	state.pairs <- combn(states,2)
	
	df <- data.frame(t(state.pairs), 
		apply(state.pairs, 2, function(x) Q[x[1],x[2]]),
		apply(state.pairs, 2, function(x) xx[x[1],x[2]]),
		apply(state.pairs, 2, function(x) mean(xx[x[1],x[1]], xx[x[2],x[2]]))
	)

	colnames(df) <- c('s1','s2','rate','mppd','inner')
	df$label = apply(df[,c('s1','s2')], 1, function(x) paste(x, collapse='_'))
	#label.dat <- df[df$rate != 0,]
	label.dat <- df

	ggplot(data=df, aes(x=mppd, y=rate, color=inner)) + 
    	geom_point() +
    	geom_text(data=label.dat,aes(x=mppd, y=rate, label=label),hjust=-.2) +
    	scale_color_gradient(low="blue", high="red")
}

mppd <- function(pdist, subtypes) {
	# Get median patristic distances

	states <- sort(levels(subtypes))
	n <- length(states)
	state.groups <- split(names(subtypes), subtypes)
	state.pairs <- combn(states,2)
	keyfunc <- function(t1,t2) {
		paste( sort(c(t1,t2)), collapse='__')
	} 
	pdist$key <- apply(pdist, 1, function(r) { keyfunc(r['t1'],r['t2']) })
	mppd = matrix(NA, nrow=n, ncol=n)
	rownames(mppd) = colnames(mppd) = states

	# Compute median pairwise patristic distance between state subtrees
	for(i in 1:ncol(state.pairs)) {
		st <- sort(state.pairs[,i])

		grp1 <- state.groups[[st[1]]]
		grp2 <- state.groups[[st[2]]]

		patdist <- array(NA, dim=length(grp1)*length(grp2))
		j <- 1

		for(g1 in grp1) {
			for(g2 in grp2) {
				k <- keyfunc(g1,g2)
				if(k %in% pdist$key) {
					patdist[j] <- pdist$distance[pdist$key == k]
					j <- j + 1
				} else {
					stop(paste('patristic distance not found for inter-leaf pair: ',g1,', ',g2,sep=''))
				}

			}
		}

		med <- median(patdist)
		mppd[st[1],st[2]] <- med
	}

	# Compute median pairwise patristic distance within each subtype
	for(i in 1:n) {

		st <- states[i]
		grp <- state.groups[[st]]
		ng <- length(grp)

		if(ng == 1) {
			mppd[st,st] <- NA

		} else {
			l <- ng*(ng-1)/2
			patdist <- array(NA, dim=l)
			j <- 1

			for(p in 1:(ng-1)) {
				for(q in (p+1):ng) {
					g1 <- grp[p]
					g2 <- grp[q]
				
					k <- keyfunc(g1,g2)
					if(k %in% pdist$key) {
						patdist[j] <- pdist$distance[pdist$key == k]
						j <- j + 1
					} else {
						stop(paste('patristic distance not found for intra-leaf pair: ',g1,', ',g2,sep=''))
					}
				}
			}

			med <- median(patdist)
			mppd[st,st] <- med
		}
	}

	mppd[lower.tri(mppd)] <- t(mppd)[lower.tri(mppd)]

	return(mppd)
}

transitions <- function(tree, subtypes, plot.name=NULL) {
# Identify subtrees where all patristic distances between leaves are below 
# specified threshold
#
# Args:
#  tree: phylo object containing subtype tree
#  patristic.distance: patristic distance threshold
#  plot.name[OPTIONAL]: Plot output filename
#
# Returns:
#   factor of tip.label names assigned to subtype
#

	ntips<-Ntip(tree)
	nnodes<-Nnode(tree)
	subtree.states <- list()
	subtree.distances <- list()
	states <- sort(levels(subtypes))

	# DFS tree
	igraph.tree <- graph.edgelist(tree$edge) # tree in igraph form
	dfs <- graph.dfs(igraph.tree,root=ntips+1,neimode='out', order=TRUE,dist=TRUE)

	set.count <- function(state, val=1) {
		counts <- rep(0,length(states))
		names(counts) <- states
		counts[state] = val

		return(counts)
	}

	# Travese the tree in post-order
	ord <- rev(dfs$order)
	for(i in 1:length(ord)) {
		node <- ord[i]
		cnode <- as.character(node)
    	
    	if(node <= ntips) {
    		# Leaf node
    		# Save state

    		tip.label = tree$tip.label[node]
    		st <- as.character(subtypes[tip.label])
    		subtree.states[[cnode]] = set.count(st)
    		edge <- tree$edge[,2] == node
    		if(sum(edge) != 1) {
    			stop(paste('Cant find correct edge for tip',node))
    		}
    		tj <- tree$edge.length[tree$edge[,2] == node]
    		subtree.distance[[cnode]] = set.count(st, tj)

    	} else {
    		# Internal node
    		# Record states

    		children = Children(tree, node)
    		if(!all(as.character(children) %in% names(subtree.states))) {
    			stop('Error in dfs traversal of tree')
    		}

    		edge <- tree$edge[,2] == node
    		tj <- 0
    		if(any(edge)) {
    			tj <- tree$edge.length[tree$edge[,2] == node]
    		}

    		node.states = rep(0,length(states))
    		node.distances = rep(tj, length(states))
			names(node.states) <- states
			names(node.distances) <- states
    		sapply(children, function(j) { 
    			cc = subtree.states[[as.character(j)]];
    			node.states <<- node.states + cc
    		})

    		subtree.states[[cnode]] = node.states
    	}
    }
    
    # subtree.labels = sapply(subtree.states, function(l) paste(sort(l), collapse=''))
    # o <- order(as.numeric(names(subtree.labels)))
    # subtree.labels <- subtree.labels[o]
    # plot(tree)
    # tiplabels(subtree.labels[1:ntips])
    # #nodelabels(subtree.labels[(ntips+1):(ntips+nnodes)])

    return(subtree.states)

}



plot.liks <- function(fit, index.state) {

	stopifnot(index.state %in% fit$states)

	# Covert to long format data
	df <- melt(as.data.frame(fit$lik.anc), id.vars=index.state)

	# Plot
	ggplot(data=df, aes_string(x=index.state, y="value", color="variable")) + 
    	geom_point() +
    	scale_fill_discrete()
	
}


isolate.rate <- function(states, p1, p2) {
	n <- length(states)
	qp <- matrix(1, nrow=n, ncol=n, dimnames=list(states, states))
	diag(qp) <- 0
	qp[p1,p2] = qp[p2,p1] = 2

	qp
}

transition.rate.parameters <- function(tree, subtypes, zero=TRUE, nbins=3:12) {
# Iteratively sets each pair of states as free model parameter and
# records transition rate.  Bins rates from first step into categories and 
# assigns rate-fitting parameter to each bin to create new transition matrix model.
#  
# In tests, symmetric matrices have improved predictive performance over equal rate
# models, however, in large subtype sets, 1) a symmetric model computation time is 
# prohibitive 2) many of the rates are similar.  This method attempts to achieve the 
# predictive performance gains while reducing the number of parameters compared to the 
# symmetric model.
#
# Args:
#  tree: phylo object containing subtype tree
#  subtypes: factor list of subtype assignments for tree$tip.label entries.
#  zero: Boolean, if TRUE set rate parameter as 0 for estimated rates that are zero instead of free parameter
#  nbins: Number of free parameters to assign to estimated rates
#
# Returns:
#   Q parameter matrix model
#

	states <- sort(as.character(levels(subtypes)))
	n <- length(states)
	state.pairs <- combn(states,2)
	rates <- cbind.data.frame(t(state.pairs), as.numeric(0), apply(state.pairs, 2, paste, collapse='__'), stringsAsFactors=FALSE)
	colnames(rates) <- c('state1','state2','rate','key')

	for(i in 1:ncol(state.pairs)) {
		s1 <- state.pairs[1,i]
		s2 <- state.pairs[2,i]
		model.params <- isolate.rate(states, s1, s2)
		fit <- fitMk(tree, subtypes, model=model.params)
		rate <- fit$rates[2]

		k <- paste(c(s1, s2), collapse='__')
		rates$rate[rates$key == k] <- rate
		#print(paste('completed iteration:',k))
	}

	cl <- Mclust(rates$rate, G=nbins)

	qp <- matrix(0, n, n)
	rownames(qp) = colnames(qp) = states
	qp[lower.tri(qp)] = cl$classification
	qp[upper.tri(qp)] = t(qp)[upper.tri(qp)]

	return(list('estimated.rates'=rates, 'clustering'=cl, 'model'=qp))
}

transition.rate.parameters2 <- function(tree, subtypes) {
# Sets all subtypes pairs with a small median patristic distance as a free parameter.
#  
# This method tries to identify a transition rate parameter formulation that has 
# improved performance over the equal rates model, but fewer parameters than the
# symmetric model.
#
# Args:
#  tree: phylo object containing subtype tree
#  subtypes: factor list of subtype assignments for tree$tip.label entries.
#
# Returns:
#   Q parameter matrix model
#

	states <- sort(levels(subtypes))
	n <- length(states)
	qp <- matrix(0, n, n)
	rownames(qp) = colnames(qp) = states

	bdist = branchDistances(tree, subtypes, plot.name=NULL)

	# Compute mppd
	mppds <- mppd(bdist, subtypes)
	intermppds <- mppds[lower.tri(mppds)]

	# Set all pairs in lower quantiles as free parameter
	# Leave pairs in higher quantiles as same parameter
	cutoffs <- quantile(intermppds, c(.25,.75))
	clusters <- rep(1, length(intermppds))
	clusters[intermppds >= cutoffs[2]] <- 2
	nfree <- sum(intermppds < cutoffs[1])
	#clusters[intermppds < cutoffs[1]] <- 3:(nfree+2)
	clusters[intermppds < cutoffs[1]] <- 3

	qp[lower.tri(qp)] = clusters

	qp[upper.tri(qp)] = t(qp)[upper.tri(qp)]
	
	return(list('model'=qp))
}


reassign.subtypes <- function(tree, subtypes) {
# Sometimes subtypes are mislabelled. This method finds subtrees that have
# max patristic distance less than a given threshold (indicating 0.4 probability
# distances equal to or greater then distance cutoff belong to intra-subtype distribution)
# and reassigns all tips in subtree as the max represented subtype in that tree
#
# Args:
#  tree: phylo object containing subtype tree
#  subtypes: factor list of subtype assignments for tree$tip.label entries.
#
# Returns:
# 	data.frame with 'new' and 'prev' subtype factors
#

	rs = subtype.diameter(tree, subtypes, 0.4)

	# Find subtrees that fall under this distance threshold
	print(rs$diameter)
	subt = assignSubtrees(rs$patristic.distance, tree, rs$diameter, plot.name=NULL)
	newsubt = subt

	for(s in levels(subt)) {
		grp <- names(subt)[subt == s]
		oldsubtypes = as.character(subtypes[grp])

		newsubtype = names(which(table(oldsubtypes) == max(table(oldsubtypes))))
		
		# Rename
		levels(newsubt)[levels(newsubt)==s] <- paste(newsubtype,collapse='_')

	}	

	ord <- names(newsubt)[order(newsubt)]
	return(data.frame(new=newsubt[ord], prev=subtypes[ord]))
}

subtype.diameter <- function(tree, subtypes, p) {
# Determine a distance cutoff based on the intra-subtype distance
# distribution for a given probability value (indicating likelihood distances 
# equal to or greater then distance cutoff belong to intra-subtype distribution)
#
# Args:
#  tree: phylo object containing subtype tree
#  subtypes: factor list of subtype assignments for tree$tip.label entries.
#  p: probability
#
# Returns:
# 	diameterquantile value for given probability
#
	bdist = branchDistances(tree, subtypes, plot.name=NULL)

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
	# file <- 'subtype_patristic_distribution_fit'
	# fn = file.path(output_dir, paste(file, '.png', sep=''))
	# graphics.off()
	# png(filename=fn)
	# plot(fit)
	# graphics.off()

	qmethod = paste('q',distr,sep='')
	param = c(p=p, as.list(fit$estimate), lower.tail=FALSE)
	dcutoff = do.call(qmethod, param)

	return(list(diameter=dcutoff,patristic.distances=bdist))
}


monophyletic.subtypes <- function(tree, subtypes) {
# Assigns each subtree under a given max patristic threshold (indicating 0.05 probability
# distances equal to or greater then distance cutoff belong to intra-subtype distribution)
# a unique subtype designation. New subtype names will be merged or altered names from
# old subtype naming scheme.
#
# Args:
#  tree: phylo object containing subtype tree
#  subtypes: factor list of subtype assignments for tree$tip.label entries.
#
# Returns:
# 	data.frame with 'new' and 'prev' subtype factors
#

	rs = subtype.diameter(tree, subtypes, 0.05)

	# Find subtrees that fall under this distance threshold
	subt = assignSubtrees(rs$patristic.distance, tree, rs$diameter, plot.name=NULL)

	# Generate new subtype names based old subtype naming scheme

	subsubtypes <- list()
	subsubtype_names <- as.character(1:26)
	if(grepl('\\d$', levels(subtypes)[1], perl=TRUE)) {
		subsubtype_names <- letters
	}
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

				grp_subtypes[os] = paste(os, subsubtype_names[st], sep='')
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

	ord <- names(subt)[order(subt)]
	return(data.frame(new=subt[ord], prev=subtypes[ord]))

}






