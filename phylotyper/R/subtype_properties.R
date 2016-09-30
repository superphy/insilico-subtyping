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
subtypeBranchDistances <- function(tree, subtypes, do.plot=TRUE) {
# Compute patristic distance between tips in same vs different subtype groups
#
# Args:
#  tree: phylo object containing subtype tree
#  subtypes: factor list of subtype assignments
#  do.plot: Generate plots
#
# Returns:
#     data.frame of same/different subtype classification and corresponding distance 
#

	d = cophenetic(tree)
	tips = tree$tip.label
	pairs = combn(tips,2)
	same = matrix(NA, ncol=ncol(pairs), nrow=nrow(pairs)+1)
	diff = matrix(NA, ncol=ncol(pairs), nrow=nrow(pairs)+1)
	si = 1
	di = 1

	for(c in 1:ncol(pairs)) {
		pair = pairs[,c]

		a = pair[1]
		b = pair[2]

		if(subtypes[a] == subtypes[b]) {
			same[,si] = c(a,b,d[a,b])
			si = si + 1
		}
		else {
			diff[,di] = c(a,b,d[a,b])
			di = di + 1
		}
	}

	same = same[,apply(same,2,function(x) all(!is.na(x)))]
	diff = diff[,apply(diff,2,function(x) all(!is.na(x)))]

	df = data.frame("subtype"=c(rep('same',ncol(same)), rep('different',ncol(diff))),
			"distance"=as.numeric(c(same[3,], diff[3,])),
			"t1"=as.character(c(same[1,], diff[1,])), "t2"=as.character(c(same[2,], diff[2,])))

	# Plot with density and rug
	if(do.plot) {
		ggplot(df, aes(x = distance, colour = subtype)) +
			geom_freqpoly(binwidth = .01)
	}

	return(df)
}


fitSubtypeDistribution <- function(sdist, do.plot=FALSE) {
# Identify the branch length distribution that correspond to 
# tips in the same subtype
#
# Args:
#  sdist: data.frame with branch distance (named $distance)
#
# Returns:
#   list:
#     fit = Mclust object
#     same = subset of sdist for tips that are in same subtype
#

	fit = Mclust(sdist$distance)

	if(do.plot) {
		layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
		plot(fit, what="density")
		plot(fit, what="BIC")
		plot(fit, what="uncertainty")
		plot(fit, what="classification")
	}

	# Find the smallest distribution
	comp = as.character(names(which.min(fit$parameters$mean)))

	# Check the properties of the distribution

	# Retrieve subset of tip pairs in same subtype
	same = sdist[fit$classification %in% comp,]

	return(list(fit=fit, same=same))
}


assignSubtypes <- function(sdist, tree, do.plot=FALSE, verbose=FALSE) {
# Using tips that fall into subtype branch length distribution, assign 
# tree tips to subtypes
#
# Args:
#  sdist: data.frame with tip.labels pairs in same subtype (named $t1 and $t2)
#  tree: phylo object containing subtype tree
#
# Returns:
#   factor of tip.label names assigned to subtype
#

	cgroup = 1
	tips = tree$tip.label
	subtypes = integer(length=length(tips))
	names(subtypes) = tips

	sdist = sdist[order(sdist$distance),]
	grouped = as.character(apply(sdist[,c(3,4)], 1, mykeyFunc))

	for(r in 1:nrow(sdist)) {
		t1 = as.character(sdist$t1[r])
		t2 = as.character(sdist$t2[r])
		s1 = as.character(subtypes[t1])
		s2 = as.character(subtypes[t2])

		if(verbose) print(paste('Branch:',t1,t2,', distance:',sdist[r,2]), sep=' ')
		
		if(s1 > 0) {
			# t1 already assigned
			if(s2 > 0) {
				# t2 already assigned

				if(s1 != s2) {
					# Disjoint subtypes
					# Check if they can be merged
					g1 = tips[subtypes == s1]
					g2 = tips[subtypes == s2]
					res = mergeable(grouped, g1, g2, subtypes, cgroup)

					if(res$result) {
						subtypes <- res$subtypes
						cgroup <- res$cgroup
						if(verbose) print(paste('merging',s1,'&',s2,'into',cgroup-1))
					} else {
						if(verbose) print(paste(s1, "&", s2, "groups cannot be merged. The intersection of groups do not form a clique."))
					}
				}

			} else {
				# t2 not assigned
				# Check if it belongs in group
				g1 = tips[subtypes == s1]
				res = mergeable(grouped, g1, t2, subtypes, cgroup)

				if(res$result) {
					subtypes <- res$subtypes
					cgroup <- res$cgroup
					if(verbose) print(paste('assigning',t2,'to',s1))
				} else {
					if(verbose) print(paste(t2, " cannot be added to ",s1,". It does not form a clique."))
				}

			}

		} else {
			# t1 not assigned

			if(s2 > 0) {
				# t2 already assigned
				# Check if it belongs in group
				g2 = tips[subtypes == s2]
				res = mergeable(grouped, t1, g2, subtypes, cgroup)

				if(res$result) {
					subtypes <- res$subtypes
					cgroup <- res$cgroup
					if(verbose) print(paste('assigning',t1,'to',s2))
				} else {
					if(verbose) print(paste(t1, "cannot be added to ",s2,". It does not form a clique."))
				}
				

			} else {
				# t2 not assigned
				# Create new group
				subtypes[t1] = cgroup
				subtypes[t2] = cgroup
				cgroup = cgroup + 1
				if(verbose) print(paste('creating new group',cgroup-1))

			}

		}
	}

	# Assign remaining to individual subtype groups
	for(i in 1:length(subtypes)) {
		if(subtypes[i] == 0) {
			subtypes[i] = cgroup
			cgroup = cgroup + 1
		}
	}

	# Rename
	st = levels(factor(subtypes))
	for(i in 1:length(st)) {
		subtypes[subtypes == st[i]] = i
	}

	return(factor(subtypes))

}


mergeable = function(grouped, g1, g2, subtypes, cgroup) {
	# Internal function

	for(t1 in g1) {
		for(t2 in g2) {
			if(!(mykeyFunc(c(t1,t2)) %in% grouped)) {
				return(list(result=FALSE))
			}
		}
	}

	# All members of two groups are linked
	for(t1 in g1) {
		subtypes[t1] = cgroup
	}
	for(t2 in g2) {
		subtypes[t2] = cgroup
	}

	return(list(result=TRUE, subtypes=subtypes, cgroup=cgroup+1))
}

mykeyFunc = function(trow) {
	# Internal function

	return(paste(as.character(sort(trow)),collapse='__'))
}


