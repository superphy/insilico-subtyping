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


fitSubtypeDistribution <- function(sdist, plot.name=NULL) {
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

	fit = Mclust(sdist$distance)

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

	# Find the smallest distribution
	comp = as.character(names(which.min(fit$parameters$mean)))

	# Check the properties of the distribution

	# Retrieve subset of tip pairs in same subtype
	same = sdist[fit$classification %in% comp,]

	return(list(fit=fit, same=same))
}


assignSubtypes <- function(sdist, tree, plot.name=NULL, verbose=FALSE) {
# Using tips that fall into subtype branch length distribution, assign 
# tree tips to subtypes
#
# Args:
#  sdist: data.frame with tip.labels pairs in same subtype (named $t1 and $t2)
#  tree: phylo object containing subtype tree
#  plot.name[OPTIONAL]: Plot output filename
#
# Returns:
#   factor of tip.label names assigned to subtype
#

	cgroup = 1
	tips = tree$tip.label
	subtypes = integer(length=length(tips))
	names(subtypes) = tips

	sdist = sdist[order(sdist$distance),]
	grouped = as.character(apply(sdist[,c('t1','t2')], 1, mykeyFunc))

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

	subtypes = factor(subtypes)
	if(!is.null(plot.name)) {
		graphics.off()
		png(filename=plot.name,
	        width=12*72,height=12*72
	    )
		phylotyper$plot.subtype(tree, subtypes)
		graphics.off()
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


