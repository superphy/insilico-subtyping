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

phylotyper$makePriors = function(tree, x) {
	#
	#
	#
	#
	#

	
	# Find leaves with no subtype
	undefined = setdiff(tree$tip.label rownames(x))

	# Build prior matrix
	priors <- to.matrix(x,seq=levels(x))

	# Add flat priors for leaves with no subtype
	priors <- rbind(prior, matrix(1/length(levels(x)),
		length(undefined), 2, dimnames=list(undefined)))


	return(prior.matrix=priors, untyped=undefined)
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
	if((! is.null(libloc)) && libloc) {
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
	
	NULL
}

phylotyper$test_tree = function() {

	## simulate stochastic pure-birth tree
	tree<-pbtree(n=26,tip.label=LETTERS,scale=1)
	## generate character transition matrix
	Q<-matrix(c(-1,1,1,-1),2,2)
	rownames(Q)<-colnames(Q)<-letters[1:2]
	Q


}




while("phylotyper" %in% search())
  detach("phylotyper")
attach(phylotyper)

