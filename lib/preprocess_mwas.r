# Preprocessing the original OTU table 
#
# Contributors: PJ
# -----
# Input: 
#               otu : OTU table 
#   minOTUInSamples : 
# ----
# Output:
#      otu : filtered otu table
# ----- 
# Last update: 10/27/2014
#

# temporary place holder for use by PJ until the real preprocess function is written
"preprocess" <- function(otu, minOTUInSamples=.001, filter.kegg=FALSE)
{
	source("util.load.r")
		
	otu <- (remove.nonoverlapping.samples(map = map, otus = otus))$otus

	otu <- sweep(otu, 1, rowSums(otu), '/')
	otu <- otu[, colMeans(otu) > minOTUInSamples, drop=FALSE]
	otu <- asin(sqrt(otu))
	ret <- collapse.by.correlation(otu, .95)
	otu <- otu[, ret$reps]
	
	if(filter.kegg)
	{
		keep.pathways.file <- parse.params("preprocess", "keggfilterlist")
		kegg_pathways <- NULL
		# filter the kegg pathways by what's in the file
		kegg <- filter.pathways(kegg, keep.pathways.file)
		# return the next level up for displaying
		next.kegg <- get.next.kegg(kegg)
		names(next.kegg) <- names(kegg)
		kegg_pathways<-next.kegg
	}

	return(otu)
}

# This function filters the kegg descriptions by a vector of unknown-level kegg pathways
# kegg: vector of kegg descriptions named by current kegg level
# TODO: mention we havent implemented specifying current level filtering
"filter.pathways" <- function(kegg, keep.pathways.file)
{
	#filter by comparing to the entire database
#	kegg.db.file<- paste(Sys.getenv("MWAS.HEATMAP"),'/data/ko.to.pathways.txt', sep='')	
#	kegg.db <- read.delim(kegg.db.file, header=T)

	pathways.to.keep <- read.delim(keep.pathways.file, header=F)[,1]
		
	keggnames <- names(kegg)
	kegg <- as.character(kegg)
	names(kegg) <- keggnames
	# split the kegg descriptions of this pathway table into L1, L2, L3
	keggsplit <- lapply(kegg, function(xx) unlist(strsplit(xx, split="; ", fixed=TRUE)))
	keggmatrix <- do.call(rbind, keggsplit)
	
	# remove any L1 pathways that begin with Unclassified;
	kegg <- kegg[keggmatrix[,1]!="Unclassified"]
	keggmatrix <- keggmatrix[keggmatrix[,1]!="Unclassified",]
	kegg.level <- 0
	for(i in 3:1){
		if(sum(keggmatrix[,i] %in% pathways.to.keep) > 0)
			kegg.level <- i
 	}
	filtered.kegg <- kegg[keggmatrix[, kegg.level] %in% pathways.to.keep]

	as.factor(filtered.kegg)
}

# This function returns the next kegg level up from the current kegg level
# kegg: vector of kegg description named by current kegg level
"get.next.kegg" <- function(kegg)
{
	# split the kegg descriptions of this pathway table into L1, L2, L3
	keggsplit <- lapply(as.character(kegg), function(xx) unlist(strsplit(xx, split="; ", fixed=TRUE)))
	keggmatrix <- do.call(rbind, keggsplit)
	kegg.level <- 0
	for(i in 3:1){
		if(identical(keggmatrix[,i], names(kegg)))
			kegg.level <- i
 	}
 	if(kegg.level == 1)
 		stop("ERROR: current KEGG level cannot be 1")
 	if(kegg.level == 0)
 		stop("ERROR: check your kegg vector names. level not found.")
 		
 	kegg.level.up.one <- kegg.level-1
	kegg_pathways <- keggmatrix[, kegg.level.up.one]
	names(kegg_pathways) <- rownames(kegg)	
	kegg_pathways
}