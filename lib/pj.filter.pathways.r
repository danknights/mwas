# This function filters the kegg descriptions by a vector of unknown-level kegg pathways
# kegg: vector of kegg descriptions named by current kegg level
# specifically mention we havent implemented specifying current level filtering
filter.pathways <- function(kegg, keep.pathways.file)
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