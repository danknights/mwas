# This function returns the next kegg level up from the current kegg level
# kegg: vector of kegg description named by current kegg level
get.next.kegg <- function(kegg)
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