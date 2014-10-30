# THIS FILE IS USED FOR TESTING THE HEATMAP LIB FUNCTIONS ONLY
# this is what an MWAS user code should look like for running mwas.heatmap
#
# ----------------
# Contributors: PJ
# ----------------
# Input: 
#              otufile: name of file containing OTU Table or PICRUST pathway table
#              mapfile: name of file containing mapping metadata
#   keep.pathways.file: [pathway analysis only] file containing pathways to filter by (keep); can be any level (L1-L3)
#        heatmap.title: title to display on the heatmap 
# -------
# output:
#        outputfile: heatmap output file name 
#       cluster.var: vector of column names to cluster samples by (e.g. can be a single column for simple heatmaps, such as Treatment)
#        color.list: list of color vectors, one vector per color.var to specify the colors used
#       constraints: optional. list of vectors, each named after a map column; vector values represent values to keep
# -------
#  Last updates: 10/29/2014
#

"run.heatmap" <- function(otufile, mapfile)
{
	source("util.load.r")
	source("preprocess_mwas.r")
		
	# read in data files	
	map <-  load.qiime.mapping.file(mapfile)   # mapping file
	results <- load.qiime.otu.table(otufile, include.lineages=TRUE) # flag to save kegg column
	otu <- results$otu
	kegg <- setNames(results$lineages, rownames(otu))

	# intersect, filter, and transform (AND filter kegg pathways here)
	# *** TODO replace this chunk with the real preprocess function when it's available
	otu <- preprocess(map=map, otu=otu, filter.kegg=TRUE)
	
	# burden the user with some custom subsetting 
	map <- subset(map, !is.na(Diabetes) & Location=="fecal" & Week==6 & Sex %in% c("M") & Treatment %in% c("Control","PAT"))
	otu <- otu[rownames(map),]

	# call the stats function to test for differentiated features
	# offload to user: create extra variable for combination vars if > 1
	# combine cluster variable values to form groups to cluster samples by
	new.treatments <- create.new.treatments(map, cluster.var=c("Sex", "Treatment"))
	diff.features <- differentiation.test(otu, category=new.treatments)$features

	# make the heatmap
	# TODO: replace this with the parent plot function once it's ready
	mwas.heatmap(otu, map, diff.features, cluster.var=c("Sex", "Treatment"), color.var=names(color.list), color.list, kegg_pathways=kegg_pathways, heatmap.title=heatmap.title, outputfile=outputfile)
}
