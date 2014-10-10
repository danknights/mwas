# THIS FILE IS USED FOR TESTING THE HEATMAP LIB FUNCTIONS ONLY
# this is what an MWAS user code should look like for running mwas.heatmap
# otufile: name of file containing OTU Table or PICRUST pathway table
# mapfile: name of file containing mapping metadata
# keep.pathways.file: [pathway analysis only] file containing pathways to filter by (keep); can be any level (L1-L3)
# heatmap.title: title to display on the heatmap 
# outputfile: heatmap output file name 
# cluster.var: vector of column names to cluster samples by (e.g. can be a single column for simple heatmaps, such as Treatment)
# color.list: list of color vectors, one vector per color.var to specify the colors used
# constraints: optional. list of vectors, each named after a map column; vector values represent values to keep
"run.heatmap" <- function(otufile, mapfile)
{
	source("I-methods.r")
	source("pj.preprocess.r")
		
	# read in data files	
	map <-  load.qiime.mapping.file(mapfile)   # mapping file
	if (grep(".biom",mapfile)) {
		biom_table <- read_biom(otufile)         # OTU table - biom format
		otus <- t(as.matrix(biom_data(biom_table)))        # OTU table - classic format
	}
	else {
		trycatch(otus <- read.delim(otufile, sep='\t',
		comment='',head=T,row.names=1,check.names=F),error = function(err) 
			print("Couldn't parse OTU table. If BIOM format, use .biom extension"))
	}

	# intersect, filter, and transform
	# *** TODO replace this chunk with the real preprocess function when it's available
	otu <- preprocess(map=map, otu=otus)
	
	#TODO figure out kegg stuff later


	# burden the user with some custom subsetting 
	map <- subset(map, !is.na(Diabetes) & Location=="fecal" & Week==6 & Sex %in% c("M") & Treatment %in% c("Control","PAT"))
	otu <- otu[rownames(map),]

	# prior to testing for differentiation make sure to filter KEGG features

	# call the stats function to test for differentiated features
	# offload to user: create extra variable for combination vars if > 1
	# combine cluster variable values to form groups to cluster samples by
	new.treatments <- create.new.treatments(map, cluster.var=c("Sex", "Treatment"))
	diff.features <- differentiation.test(otu, category=new.treatments)$features

	# make the heatmap
	# TODO: replace this with the parent plot function once it's ready
	mwas.heatmap(otu, map, diff.features, cluster.var=c("Sex", "Treatment"), color.var=names(color.list), color.list, kegg_pathways=kegg_pathways, heatmap.title=heatmap.title, outputfile=outputfile)
}


# combines cluster variable values to create new groupings to cluster by
"create.new.treatments"<-function(map, cluster.var)
{
	if(length(cluster.var) > 1){
		new.treatments <- as.factor(apply(map[,cluster.var], 1, paste0, collapse="."))
	} else {
		new.treatments <- map[,cluster.var]
		# cluster.var must be a factor in order to cluster properly, apply original levels
		new.treatments <- as.factor(new.treatments)
		levels(new.treatments) <- names(color.list[[cluster.var]])
	}
	names(new.treatments) <- rownames(map)	
	new.treatments
}

