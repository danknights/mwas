# Creates a heatmap to display differentiated kegg pathways, with additional color bars at the top
# NOTE: uses pearson's correlation as distance function due to inherent limitations with Bray-Curtis
# -------------
# Contributor: PJ
# -------------
# Input:
#         otu: normalized and filtered otu table
#         map: mapping table that has already been subsetted to include samples of interest
# cluster.var: vector of column names to cluster samples by (e.g. can be a single column for 
#              simple heatmaps, such as Treatment)
#   color.var: vector of column names, color bars added to top of heatmap (suggested usage: 
#              set color.var=cluster.var) note that factors will be treated as distinct colors, 
#              and numerics will be transformed into gradients
#  color.list: list of color vectors, one vector per color.var to specify the colors used
# kegg_pathways: vector of pathways, useful for highlighting different levels
# heatmap.title: title to display on the heatmap 
# -------------
# Output:
#   outputfile: heatmap output file name 
# -------------
# Last update:10/29/2014
#

# hclusts the columns of an OTU table according to groupings of samples
cluster.columns<-function(otu, distfun, new.treatments)
{
	col.labels <- NULL
	# separate samples by unique values of new.treatments that we're interested in clustering by
	unique.new.treatments <- levels(new.treatments)
	for(i in 1:length(unique.new.treatments)){
		otu.by.clustergroup <- otu[new.treatments==unique.new.treatments[i], ]
		data.dist.group <- distfun(otu.by.clustergroup)
		col.clus.group <- hclust(data.dist.group)
		d <- as.dendrogram(col.clus.group)
		# since the index order is per group only, we'll have to add the current length each time 
		# to get correct indices along the combined dendrogram
		col.labels <- c(col.labels, labels(d))
	}
	
	col.labels
}
# hclusts the rows 
cluster.rows<-function(otu, distfun)
{
	# do the row-clustering (note that we need to explicitly reorder the dendrogram here)
	Rowv <- rowMeans(otu)
	hcr <- hclust(distfun(otu))
	ddr <- as.dendrogram(hcr)
	ddr <- reorder(ddr, Rowv)
	rowInd <- order.dendrogram(ddr)
#		row.labels <- labels(ddr)[rowInd] #THIS IS CRITICAL?
	row.labels <- labels(ddr)
#print(labels(ddr)[rowInd])

	list(ddr=ddr, row.labels=row.labels)
}

#"shorten.taxonomy" <- function(ids,delim=';'){
#	ids <- gsub('[kpcofgs]__','',ids)
#	newids <- ids
#	ids <- strsplit(ids,delim)
#	for(i in seq_along(ids)){
#		n <- length(ids[[i]])
#		j <- n
#		while(ids[[i]][j] == 'Other' || ids[[i]][j] == '') j <- j - 1
#		newids[i] <- ids[[i]][j]
#	}
# 	return(newids)
#}


# map: rows must have already been reorder by the hclust
# color.list: list of named vectors that have been named according to color.var values
create.color.bars<-function(map, color.var, color.list)
{
	col.colors <- NULL
	for(i in 1:length(color.var))
	{
		reordered <- map[,color.var[i]]
		if(is.numeric(reordered)){
			color.palette <- colorRampPalette(color.list[[color.var[i]]])
			num.breaks <- 20 #? for gradients
			colors <- color.palette(num.breaks)[as.numeric(cut(reordered,breaks = num.breaks))]
		} else {
				colors <- unname(color.list[[color.var[i]]][as.character(reordered)])
		}
		col.colors <- cbind(col.colors, colors)		
	}
	colnames(col.colors) <- color.var
	col.colors

}

create.new.treatments<-function(map, cluster.var)
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

# Creates a heatmap to display differentiated kegg pathways, with additional color bars at the top
# NOTE: uses pearson's correlation as distance function due to inherent limitations with Bray-Curtis
#
# otu: normalized and filtered otu table
# map: mapping table that has already been subsetted to include samples of interest
# cluster.var: vector of column names to cluster samples by (e.g. can be a single column for simple heatmaps, such as Treatment)
# color.var: vector of column names, color bars added to top of heatmap (suggested usage: set color.var=cluster.var)
#			note that factors will be treated as distinct colors, and numerics will be transformed into gradients
# color.list: list of color vectors, one vector per color.var to specify the colors used
# kegg_pathways: vector of pathways, useful for highlighting different levels
# heatmap.title: title to display on the heatmap 
# outputfile: heatmap output file name 
heatmap.mwas <- function(otu, map, diff.features, cluster.var, color.var, color.list, kegg_pathways, heatmap.title="", outputfile)
#TODO add some extra formatting parameters here
# OPTIONAL: distance matrix
{
	source("heatmap.3.r")
    #source("util.r")
	
    #require(vegan)
    #require(RColorBrewer)
    #require(gplots)
	
	use.kegg <- as.logical(length(kegg_pathways))
	
	otu <- otu[,diff.features]	
	
	# use pearson's correlation for the distance function (bray-curtis has scaling issues)
	distfun <- function(x) as.dist((1-cor(t(x)))/2)
	
	# do row and column clustering
	col.labels <- cluster.columns(otu=otu, distfun=distfun, new.treatments=new.treatments)
	ret <- cluster.rows(otu=as.matrix(t(otu))[,col.labels], distfun=distfun)
	row.dendrogram <- ret$ddr 
	row.labels <- ret$row.labels
	
	# create each of the colored bars, reorder to previous column clustering
	col.colors <- create.color.bars(map[col.labels,], color.var, color.list)

	# color row axis labels according to KEGG
	kegg.colors <- NULL
	if(use.kegg)
	{
		# TODO add more colors here
		axis.colors <- c("orange","firebrick3","aquamarine3","steelblue3","black", "mediumpurple3")
		kegg <- kegg_pathways[row.labels]
		lookup <- axis.colors[1:length(unique(kegg))]
		names(lookup) <- sort(unique(kegg))
		# diabetes only!
		#lookup["Metabolism of Other Amino Acids"] <- lookup["Amino Acid Metabolism"] 
		kegg.colors <- lookup[as.character(kegg)] 
		names(kegg.colors) <- names(kegg) #change these names so it's easier to reference later
	}	
	
	# do all formatting here
	## consider redrawing on a really large PDF - extra spacing is okay
	## allow for max cap on # of features --> set this from the command line (with disclaimer)
	
	lwid = c(1.9,9,4) #col widths
	lhei = c(1.1,.3+.1*length(color.var),7.1,1.6) #row heights
	lmat = rbind(c(3,4,0), c(0,1,0), c(0,2,0), c(0,5,0))
	margins=c(4,6)
	legend1.x=.601 #color bar legend 
	legend1.y=1.011
	width=7 #pdf size
	height=7 
	hm.colors <- colorRampPalette(c("blue","black","yellow"))(27)
	main <- heatmap.title
	
	if(!use.kegg) #update some params if showing only taxa
	{
		lhei = c(.7,.3,8,1) #row heights
		margins=c(.1,.1)
		legend1.x=.735
		legend1.y=1
		width=10
		height=12
	}
	
	# standardize each feature so it stands out more		
	otu <- apply(otu, 2, function(x) (x-mean(x))/sd(x))
	
	# make the heatmap
	pdf(file=outputfile, width=width, height=height)	

	otu.for.heatmap <- as.matrix(t(otu))
	if(!use.kegg){
		rownames(otu.for.heatmap) <- shorten.taxonomy(rownames(otu.for.heatmap))
	}
	
	heatmap.3(otu.for.heatmap[,col.labels], Rowv = row.dendrogram, Colv = NA, dendrogram='none',
	col = hm.colors, ColSideColors = col.colors, trace = "none", density.info = "none", rowAxisColors=kegg.colors,
	main = main, key=TRUE, keysize=5, lmat=lmat, lwid=lwid, lhei=lhei, labCol=NA, 
	margins=margins,cexCol=.7, scale="none"
	)

	# add the legend for the colored bars
	legend(x=legend1.x, y=legend1.y,legend=unname(unlist(lapply(color.list, names))), 
	fill=unlist(color.list),border=FALSE, bty="n", y.intersp = 0.8, cex=0.7, ncol=3, text.width=c(0, .04, .06 ) )#
																						#= strwidth("1,000")
	if(use.kegg) # add the legend for the row axis colors
	{
		legend(x=.59, y=.12,legend=names(lookup), 
		fill=lookup, border=FALSE, bty="n", y.intersp = 0.8, cex=0.7, text.col="dimgray")
	}

	dev.off()
		
}