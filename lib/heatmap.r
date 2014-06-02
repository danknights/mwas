# assumes samples are in the rows
"get.per.group.clustering" <- function(x, category){
	require(vegan)
	category <- as.factor(as.character(category))
	column.order <- NULL
	for(level in levels(category)){
		ix <- category == level
		feature.ix <- colSums(x[ix,] > 0) > 0

		# cols are the samples, let's hclust within each treatment
# 		d <- vegdist(x[ix,feature.ix], method = "bray")
		d <- dist(x[ix,feature.ix])
		clus <- hclust(d, "aver")
		dg <- as.dendrogram(clus)
		column.order <- c(column.order, which(ix)[order.dendrogram(dg)])
	}
	return(column.order)
}

"make.heatmap" <- function(x,category,column.color.by=NULL,filename='heatmap.pdf')
{
	require(vegan)
	library(RColorBrewer)
	library(gplots)
	
	# Cluster the taxa/pathways all together
	data.dist <- vegdist(t(x), method = "bray")
	row.clus <- as.dendrogram(hclust(data.dist, "complete"))
	
	# cluster within each factor level		
	category <- as.factor(as.character(category))
	count <- 0
	column.order <- get.per.group.clustering(x, category)	
	hm.colors <- colorRampPalette(c("blue","black","yellow"))(27)

	# log and normalize after clustering	
	x <- log10(x + min(x[x>0])/2)
 	x <- sweep(x,2,colMeans(x),'-')
 	x <- sweep(x,2,apply(x,2,sd),'/')
	
	# additional sample color bar
	nbreaks <- 20
	# note: we use the transpose here since we want taxa as rows
	if(!is.null(filename)) pdf(file=filename,width=8.5,height=11)

	sepwidth=c(0,0)
	if(!is.null(column.color.by)){
		dw.colors <- rep('black',nrow(x))
		dw.color.palette <- colorRampPalette(c(brewer.pal(9,'Set1')[1],'snow3'))(nbreaks)
		color.ix <- as.numeric(cut(column.color.by,breaks = nbreaks))
		for(i in seq_along(column.color.by)){
			if(!is.na(column.color.by[i])){
				dw.colors[i] <- dw.color.palette[color.ix[i]]
			}
		}
		dw.colors <- dw.colors[column.order]
		obj <- heatmap.2(as.matrix(t(x))[,column.order],
			Rowv = row.clus,
			Colv = NA, dendrogram='none',
			col = hm.colors, ColSideColors = dw.colors, trace = "none", density.info = "none", 
			main = filename, margins = c(15,13),
			scale='none',
			key=TRUE,
			symbreaks=TRUE,
			sepwidth=sepwidth
		)
	} else {
		obj <- heatmap.2(as.matrix(t(x))[,column.order],
			Rowv = row.clus,
			Colv = NA, dendrogram='none',
			col = hm.colors, trace = "none", density.info = "none", 
			main = filename, margins = c(15,13),
			scale='none',
			key=TRUE,
			symbreaks=TRUE,
			sepwidth=sepwidth
		)
	}
	
	if(!is.null(filename)) dev.off()
	invisible(obj)
}
