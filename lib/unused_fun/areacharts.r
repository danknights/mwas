# plots an area chart of communities
# sorts by jensen-shannon divergence
# or by given order
# order.by = 'pearson', or indices, or NULL
# if colors is null, uses brewer Set1
"areachart" <- function(x,order.by=NULL, col=NULL){
	if(is.null(col)) {
		require('RColorBrewer')
		col <- brewer.pal(12,'Paired')
	}
	nr <- nrow(x)
	nc <- ncol(x)
	x <- x[,rev(order(colMeans(x)))]
	col <- rep(col,ceiling(nc / length(col)))
	col <- c(col[1:ncol(x)],'gray')

	# reorder rows	
	if(!is.null(order.by)){
		if(length(order.by) == 1){
			if(order.by == 'spearman'){
				d <- 2 - cor(t(x),method='spear')
				pc <- cmdscale(d,1)
				x <- x[order(pc),]
			} else if(order.by == 'mean'){
				x <- x[order(rowSums(x)),]
			}
		} else {
			x <- x[order.by,]
		}
	}

	x <- cbind(x,max(rowSums(x))-rowSums(x))

	# initialize
	plot(0,0,type='n',xlab='Samples',ylab='Relative abundance',
			xlim=c(0,nr),
			ylim=c(log10(min(x[x>0])/2),max(log10(rowSums(x + min(x[x>0])/2))))
# 			ylim=c(0,max(rowSums(x)))
			
			)
			
	# y holds the current max y-values for each sample
	y.coord <- rep(0,nr)
	y.coord <- rep(min(x[x>0])/2,nr)
	x.coord <- c(1:nr,nr:1)
	for(i in 1:nc){
# 		print(cbind(x.coord,log10(c(y.coord + x[,i],rev(y.coord)))))
		polygon(x.coord,log10(c(y.coord + x[,i],rev(y.coord))),col=col[i],border=NA)
		y.coord <- y.coord + x[,i]
	}
}


# beeswarm for multiple features
# makes a beeswarm for all columns of x, each a diff color
# if zero.color not null, uses that for zero values
"beeswarm.multi" <- function(x,category,col=NULL,zero.color=NULL,pt.alpha='99'
	){

	if(is.null(col)) {
		require('RColorBrewer')
		col <- brewer.pal(12,'Paired')
	}
	nr <- nrow(x)
	nc <- ncol(x)
	col <- rep(col,ceiling(nc / length(col)))
	col <- sprintf('%s%s',col,pt.alpha)
	# build a single large vector for beeswarm
	yplot <- NULL
	pwcol <- NULL
	xplot <- NULL
	for(i in 1:ncol(x)){
		yplot.i <- x[,i]
		pwcol.i <- rep(col[i],length(yplot.i))
		xplot.i <- category
		if(!is.null(zero.color)) pwcol.i[yplot.i == 0] <- zero.color
# 		keepix <- yplot.i > 0
# 		yplot.i <- yplot.i[keepix]
# 		xplot.i <- xplot.i[keepix]
# 		pwcol.i <- pwcol.i[keepix]
		yplot <- c(yplot,yplot.i)
		pwcol <- c(pwcol,pwcol.i)
		xplot <- c(xplot,xplot.i)
	}
	beeswarm(yplot ~ xplot,method='swarm',corral='random',pwcol=pwcol,pwbg=pwcol,pch=16)
	bxplot(yplot ~ xplot,add=T)

}
