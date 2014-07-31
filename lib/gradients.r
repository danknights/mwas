# Shows a PCoA plot colored by taxon gradients
# x is the taxon
# pc contains the two pc axes
# filename: if not NULL, must be a string of the pdf output file name
# scale.axes: e.g. c(1,-1) means to flip the second axis, not the first
"show.gradients" <- function(x, pc, filename=NULL,
                                incl.legend=FALSE,
                                pt.alpha='FF',bin.type=c('raw','log','quantile')[3],
                                title.text='',
                                axis.labels=c('PC1','PC2')){
    gradcols <- rev(rainbow(100,start=0,end=2/3))
    gradcols <- sprintf('%s%s',sapply(gradcols,function(x) substring(x,1,7)),pt.alpha)
	# reorder points randomly for drawing order
    scrambleix <- sample(nrow(pc))
    x <- x[scrambleix]
    pc <- pc[scrambleix, ]
    
    pdf.width <- 5
    if(incl.legend) pdf.width <- 6.5
    if(!is.null(filename)) pdf(filename, width=5, height=5)
    
    par(cex.main=.8)

    xlim <- range(pc[,1]) * 1.05
    ylim <- range(pc[,2]) * 1.05
    if(incl.legend) xlim[2] <- xlim[2] + diff(range(xlim)) * .25
	newx <- x

    if(bin.type=='quantile'){
	    # quantiles require 11 gradient points, to include 0%, 10%, ... 100%
        legend.gradcols <- rev(rainbow(11,start=0,end=2/3))
    	quantiles <- quantile(as.numeric(newx),probs=seq(0,1,length=length(gradcols)))
        gradixs <- sapply(newx, function(xx) min(which(quantiles >= xx)))
	} else if(bin.type=='log'){
		legend.gradcols <- rev(rainbow(10,start=0,end=2/3))
		newx[newx==0] <- min(newx[newx>0])/2
		newx <- log(newx,10)
		gradixs <- get.gradient.ixs(newx,length(gradcols))
	} else if(bin.type=='raw'){
		legend.gradcols <- rev(rainbow(10,start=0,end=2/3))
		gradixs <- get.gradient.ixs(newx,length(gradcols))
	}
    plot(pc[,1], pc[,2],type='p',pch=21,
                 col='#00000099', bg=gradcols[gradixs],
                 xlab=axis.labels[1],ylab=axis.labels[2],main=title.text,
                 xlim=xlim, ylim=ylim)
    if(incl.legend){
		legend.pos <- best.legend.location(pc[,1],pc[,2])
		if(bin.type=='quantile'){
			quantiles <- quantile(as.numeric(newx),probs=seq(0,1,length=length(legend.gradcols)))
			legend.text <- sprintf('%.3f', quantiles)
		} else {
			legend.text <- sprintf('%.2f',seq(min(newx), max(newx),length=length(legend.gradcols)))
		}
		legend.pos <- 'right'
		legend(legend.pos, legend.text, fill=legend.gradcols, cex=.65, y.intersp=1)
	}

    if(!is.null(filename)) dev.off()
}



# Shows a PCoA plot colored by taxon gradients
# x is a factor to color by
# pc contains the two pc axes
# filename: if not NULL, must be a string of the pdf output file name
# scale.axes: e.g. c(1,-1) means to flip the second axis, not the first
"show.metadata" <- function(y, pc, filename=NULL,
                                incl.legend=TRUE,
                                pt.alpha='FF',
                                title.text='',
                                axis.labels=c('PC1','PC2')){
	library('RColorBrewer')
	y.colors <- c(brewer.pal(9,'Set1'),brewer.pal(9,'Pastel1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))
	y.colors <- sprintf('%s%s',y.colors,pt.alpha)
	
	# reorder points randomly for drawing order
    scrambleix <- sample(nrow(pc))
    y <- y[scrambleix]
    pc <- pc[scrambleix, ]
    
    y <- as.factor(as.character(y))
    
    pdf.width <- 5
    if(incl.legend) pdf.width <- 6.5
    if(!is.null(filename)) pdf(filename, width=5, height=5)
    color.ix <- as.numeric(y)
    if(max(color.ix) > length(y.colors)) stop(paste('y has more than',y.colors,'values.',sep=' '))
    par(cex.main=.8)

    xlim <- range(pc[,1]) * 1.05
    ylim <- range(pc[,2]) * 1.05
    if(incl.legend) xlim[2] <- xlim[2] + diff(range(xlim)) * .25

    plot(pc[,1], pc[,2],type='p',pch=21, cex=1.5,
                 col='#00000099',bg=y.colors[color.ix],
                 xlab=axis.labels[1],ylab=axis.labels[2],main=title.text,
                 xlim=xlim, ylim=ylim)
    if(incl.legend){
		legend.pos <- 'right'
		legend(legend.pos, legend=levels(y), fill=y.colors, cex=.65, y.intersp=1)
	}

    if(!is.null(filename)) dev.off()
}


# converts real values to indices between 1 and n 
"get.gradient.ixs" <- function(x,n=10){
    stdx <- x-min(x)
    stdx <- stdx/abs(diff(range(stdx)))
    gradientix <- 1 + round((n-1)*stdx)
    return(gradientix)
}


# returns the best location from ('bottomleft','bottomright','topleft','topright')
"best.legend.location" <- function(x,y,xlim=NULL,ylim=NULL){
    if(is.null(xlim)) xlim <- range(x)
    if(is.null(ylim)) ylim <- range(y)
    
    # get boundaries for a 3x3 grid
    xbounds <-seq(xlim[1],xlim[2],length=4)
    ybounds <-seq(ylim[1],ylim[2],length=4)
    density <- numeric(4)
    # bottom left
    density[1] <- sum(x < xbounds[2] & y < ybounds[2])
    # bottom right
    density[2] <- sum(x > xbounds[3] & y < ybounds[2])
    # top left
    density[3] <- sum(x < xbounds[2] & y > ybounds[3])
    # top right
    density[4] <- sum(x > xbounds[3] & y > ybounds[3])
    return(c('bottomleft','bottomright','topleft','topright')[which.min(density)])
}
