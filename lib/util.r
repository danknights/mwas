# loads experimental data into a list
# required is map
# otu table and taxa.fp are optional
# min.presence is the fraction of samples in which a feature must be present to keep
# set collapse.threshold 0..1 to collapse taxa/otus by correlation (lower threshold
#   means more collapsing)
"load.experiment" <- function(map.fp, otu.fp=NULL, taxa.fp=NULL,
		alpha.fp=NULL, beta.fp=NULL, min.presence=0.1, collapse.threshold=NULL){

	otus <- NULL
	taxa <- NULL
	adiv <- NULL
	bdiv <- NULL
	pcoa <- NULL
	
	map <- read.table(map.fp,sep='\t',head=T,row=1,check=F,comment='')

	if(!is.null(otu.fp)){
		otus <- read.table(otu.fp,sep='\t',head=T,row=1,check=F,comment='',skip=1)
		lineages <- otus[,'taxonomy']
		otus <- otus[,-ncol(otus)]
		otus <- t(otus)	
	}

	if(!is.null(taxa.fp)){
		taxa <- list()
		for(i in seq_along(taxa.fp)){
			taxa[[taxa.fp[i]]] <- t(read.table(taxa.fp[i],sep='\t',head=T,row=1,check=F,comment=''))
		}
	}
	
	# keep only overlapping samples
	sample.names <- rownames(map)
	if(!is.null(otus)) sample.names <- intersect(sample.names, rownames(otus))
	if(!is.null(taxa)) sample.names <- intersect(sample.names, rownames(taxa[[1]]))
	map <- droplevels(map[sample.names,])

	# drop rare taxa
	if(!is.null(otus)) {
		otus <- otus[sample.names,]
		lineages <- lineages[colMeans(otus > 0) > min.presence]
		otus <- otus[,colMeans(otus > 0) > min.presence]
	}
	if(!is.null(taxa)) {	
		for(i in seq_along(taxa)){
			taxa[[i]] <- taxa[[i]][sample.names,]
			taxa[[i]] <- taxa[[i]][,colMeans(taxa[[i]] > 0) > min.presence]
		}
	}
	
	# collapse OTUs, taxa 
	if(!is.null(collapse.threshold)){
		if(!is.null(otus)){
			otus <- otus[,collapse.by.correlation(otus, min.cor=collapse.threshold)$reps]
		}
		if(!is.null(taxa)){
			for(i in seq_along(taxa)){
				taxa[[i]] <- taxa[[i]][,collapse.by.correlation(taxa[[i]], min.cor=collapse.threshold)$reps]
			}
		}
	}

	# load bdiv, adiv
	if(is.null(alpha.fp)){
		adiv <- NULL
	} else {
		adiv <- read.table(alpha.fp,sep='\t',head=T,row=1,check=F)
		adiv <- adiv[sample.names,]
	}
	
	if(is.null(beta.fp)){
		bdiv <- NULL
	} else {
		bdiv <- list()
		pcoa <- list()
		for(i in seq_along(beta.fp)){
			bdiv[[beta.fp[i]]] <- read.table(beta.fp[i],sep='\t',head=T,row=1,check=F)
			bdiv[[beta.fp[i]]]  <- bdiv[[beta.fp[i]]][sample.names,sample.names]
			pcoa[[beta.fp[i]]] <- cmdscale(bdiv[[beta.fp[i]]],3)
		}
	}
			
	return( list(
		map=map,
		otus=otus,
		taxa=taxa,
		bdiv=bdiv,
		adiv=adiv,
		pcoa=pcoa
	))
}



# linear test
"linear.test" <- function(x, y, covariates=NULL){
	if(!is.null(covariates)){
		covariates <- as.data.frame(covariates)
		covariates <- cbind(y, covariates)
		covariates <- droplevels(covariates)
		design <- model.matrix(~ ., data=covariates)		
	} else {
		design <- model.matrix(~y)
	}
	pvals <- apply(x, 2, function(xx) summary(lm(xx ~ ., data=covariates))[[4]][2,4])
	return(pvals)
}

"t.test.wrapper" <- function(x, y, use.fdr=TRUE){
	y <- as.factor(y)
	ix1 <- y == levels(y)[1]
	pvals <- apply(x,2,function(xx) t.test(xx[ix1],xx[!ix1])$p.value)
	if(use.fdr) pvals <- p.adjust(pvals,'fdr')
	return(pvals)
}


# returns vector of cluster ids for clusters with internal
# complete-linkage correlation of min.cor
"cluster.by.correlation" <- function(x, min.cor=.5){
#     library('fastcluster')
    cc <- cor(x,use='pairwise.complete.obs',method='pear')
    if(ncol(x) == 379) browser()
    cc <- as.dist(1-cc)
    hc <- hclust(cc)
    res <- cutree(hc,h=1-min.cor)
    names(res) <- colnames(x)
    return(res)
}

# returns vector of cluster ids for clusters with internal
# complete-linkage correlation of min.cor
#
# by default, chooses cluster reps as highest-variance member
# if select.rep.fcn=mean
"collapse.by.correlation" <- function(x, min.cor=.5, select.rep.fcn=c('var','mean','lowest.mean',
			'longest.name', 'shortest.name')[2],
        verbose=FALSE){
    if(verbose) cat('Clustering',ncol(x),'features...')
    gr <- cluster.by.correlation(x, min.cor=min.cor)
    if(verbose) cat('getting means...')
    if(select.rep.fcn == 'mean'){
        v <- apply(x,2,function(xx) mean(xx,na.rm=TRUE))
    } else if(select.rep.fcn == 'lowest.mean'){
        v <- apply(x,2,function(xx) -mean(xx,na.rm=TRUE))
    } else if(select.rep.fcn == 'longest.name'){
        v <- nchar(colnames(x))
    } else if(select.rep.fcn == 'shortest.name'){
        v <- -nchar(colnames(x))
    } else {
        v <- apply(x,2,function(xx) var(xx,use='complete.obs'))
    }
    if(verbose) cat('choosing reps...')
    reps <- sapply(split(1:ncol(x),gr),function(xx) xx[which.max(v[xx])])
    if(verbose)
        cat(sprintf('collapsed from %d to %d.\n',ncol(x), length(reps)))
    return(list(reps=reps, groups=gr))
}


"inspect.env" <- function(){
	objs <- ls(1)
	objs <- objs[!is.null(objs)]
	res <- split(objs, sapply(objs, function(xx) class(get(xx))))
	for(i in seq_along(res)){
		cat('\n')
		cat(names(res)[i], ':\n',sep='')
		print(res[[i]])
	}
	invisible(res)
}


"shorten.taxonomy" <- function(ids,delim=';',num.levels=1,must.include.level=7){
	ids <- gsub('[kpcofgs]__$','Other',ids)
	ids <- gsub('[kpcofgs]__','',ids)
	newids <- ids
	ids <- strsplit(ids,delim)
	for(i in seq_along(ids)){
		n <- length(ids[[i]])
		j <- n
		while(ids[[i]][j] == 'Other' || ids[[i]][j] == '') j <- j - 1
		start.level <- min(must.include.level,j-num.levels+1)
		start.level <- max(1,start.level)
		newids[i] <- paste(ids[[i]][start.level:j],collapse=' ')
		if(j < n) newids[i] <- paste('Uncl. ',newids[i],sep='')
	}
	
	# add indices to duplicate names
	counts <- table(newids)
	if(max(counts) > 1){
		for(name in names(counts)){
			if(counts[name] > 1){
				ix <- which(newids == name)
				for(i in seq_along(ix)){
					newids[ix[i]] <- paste(name,i,sep=' ')
				}
			}
		}
	}
	
	return(newids)
}

#QQ plot
# http://GettingGeneticsDone.blogspot.com/
# See http://gettinggeneticsdone.blogspot.com/p/copyright.html
 
# Define the function
qqplot.pvals	 = function(pvector, main=NULL, ...) {
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    plot(e,o,pch=19,cex=1, main=main, ...,
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        xlim=c(0,max(e)), ylim=c(0,max(o)))
    lines(e,e,col="red")
}



# Get balanced folds where each fold has close to overall class ratio
"balanced.folds" <- function(y, nfolds=10){
    folds = rep(0, length(y))
	y <- as.factor(y)
    classes = levels(y)
    # size of each class
    Nk = table(y)
    # -1 or nfolds = len(y) means leave-one-out
    if (nfolds == -1 || nfolds == length(y)){
        invisible(1:length(y))
    }
    else{
    # Can't have more folds than there are items per class
    nfolds = min(nfolds, max(Nk))
    # Assign folds evenly within each class, then shuffle within each class
        for (k in 1:length(classes)){
            ixs <- which(y==classes[k])
            folds_k <- rep(1:nfolds, ceiling(length(ixs) / nfolds))
            folds_k <- folds_k[1:length(ixs)]
            folds_k <- sample(folds_k)
            folds[ixs] = folds_k
        }
        invisible(folds)
    }
}



