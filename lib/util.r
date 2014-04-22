"load.experiment" <- function(map.fp, otu.fp=NULL, taxa.fp=NULL,
		alpha.fp=NULL, beta.fp=NULL, min.avg=0.1){
	otus <- read.table(otu.fp,sep='\t',head=T,row=1,check=F,comment='',skip=1)
	lineages <- otus[,'taxonomy']
	otus <- otus[,-ncol(otus)]
	otus <- t(otus)
	taxa <- list()
	for(i in seq_along(taxa.fp)){
		taxa[[taxa.fp[i]]] <- t(read.table(taxa.fp[i],sep='\t',head=T,row=1,check=F,comment=''))
	}
	
	map <- read.table(map.fp,sep='\t',head=T,row=1,check=F,comment='')
	ix <- intersect(rownames(map),rownames(otus))
	map <- droplevels(map[ix,])
	otus <- otus[ix,]

	# drop rare taxa
	min.avg = 0.1
	lineages <- lineages[colMeans(otus > 0) > min.avg]
	otus <- otus[,colMeans(otus > 0) > min.avg]
	for(i in seq_along(taxa)){
		taxa[[i]] <- taxa[[i]][ix,]
		taxa[[i]] <- taxa[[i]][,colMeans(taxa[[i]] > 0) > min.avg]
	}

	# load bdiv, adiv
	adiv <- read.table(alpha.fp,sep='\t',head=T,row=1,check=F)
	bdiv <- list()
	for(i in seq_along(beta.fp)){
		bdiv[[beta.fp[i]]] <- read.table(beta.fp[i],sep='\t',head=T,row=1,check=F)
		bdiv[[beta.fp[i]]]  <- bdiv[[beta.fp[i]]][ix,ix]
		bdiv[[beta.fp[i]]] <- cmdscale(bdiv[[beta.fp[i]]],3)
	}
	adiv <- adiv[ix,]
		
	return( list(
		map=map,
		otus=otus,
		taxa=taxa,
		bdiv=bdiv,
		adiv=adiv
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