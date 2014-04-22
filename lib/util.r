# loads experimental data into a list
# required is map
# otu table and taxa.fp are optional
# min.presence is the fraction of samples in which a feature must be present to keep
"load.experiment" <- function(map.fp, otu.fp=NULL, taxa.fp=NULL,
		alpha.fp=NULL, beta.fp=NULL, min.presence=0.1){

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