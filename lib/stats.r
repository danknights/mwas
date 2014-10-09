# returns p-values, q-values, 
# indices of those below alpha, and renamed features with '*' etc.
# and subset of data if requested
# add.stars adds stars to the column names that are significant at .05, .01, .001, .0001
# if filename is provided, saves result of statistical test to file
# if parametric, uses t-test/ANOVA; else uses mann-whitney U/Kruskal-wallis
"differentiation.test" <- function (x,category, alpha=0.05, parametric=TRUE,
		include.subset=FALSE){
	category <- as.factor(as.character(category))
	if(length(unique(category)) < 2) stop('Category only has one level')
	if(parametric){
		pvals <- apply(x,2, function(taxon){
				if(var(taxon) == 0){
					NA
				} else {
					summary(lm(taxon~category))[[4]][2,4]
				}
			})
		stats <- apply(x,2, function(taxon){
				if(var(taxon) == 0){
					NA
				} else {
					summary(lm(taxon~category))[[4]][2,1]
				}
			})
	} else {
		if(length(levels(category)) == 2){
			ix1 <- category == levels(category)[1]
			ix2 <- category == levels(category)[2]
			pvals <- apply(x,2,function(taxon) wilcox.test(taxon[ix1],taxon[ix2],exact=FALSE)$p.value)
			stats <- apply(x,2,function(taxon) wilcox.test(taxon[ix1],taxon[ix2],exact=FALSE)$statistic)
		} else {
			pvals <- apply(x,2,function(taxon) kruskal.test(taxon ~ category)$p.value)
			stats <- apply(x,2,function(taxon) kruskal.test(taxon ~ category)$statistic)
		}
	}
	na.ix <- is.na(pvals)
	
	adj.pvals <- rep(NA,length(pvals))
	names(adj.pvals) <- names(pvals)
	adj.pvals[!na.ix] <- p.adjust(pvals[!na.ix],'fdr')
	keep.ix <- adj.pvals < alpha
	keep.ix[is.na(keep.ix)] <- FALSE
	if(!any(keep.ix)) stop('select.features failed to find any features.')

	# add stars to column names based on significance
	annotations <- colnames(x)
	thresholds <- c(.05, .01, .001, .0001)
	for(i in seq_along(thresholds)){
		
		star.ix <- adj.pvals[!na.ix] <= thresholds[i]
		
		if(any(star.ix)){
			for(j in which(star.ix)){
				annotations[!na.ix][j] <- paste(annotations[!na.ix][j],'*',sep='')
			}
		}
	}
	result <- list()
	result$annotations <- annotations
	result$features <- which(keep.ix)
	result$qvalues <- adj.pvals
	result$pvalues <- pvals
	result$stats <- stats
	
	# classwise means
	result$classwise.means <- t(apply(x,2,function(xx) sapply(split(xx,category),mean)))
	colnames(result$classwise.means) <- sprintf('%s mean',colnames(result$classwise.means))
	
	if(include.subset){
		result$subset <- x[,keep.ix,drop=F]
		colnames(result$subset) <- annotations[keep.ix]
	}
	
	
	return(result)
}


# saves list of results from differentiation.test to file (or prints)
"write.differentiation.test.results" <- function(results, filename='differentiated.features.txt'){
	if(!is.null(filename)){
		scipen.save <- options('scipen')
		options(scipen=20)
		hits <- cbind(results$pvalues, results$qvalues)
		hits <- cbind(hits, results$classwise.means)
		colnames(hits)[1:2] <- c('pvalue','qvalue')
		hits <- hits[!is.na(hits[,1]),]
		hits <- hits[order(hits[,1]),]
		sink(filename)
		cat('Feature\t')
		write.table(hits,quote=F,sep='\t')
		sink(NULL)
		options(scipen=scipen.save)
	}
}

