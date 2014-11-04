# Statistical testing function.
# returns p-values, q-values, 
# indices of those below alpha, and renamed features with '*' etc.
# and subset of data if requested
# add.stars adds stars to the column names that are significant at .05, .01, .001, .0001
# if filename is provided, saves result of statistical test to file
# if parametric, uses t-test/ANOVA; else uses mann-whitney U/Kruskal-wallis
# -------------
# Contributors: Dan
# -------------
# Input:
#     
# ------
# Output:
#   
# ------
#  Last update: 10/28/2014
"model.statistical.test.mwas" <- function(data, ...){
  
  options <- list(...)
  
  if (class(data) == "mwas"){
    x <- data$x 
    pc <- data$pc 
    fp <- data$fp 
    m <- data$m 
    taxon.names <- data$taxon.names
    category <- data$category
    is.multiple_axes <- data$is.multiple_axes
    category_order <- data$category_order
    is.sort_by_abundance <- data$is.sort_by_abundance
    num_taxa <- data$num_taxa
    alpha <- data$alpha
    x_axis_label <- data$x_axis_label
    out.dir <- data$outdir
    plot.type <- data$plottype
  } else{
    x <- data 
    pc <- options$pc 
    fp <- options$fp 
    m <- options$m 
    taxon.names <- options$taxon.names
    category <- options$category
    is.multiple_axes <- options$is.multiple_axes
    category_order <- options$category_order
    is.sort_by_abundance <- options$is.sort_by_abundance
    num_taxa <- options$num_taxa
    alpha <- options$alpha
    x_axis_label <- options$x_axis_label
    out.dir <- options$outdir
    plot.type <- options$plottype
  }
  switch(plot.type,
         beeswarm = {
           plot.differentiated.taxa(x=x, m=m, category=category, 
                                    num_taxa=num_taxa, 
                                    alpha=alpha, out.dir=out.dir, 
                                    is.sort_by_abundance=is.sort_by_abundance, 
                                    taxon.names=taxon.names,
                                    category_order=category_order, 
                                    x_axis_label=x_axis_label)
         },
         gradients = {
           processed.data <- preprocess.mwas(data)
           plot.gradients(x=processed.data$otu, 
                          pc=pc, fp=fp, m=m,
                          taxon.names=taxon.names, 
                          category=category,
                          is.multiple_axes=is.multiple_axes)
         },
         heatmap = { # need to fix
           heatmap.mwas(x, map, diff.features, cluster.var=c("Sex", "Treatment"), 
                        color.var=names(color.list), color.list, 
                        kegg_pathways=kegg_pathways, 
                        heatmap.title=heatmap.title, 
                        outputfile=fp)
         }, 
         scatter = {
           
         },
         stop("Please assign the correct plot type!(Optioins: beeswarm, graidents, heatmap, scatter.")
  )
}


"differentiation.test" <- function (x, category, alpha=0.05, parametric=FALSE, include.subset=FALSE){
	# category here is a column, not a name.
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
