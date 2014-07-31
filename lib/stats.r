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


# performs a linear test of association of predictors in X (e.g. clinical data)
# with each column in Y (e.g. bacterial taxa).
# includes option for mixed linear model with a random effect for e.g. multiple samples
# in the same subject.
# rv.id is the ID of the random variable (subject ID column header is the intention)
# use drop.outlier.range=NA to not drop outliers
"association.tests.mixed" <- function(X,Y,random.id=NULL,mixed=FALSE,
									 remove.vowels.from.taxon.names=TRUE,
									 include.residuals=FALSE,
									 drop.outlier.range=3,
									 use.qvalue=FALSE,
									 verbose=FALSE,
									 test.type=c('linear','np')[1]
                            ){
    if(mixed && test.type == 'np') stop('Cannot perform mixed modeling with nonparametric test')
	rY <- Y # residuals
	rY[,] <- NA
	.X.72828 <<- as.data.frame(X)
	na.ix <- rowSums(is.na(.X.72828)) > 0
	.X.72828 <- droplevels(.X.72828[!na.ix,,drop=F])
	Y <- Y[!na.ix,]
    
	fixed.ids <- setdiff(colnames(X),random.id)
	# assumes dependent var will be "y"
	fixed.formula <- as.formula(paste("y.72828", paste(fixed.ids, collapse=" + "), sep=" ~ "))
	if(mixed){
		random.formula <- as.formula(sprintf('~ 1 | %s', random.id))
	}

	index <- 1

	# make one model matrix to determine number of covariates
	y.72828 <<- Y[,1]
	n.covar <- ncol(model.matrix(fixed.formula, data=.X.72828)) - 1
	res <- matrix(0,nrow=n.covar * ncol(Y), ncol=3); # will be numeric
	desc <- matrix("", nrow=nrow(res), ncol=2) # will be character
	colnames(res) <- c('pvalue','qvalue','coefficient')
	colnames(desc) <- c('Covariate','Taxon')
	rx <- Y # residuals
	rx[,] <- NA
    count <- 1
	for(i in 1:ncol(Y)){
		y.72828 <<- Y[,i]
		if(is.na(drop.outlier.range)){
			outlier.72828 <- rep(FALSE,length(y.72828))
		} else {
			outlier.72828 <<- is.outlier(y.72828,range=drop.outlier.range)
		}
		if(var(y.72828[!outlier.72828]) == 0){
			if(verbose) cat('No variance:',colnames(Y)[i],'\n')
		} else if(mean(outlier.72828) > .1){
			if(verbose) cat('Outliers:',mean(outlier.72828),'for',colnames(Y)[i],'\n')
		} else {
			if(mixed){
				m <- try(lme(fixed=fixed.formula,
						random=random.formula,
						data=.X.72828,
						subset=!outlier.72828))
				if(class(m)=='try-error'){
					res[ix,c('pvalue','coefficient')] <- c(NA,NA)
					rY[rownames(Y),i][!outlier.72828] <- NA
					cat('Error fitting',colnames(Y)[i],'\n')
				} else {
					res.i <- summary(m)[[20]]
					ix <- count:(count + nrow(res.i) - 2)
					res[ix,c('pvalue','coefficient')] <- res.i[-1,c('p-value','Value')]
					rY[rownames(Y),i][!outlier.72828] <- resid(m)
				}
			} else {
				m <- lm(fixed.formula, data=.X.72828, subset=!outlier.72828)
				res.i <- summary(m)[[4]]
				ix <- count:(count + nrow(res.i) - 2)
				res[ix,c('pvalue','coefficient')] <- res.i[-1,c('Pr(>|t|)','Estimate')]
				rY[rownames(Y),i][!outlier.72828] <- resid(m)
			}
			desc[ix,'Taxon'] <- colnames(Y)[i]
			desc[ix,'Covariate'] <- rownames(res.i)[-1]
			count <- max(ix) + 1
		}
	}
	# if res is not full, drop empty rows
	if(count <= nrow(res)){
	    res <- res[-(count:nrow(res)),]
	    desc <- desc[-(count:nrow(desc)),]
	}
	desc <- desc[!is.na(res[,'pvalue']),]
	res <- res[!is.na(res[,'pvalue']),]
    if(use.qvalue){
        require('qvalue')
        qvals <- try(qvalue(res[,'pvalue'])$qvalue,silent=TRUE)
    }
    if(!use.qvalue || class(qvals)=='try-error'){
		qvals <- p.adjust(res[,'pvalue'],'fdr')
    }
    res[,'qvalue'] <- qvals

	desc <- desc[order(res[,'pvalue']),]
	res <- res[order(res[,'pvalue']),]
	res.df <- as.data.frame(desc)
	res.df <- cbind(res.df,res)
	
	# make taxon names shorter
	taxon.names <- res.df$Taxon
	for(to.replace in c('k__',' p__',' c__',' o__',' f__',' g__',' s__')){
		taxon.names <- gsub(to.replace,'',taxon.names)
	}
	if(remove.vowels.from.taxon.names){
		taxon.names <- gsub('[aeiou]','',taxon.names)
	}
	res.df$Taxon_abbrev <- taxon.names

	if(include.residuals){
		return(list(res=res.df, residuals=rY))
	} else {
		return(res.df)
	}
}

"association.tests.mixed.with.partial.residuals" <-function(X,Y,random.id=NULL,mixed=FALSE,
									 remove.vowels.from.taxon.names=FALSE, drop.outlier.range=3,
									 use.qvalue=FALSE, verbose=FALSE){
	# drop constant columns
	X <- X[,apply(X,2,function(xx) length(unique(xx)) > 1)]

	full.res <- association.tests.mixed(X,Y,random.id=random.id,mixed=mixed,
					remove.vowels.from.taxon.names,include.residuals=TRUE, drop.outlier.range=drop.outlier.range, use.qvalue=use.qvalue)
	partial.residuals <- list()
	random.var.ix <- rep(FALSE,ncol(X))
	if(mixed) random.var.ix <- colnames(X) == random.id
	for(i in (1:ncol(X))[!random.var.ix]){
		if(verbose) cat('Getting partial residuals for covariate ',colnames(X)[i],'...\n',sep='')
		partial.residuals[[colnames(X)[i]]] <- 
			association.tests.mixed(X[,-i],Y,random.id=random.id,mixed=mixed,
				remove.vowels.from.taxon.names,include.residuals=TRUE,use.qvalue=use.qvalue, drop.outlier.range=drop.outlier.range)$residuals
	}
	return(list(res=full.res$res,
			residuals=full.res$residuals,
			partial.residuals=partial.residuals
	))
}

# tests whether each column of Y (e.g. taxa) can be predicted by each column of X
# (e.g. genetics), while controlling for columns of Z
"mwas.xwas.association.tests" <- function(X,Y,Z,
			drop.outlier.range=3,
			include.residuals=FALSE,
			verbose=FALSE
			){
	gx.res.list <- list()
	gx.stats <- NULL
	for(i in 1:ncol(X)){
		gx.id <- colnames(X)[i]
		if(verbose) cat(i,': ', gx.id,'\n',sep='');
		ix.i <- !is.na(X[,i])
		XX <- cbind(gx=X[ix.i,i],Z[ix.i,])
		XX <- droplevels(XX[,apply(XX,2,function(xx) length(unique(xx))>1)])
		colnames(XX)[1] <- gx.id
		res <- association.tests.mixed(XX,Y[ix.i,],drop.outlier.range=drop.outlier.range)
		if(include.residuals){
			res.resid <- association.tests.mixed(Y[,-1],asin(sqrt(x[ix.i,taxon.ix])),include.residuals=TRUE,drop.outlier.range=drop.outlier.range)
			gx.res.list[[gx.id]] <- res.resid$residuals
		}
		gx.stats <- rbind(gx.stats,res[grep(gx.id,res[,'Covariate']),])
	}

	gx.stats$qvalue <- p.adjust(gx.stats$pvalue,'fdr')
	gx.stats <- gx.stats[order(gx.stats$pvalue),]

	if(include.residuals){
		return(list(result=gx.stats, residuals=gx.res.list))
	} else {
		return(gx.stats)
	}
}



# tests whether each column of Y (e.g. taxa) can be predicted by each column of X
# (e.g. genetics), while controlling for columns of Z
# uses monte carlo-based test to compare coefficients to permuted data
"mwas.xwas.association.tests.mc" <- function(X,Y,Z,
			drop.outlier.range=3,
			verbose=FALSE,
			nperm=999
			){
	# get actual stats		
	gx.stats <- mwas.xwas.association.tests(X,Y,Z,drop.outlier.range=drop.outlier.range,
			include.residuals=FALSE,
			verbose=verbose)
	
	# get permuted stats, lots of times
	mc.gx.stats <- replicate(nperm,
			mwas.xwas.association.tests(X[sample(nrow(X)),,drop=F],Y,Z,drop.outlier.range=drop.outlier.range,
			include.residuals=FALSE,
			verbose=verbose)
		)

	
	# calculate p-value for each column of X for each column of Y
	if(verbose) cat('Tallying mc results...\n')
	for(i in 1:nrow(gx.stats)){
		obs.coef <- gx.stats[i,'coefficient']
		obs.taxon.name <- gx.stats[i,'Taxon']
		obs.cov.name <- gx.stats[i,'Covariate']
		mc.coef <- obs.coef
		for(j in 1:length(mc.gx.stats['coefficient',])){
			# get index of the coefficient for this observed Taxon/Covariate pair
			ix <- mc.gx.stats['Taxon',][[j]] == obs.taxon.name
			ix <- ix & mc.gx.stats['Covariate',][[j]] == obs.cov.name
			if(sum(ix) != 1) browser()
			mc.coef <- c(mc.coef, mc.gx.stats['coefficient',][[j]][ix])
		}
		gx.stats[i,'pvalue'] <- mean(abs(mc.coef) >= abs(obs.coef))		
	}
	return(gx.stats)
}




# runs cross-validation 
# if predict.fun is NULL, uses S3 predict method
# if nfolds > length(y) or nfolds==-1, uses leave-one-out cross-validation
# ...: additional parameters for train.fun
#
# value:
# y: true values
# predicted: cv predicted values
# probabilities: cv predicted class probabilities (or NULL if unavailable)
# confusion.matrix: confusion matrix (true x predicted)
# nfolds: nfolds
# params: list of additional parameters
# importances: importances of features as predictors
"rf.cross.validation.regression" <- function(x, y, nfolds=10, verbose=verbose, ...){
    if(nfolds==-1) nfolds <- length(y)
    folds <- sample(rep(1:nfolds,ceiling(length(y)/nfolds)))
    
    result <- list()
    result$y <- y
    result$predicted <- result$y
    result$importances <- matrix(0,nrow=ncol(x),ncol=nfolds)
    result$errs <- numeric(length(unique(folds)))

    # K-fold cross-validation
    for(fold in sort(unique(folds))){
        if(verbose) cat(sprintf('Fold %d...\n',fold))
        foldix <- which(folds==fold)
        model <- randomForest(x[-foldix,], factor(result$y[-foldix]), importance=TRUE, do.trace=verbose, ...)
        newx <- x[foldix,]
        if(length(foldix)==1) newx <- matrix(newx,nrow=1)
        result$predicted[foldix] <- predict(model, newx)
        probs <- predict(model, newx, type='prob')
        result$probabilities[foldix,colnames(probs)] <- probs
        result$errs[fold] <- mean(result$predicted[foldix] != result$y[foldix])
        result$importances[,fold] <- model$importance[,'MeanDecreaseAccuracy']
    }

	result$nfolds <- nfolds
    result$params <- list(...)
    result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
    return(result)    
}


# d is the data list
# d should contain at least:
# map: with 'PCDAI_06' column real values, to be predicted
# cvr: covariates table
# taxa: sample x taxon relative abundance matrix
# 
# if x is null, uses taxa object
"mc.test.taxa.prediction" <- function(d, y, x=NULL, nreps.obs=10, nreps.mc=10, nfolds=10, ntree=2000,
		covariate.names=NULL,
		subset.ix=NULL){ with(d, {
	# y <- rowMeans(map[,c('PCDAI_12','PCDAI_06')]) - map$PCDAI
	# y <- rowMeans(map[,c('PCDAI_12','PCDAI_06')])

	if(is.null(covariate.names)){
		cvr.i <- cvr
	} else {
		cvr.i <- map[,covariate.names]
	}
	
	if(is.null(subset.ix)){
		subset.ix <- 1:length(y)
	}
	
	if(!is.null(x)) taxa <- x
	obs.errs <- NULL
	probabilities <- matrix(NA,nrow=length(subset.ix),ncol=nreps.obs)
	rownames(probabilities) <- rownames(taxa)[subset.ix]
	
	for(i in 1:nreps.obs){
		cat('Obs iteration',i,'\n')
		y.i <- y
		x <- cbind(taxa, cvr.i[,colnames(cvr.i) != 'AnSubID'])
		y.i <- y.i[subset.ix]
		x <- x[subset.ix,]
		
		na.ix <- is.na(y.i) | rowSums(is.na(x)) > 0
		y.i <- y.i[!na.ix]
		x <- x[!na.ix,]

		ix <- x$Race != "Other"
		y.i <- factor(y.i[ix])
		x <- droplevels(x[ix,])

		rf.cv <- rf.cross.validation(x,y.i,
				nfolds=nfolds, verbose=FALSE,ntree=ntree,mtry=ncol(x)**.75)
		cat('OBS error:',mean(rf.cv$errs),'\n')
		probabilities[rownames(rf.cv$probabilities),i] <- rf.cv$probabilities[,2]
		if(i == 1) {
			y.to.return <- y.i
			names(y.to.return) <- rownames(x)
			predicted <- rf.cv$predicted
			names(predicted) <- rownames(x)
			importances <- rf.cv$importances
			rownames(importances) <- colnames(x)
			print(table(y.i,predicted))
			
		}
		
		obs.errs <- c(obs.errs, mean(rf.cv$errs))
	}
	cat('Error:',mean(obs.errs),'\n')
	probabilities <- rowMeans(probabilities)
	probabilities <- probabilities[!is.na(probabilities)]
	
	baseline.errs <- NULL
	mc.probabilities <- matrix(NA,nrow=length(subset.ix),ncol=nreps.mc)
	rownames(mc.probabilities) <- rownames(taxa)[subset.ix]

	if(nreps.mc > 0){
		for(i in 1:nreps.mc){
			cat('MC iteration',i,'\n')
			y.i <- y
			x <- cbind(taxa[sample(1:nrow(taxa)),], cvr.i[,colnames(cvr.i) != 'AnSubID'])
			y.i <- y.i[subset.ix]
			x <- x[subset.ix,]
		
			na.ix <- is.na(y.i) | rowSums(is.na(x)) > 0
			y.i <- y.i[!na.ix]
			x <- x[!na.ix,]

			ix <- x$Race != "Other"
	 		y.i <- factor(y.i[ix])
			x <- droplevels(x[ix,])
	
			rf.cv.mc <- rf.cross.validation(x,y.i, nfolds=nfolds, verbose=FALSE,ntree=ntree,mtry=ncol(x) ** .75)
			cat('MC error:',mean(rf.cv.mc$errs),'\n')
			baseline.errs <- c(baseline.errs, mean(rf.cv.mc$errs))
			mc.probabilities[rownames(cvr.i)[subset.ix][!na.ix][ix],i] <- rf.cv.mc$probabilities[,2]
		}
	}
	mc.probabilities <- rowMeans(mc.probabilities)
	mc.probabilities <- mc.probabilities[!is.na(mc.probabilities)]
	
	cat('\n')
	cat(sprintf('Observed error: %.3f +/- %.3f (%d reps)\n',
				mean(obs.errs),
				sd(obs.errs)/sqrt(nreps.obs),
				nreps.obs))
			
	cat(sprintf('MC error: %.3f +/- %.3f (%d reps)\n',
				mean(baseline.errs),
				sd(baseline.errs)/sqrt(nreps.mc),
				nreps.mc))
			
	pvalue <- (sum(baseline.errs <= mean(obs.errs)) / (1 + nreps.mc))
	cat(sprintf('MC p-value = %.5f\n', pvalue))

	return(list(obs.errs=obs.errs, mc.errs=baseline.errs,  pvalue=pvalue, predicted=predicted, y=y.to.return,
				importances=importances, probabilities=probabilities,
				mc.probabilities=mc.probabilities))
	
})}



# plots an ROC curve from matrix of probabilities
# only works for 2-way 
# y must be T/F
# p must be column predicting T or
# matrix where each column is a different receiver
"rf.ROC" <- function(p,y,step.size=0.01,filepath='ROC_curve.pdf',cols=NULL,do.plot=FALSE){
    if(is.null(cols)) cols <- c('black',brewer.pal(9,'Set1')[-6])
    if(is.null(dim(p))) p <- as.matrix(p,ncol=1)
    if(is.null(colnames(p))) colnames(p) <- sprintf('Classifier %d',1:ncol(p))
    fprs <- list(ncol(p))
    tprs <- list(ncol(p))
    for(i in 1:ncol(p)){
        fprs[[i]] <- sapply(seq(1,0,-abs(step.size)), function(xx) mean(p[!y,i] >= xx))
        tprs[[i]] <- sapply(seq(1,0,-abs(step.size)), function(xx) mean(p[y,i] >= xx))
        tprs[[i]][fprs[[i]] > tprs[[i]]] <- fprs[[i]][fprs[[i]] > tprs[[i]]]
    }

	if(do.plot){
		if(!is.null(filepath))	pdf(filepath,width=4,height=4)
		par(mar=c(5,4,3,1))
		plot(fprs[[1]],tprs[[1]],xlim=c(0,1),ylim=c(0,1),type='l',lwd=2,
				xlab='False positive rate',ylab='True positive rate',col=cols[1])
		abline(0,1,lty=2,col='#000000aa')
		grid()

		if(ncol(p) > 1){
			for(i in 2:ncol(p)){
				lines(fprs[[i]],tprs[[i]],xlim=c(0,1),ylim=c(0,1),type='l',lwd=2,col=cols[i-1])
			}
		}
		legend('topleft',colnames(p),lwd=2,col=cols,cex=.75)
		if(!is.null(filepath)) dev.off()
	}
	
	require('flux')
	for(i in 1:ncol(p)){
		
		auc.i <- auc(fprs[[i]],tprs[[i]])
		cat(sprintf('Column %d: %s\n',i,colnames(p)[i]))
		cat(sprintf('AUC = %f\n',auc.i))
	}
	# return last auc
	return(auc.i)
}


# use cv or jackknifing to estimate mean and confidence
# of rf.cross.validation error estimate
"rf.cross.validation.CI" <- function(x,y,nfolds=10,nfolds.inner=10,ntree=1000,
			verbose=FALSE){
    if(nfolds==-1) nfolds <- length(y)
    folds <- balanced.folds(y,nfolds=nfolds)
    errs <- numeric(length(unique(folds)))

    # K-fold cross-validation
    for(fold in sort(unique(folds))){
        if(verbose) cat(sprintf('Fold %d...\n',fold))
        foldix <- which(folds==fold)
		res <- rf.cross.validation(x[-foldix,],factor(y[-foldix]),nfolds=nfolds.inner,ntree=ntree)
		errs[fold] <- mean(res$err)
    }

    return(list(err=mean(errs),sd=sd(errs),errs=errs)) 
}

"rf.cross.validation.significance" <- function(x,y,ntree=1000,nfolds=10,nreps=100,
		confounder=NULL,
		confounder.control.type=c('linear','permute.within.levels','hold.out.by.level')[2],
		verbose=FALSE){
	y <- droplevels(as.factor(y))

	# remove linear confounder if requested
	if(!is.null(confounder) && confounder.control.type == 'linear'){
		x <- remove.linear.confounder(x,confounder)
	}

	# base err
	folds <- NULL
	base.rf <- rf.cross.validation(x, y, nfolds=nfolds, ntree=ntree)
	base.err <- mean(base.rf$err)
	
	# mc rf	
	mc.err <- numeric(nreps)
	n <- length(y)
	for(i in 1:nreps){
		if(verbose) cat('rep',i,'')
		if(!is.null(confounder)){
		 	if(confounder.control.type == 'permute.within.levels'){
				perm.ix <- permute.within.levels(confounder)
			} else if(confounder.control.type == 'hold.out.by.level'){
				perm.ix <- NULL
				# NOT FINISHED
			}
		} else {
			perm.ix <- sample(n)
		}
		rf <- rf.cross.validation(x, y[perm.ix],nfolds=nfolds,ntree=ntree)
		mc.err[i] <- mean(rf$err)
	}
	if(verbose) cat('\n')
	
	pval <- mean(mc.err <= base.err)
	if(verbose) cat('p-value:',pval,'\n')
	return(list(p.value=pval, base.err=base.err, mc.err=mc.err))
}



# runs cross-validation 
# if predict.fun is NULL, uses S3 predict method
# if nfolds > length(y) or nfolds==-1, uses leave-one-out cross-validation
# ...: additional parameters for train.fun
#
# value:
# y: true values
# predicted: cv predicted values
# probabilities: cv predicted class probabilities (or NULL if unavailable)
# confusion.matrix: confusion matrix (true x predicted)
# nfolds: nfolds
# params: list of additional parameters
# importances: importances of features as predictors
"rf.cross.validation" <- function(x, y, nfolds=10, folds=NULL, verbose=FALSE, ...){
	require('randomForest')
    if(nfolds==-1) nfolds <- length(y)
    if(is.null(folds)) folds <- balanced.folds(y,nfolds=nfolds)
    result <- list()
    result$y <- as.factor(y)
    result$predicted <- result$y
    result$probabilities <- matrix(0, nrow=length(result$y), ncol=length(levels(result$y)))
    rownames(result$probabilities) <- rownames(x)
    colnames(result$probabilities) <- levels(result$y)
    result$importances <- matrix(0,nrow=ncol(x),ncol=nfolds)
    result$errs <- numeric(length(unique(folds)))

    # K-fold cross-validation
    for(fold in sort(unique(folds))){
        if(verbose) cat(sprintf('Fold %d...\n',fold))
        foldix <- which(folds==fold)
        model <- randomForest(x[-foldix,], factor(result$y[-foldix]), importance=TRUE, do.trace=verbose, ...)
        newx <- x[foldix,]
        if(length(foldix)==1) newx <- matrix(newx,nrow=1)
        result$predicted[foldix] <- predict(model, newx)
        probs <- predict(model, newx, type='prob')
        result$probabilities[foldix,colnames(probs)] <- probs
        result$errs[fold] <- mean(result$predicted[foldix] != result$y[foldix])
        result$importances[,fold] <- model$importance[,'MeanDecreaseAccuracy']
    }
	rownames(result$importances) <- colnames(x)
	result$err <- mean(result$predicted != result$y)
    result$nfolds <- nfolds
    result$params <- list(...)
    result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
    return(result)    
}

