

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

# transforms data according to various transforms
"data.transform" <- function(x,transform_type='none'){
	if(transform_type == 'asin-sqrt'){
		x <- asin(sqrt(x))
	} else if(transform_type == 'norm-asin-sqrt'){
		x <- asin(sqrt(x))/asin(sqrt(1))
	} else if(transform_type == 'sqrt'){
		x <- x**(1/2)
	} else if(transform_type == '1.5root'){
		x <- x**(1/1.5)
	} else if(transform_type == '3root'){
		x <- x**(1/3)
	} else if(transform_type == '4root'){
		x <- x**(1/4)
	} else if(transform_type == '5root'){
		x <- x**(1/5)
	} else if(transform_type == '6root'){
		x <- x**(1/6)
	} else if(transform_type == '7root'){
		x <- x**(1/7)
	} else if(transform_type == '8root'){
		x <- x**(1/8)
	} else if(transform_type == '9root'){
		x <- x**(1/9)
	} else if(transform_type == '10root'){
		x <- x**(1/10)
	} else if(transform_type == '100root'){
		x <- x**(1/100)
	} else if(transform_type == '1000root'){
		x <- x**(1/1000)
	} else if(transform_type == '10000root'){
		x <- x**(1/10000)
	} else if(transform_type == '100000root'){
		x <- x**(1/100000)
	} else if(transform_type == '1000000root'){
		x <- x**(1/1000000)
	} else if(transform_type != 'none'){
		stop(paste('Unrecognized data transform type:',transform_type))
	}
	return(x)
}

