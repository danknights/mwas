# runs jackknife ROC and plots a smoothed ROC curve; returns AUC
# bootstrap.fraction = -1 means leave-one-out
# outcome must be T/F
"bootstrapped.ROC" <- function(predictors, outcome, nreps=10,
		bootstrap.fraction=.9,
		filename='ROC-jackknifed.pdf', do.plot=TRUE,
		title.text='ROC', include.xy=TRUE){
	require('ROCR')
	if(is.null(dim(predictors))) predictors <- matrix(predictors,ncol=1)
	N <- nrow(predictors)
	if(bootstrap.fraction == -1) {
		folds <- sapply(1:N,function(ixx) (1:N)[-ixx])
	} else {
		bootstrap.n <- max(min(round(bootstrap.fraction * N), N-1),2)
		folds <- replicate(nreps,sample(N,size=bootstrap.n))
	}
	tprs <- NULL
	fprs <- NULL
	aucs <- NULL	

	for(fold in 1:nreps){
		fold.ix <- folds[,fold]
		outcome.i <- outcome[fold.ix]
		predictors.i <- predictors[fold.ix,,drop=F]
		res <- logistic.ROC(predictors.i, outcome.i)

		tprs <- cbind(tprs, res$tprs)
		fprs <- cbind(fprs, res$fprs)
		aucs <- c(aucs, res$auc)
	}


	require('flux')
	fprs.mean=rowMeans(fprs)
	tprs.mean=rowMeans(tprs)
	auc.mean <- auc(fprs.mean,tprs.mean)

	res <- list(fprs.mean=fprs.mean, tprs.mean=tprs.mean,
				fprs=fprs, tprs=tprs,
				auc.mean=auc.mean, aucs=aucs,
				nreps=nreps,
				folds=folds)
	if(do.plot) plot.jackknifed.ROC(res, filename=filename, title.text=title.text, include.xy=include.xy)
	
	invisible(res)
}

"logistic.ROC" <- function(predictors, outcome){
	require('ROCR')
	require('flux')
	mylogit <- glm(outcome ~ predictors, family=binomial(link=logit))
	roc <- probabilities.to.ROC(fitted(mylogit), outcome)
	return(roc)
}

"plot.jackknifed.ROC" <- function(jackknifed.ROC.res, filename='ROC-jackknifed.pdf',
		title.text='ROC',include.xy=TRUE){
	require('ROCR')
	require('flux')
	require('RColorBrewer')
	fprs <- jackknifed.ROC.res$fprs
	tprs <- jackknifed.ROC.res$tprs
	fprs.mean <- jackknifed.ROC.res$fprs.mean
	tprs.mean <- jackknifed.ROC.res$tprs.mean
	auc.mean <- jackknifed.ROC.res$auc.mean
	nreps <- jackknifed.ROC.res$nreps
	
	# do plot
	if(!is.null(filename)) pdf(filename,width=4,height=4)
	plot(0,0,type='n',xlim=0:1,ylim=0:1,xlab='False positive rate',ylab='True positive rate',main=title.text)
	if(include.xy) abline(0,1,col='red',lty=2,lwd=.5)
	for(i in 1:nreps){
		lines(fprs[,i],tprs[,i],lty=3,col=sprintf('%s88',brewer.pal(9,'Set1')[2]))
	}
	lines(fprs.mean,tprs.mean,col='black',lwd=2)
	legend('bottomright',sprintf('Area under the curve = %.2f',auc.mean),
			pt.cex=.1,adj=.05,cex=.85)
	if(!is.null(filename)) dev.off()
}

# y must be T/F
# p must be vector or probabilities of T
"probabilities.to.ROC" <- function(p,y,step.size=0.01){
	require('flux')
    fprs <- sapply(seq(1,0,-abs(step.size)), function(xx) mean(p[!y] >= xx))
    tprs <- sapply(seq(1,0,-abs(step.size)), function(xx) mean(p[y] >= xx))
    tprs[fprs > tprs] <- fprs[fprs > tprs]
	auc.res <- auc(fprs, tprs)

	return(list(fprs=fprs, tprs=tprs, auc=auc.res))
}

# assuming T/F response y and numeric p, finds the threshold which maximizes OR
# When using threshold to classify samples.
# also returns confidence intervals
"maximum.odds.ratio" <- function(p,y,maximize.lower.bound=TRUE){
	# try every unique threshold in p
	thresholds <- sort(unique(p))
	
	# skip lowest and highest p
	thresholds <- thresholds[2:(length(thresholds)-1)]
	n <- length(thresholds)
	lower.CI <- numeric(n)
	upper.CI <- numeric(n)
	expectation <- numeric(n)
	
	for(i in 1:length(thresholds)){
		res.i <- odds.ratio(p, y, thresholds[i])
		lower.CI[i] <- res.i$lower.CI
		upper.CI[i] <- res.i$upper.CI
		expectation[i] <- res.i$odds.ratio
	}
	
	if(maximize.lower.bound){
		best.ix <- which.max(lower.CI)
	} else {
		best.ix <- which.max(expectation)
	}
	res <- list(odds.ratio=expectation[best.ix],
				lower.CI=lower.CI[best.ix],
				upper.CI=upper.CI[best.ix],
				threshold=thresholds[best.ix])
	return(res)
}

# returns the odds ratio of y=TRUE in group p >= threshold compared to p < threshold
# p is probabilities of y=T
# y is T/F
"odds.ratio" <- function(p, y, threshold){
	tab <- table(y, p >= threshold)
	tab[tab == 0] <- .5
	print(tab)
	p1 <- tab[2,2] / colSums(tab)[2]
	p2 <- tab[2,1] / colSums(tab)[1]
	se <- sqrt(sum(1/as.numeric(tab)))
	expectation <- p1 / p2
	lower.CI <- exp(log(expectation) - se * 1.96)
	upper.CI <- exp(log(expectation) + se * 1.96)

	
	res <- list(odds.ratio=expectation,
				lower.CI=lower.CI,
				upper.CI=upper.CI)
	return(res)
}

# Derives an additive risk index: bad bugs minus good bugs
# for y == TRUE
# x is taxon/OTU/feature table
# y is TRUE/FALSE for disease or other phenotype
# alpha determines FDR threshold for inclusion of bugs
# eps is small value to add to all zero entries in feature table
# 
# return value includes risk.index.f, a function that takes a matrix with the same
# named columns as x and returns the risk index for each row.
"get.risk.index" <- function(x, y, alpha=0.25, eps=NULL,
		transform.type=c('none','asin-sqrt','sqrt','1.5-root','3-root','10-root','100-root')[1],
		parametric=FALSE, verbose=FALSE){

	x <- data.transform(x,transform.type)
	diff.tests <- differentiation.test(x, y, parametric=parametric)
	hit.ix <- which(diff.tests$qvalues <= alpha)

	if(length(hit.ix) == 0){
		stop('No hits found.\n')
	} else if (verbose) {
		cat("There were",length(hit.ix),"taxa significant at FDR of",alpha,'\n')
	}

	# get rid of zeros (if eps > 0)
	if(!is.null(eps)) x[x==0] <- eps

	# identify good/bad bugs
	good.bugs <- diff.tests$classwise.means[,1] > diff.tests$classwise.means[,2] & diff.tests$qvalues < alpha
	bad.bugs <- diff.tests$classwise.means[,1] < diff.tests$classwise.means[,2] & diff.tests$qvalues < alpha

	sum.good.bugs <- numeric(nrow(x))
	sum.bad.bugs <- numeric(nrow(x))
	if(sum(good.bugs) > 0) sum.good.bugs <- rowSums(x[,good.bugs,drop=F])
	if(sum(bad.bugs) > 0) sum.bad.bugs <- rowSums(x[,bad.bugs,drop=F])

	# risk index is simply additive bad bugs vs good bugs
	risk.index <- sum.bad.bugs - sum.good.bugs

	# risk.index.f 
	"risk.index.f" <- function(x){
		sum.good.bugs.i <- numeric(nrow(x))
		sum.bad.bugs.i <- numeric(nrow(x))
		if(sum(good.bugs) > 0) sum.good.bugs.i <- rowSums(x[,good.bugs,drop=F])
		if(sum(bad.bugs) > 0) sum.bad.bugs.i <- rowSums(x[,bad.bugs,drop=F])
		return(sum.bad.bugs.i - sum.good.bugs.i)
	}
	
	return(list(
		risk.index=risk.index,
		risk.index.f=risk.index.f,
		diff.tests=diff.tests,
		good.bugs=good.bugs,
		bad.bugs=bad.bugs,
		hit.ix=hit.ix,
		alpha=alpha
	))
}



# performs cross-validated assessment of risk index performance
# returns AUC, maximum f1, fprs/tprs for a sequence of thresholds
#
# nfolds=-1 means leave-one-out
"cv.risk.index" <- function(x, y, alpha=.25, nfolds=10, 
		transform.type='none'){
	x <- data.transform(x, transform.type)
	folds <- balanced.folds(y,nfolds=nfolds)
	
	risk.index <- numeric(nrow(x))
	for(fold in 1:length(unique(folds))){
		fold.ix <- which(fold == folds)
		res <- get.risk.index(x[-fold.ix,,drop=F],y[-fold.ix], alpha=alpha)
		risk.index[fold.ix] <- res$risk.index.f(x[fold.ix,,drop=F])
	}
	return(risk.index)
}