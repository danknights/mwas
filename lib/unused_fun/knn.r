# performs k-nearest-neighbor classification 
# using a given distance matrix
# output is prediction vector
# if train.ix is NULL, uses leave-one-out predictions
# else uses only training indices
"knn.dist" <- function(d,y,train.ix=NULL,k=1){
	if(is.null(train.ix)){
		train.ix <- 1:length(y)
		test.ix <- 1:length(y)
		yhat <- y
	} else {
		test.ix <- (1:length(y))[-train.ix]
		yhat <- rep('',length(test.ix))
	}
	
	probs <- matrix(0,nr=length(test.ix),nc=length(unique(y)))
	rownames(probs) <- rownames(d)[test.ix]
	colnames(probs) <- sort(unique(y))	
	names(yhat) <- rownames(d)[test.ix]
	# one sample at a time
	for(i in seq_along(test.ix)){
		test.ix.i <- test.ix[i]
		train.ix.i <- setdiff(train.ix,test.ix.i)
		neighbor.order <- order(d[test.ix.i,train.ix.i])
		counts <- table(y[train.ix.i][neighbor.order[1:k]])
		# if tie, break randomly
		probs[i,names(counts)] <- counts / sum(counts)
		if(sum(counts == max(counts)) > 1){
			hits <- which(counts == max(counts))
			class.name <- names(counts)[sample(hits,1)]
		} else {
			class.name <- names(counts)[which.max(counts)]
		}
		yhat[i] <- class.name
	}

	return(list(y=y,yhat=yhat,probabilities=probs))
}


# provides cross-validated errors and predictions
# for knn
"cv.knn.dist" <- function(d,y,k=1,nfolds=10,verbose=FALSE){
    if(nfolds==-1) nfolds <- length(y)
	folds <- balanced.folds(y,nfolds=nfolds)
	result <- list()
    result$y <- as.factor(y)
    result$probabilities <- matrix(0,nr=length(y),nc=length(unique(y)))
    colnames(result$probabilities) <- sort(unique(y))
    result$predicted <- result$y
    result$errs <- numeric(length(unique(folds)))

    # K-fold cross-validation
    for(fold in sort(unique(folds))){
        if(verbose) cat(sprintf('Fold %d...\n',fold))
        foldix <- folds == fold
        knn.obj <- knn.dist(d,y,train.ix=(1:length(y))[!foldix],k=k)
		result$probabilities[foldix,] <- knn.obj$probabilities
		result$predicted[foldix] <- knn.obj$yhat
        result$errs[fold] <- mean(result$predicted[foldix] != result$y[foldix])
    }
    result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
    result$err <- mean(result$predicted != result$y)
    return(result)
}