# Generic method to evaluates the performance of an MWAS model
# Currently supported models:
# Random forests, SVM (rbf, linear, unifrac)
# KNN (unifrac, Bray-Curtis), logistic regression (penalized)
#
# Inputs:
# x (REQUIRED): data matrix (samples x features)
# y (REQUIRED): response variable (if character vector or factor, does classification)
#    (if numeric, does regression)
# model.type: string specifying the type of model: can be rf, svm, knn, glmnet
# distmat: distance matrix (for svm or knn) or R dist object
# distance.metric: string specifying the distance metric: can be bc (for Bray-Curtis), 
#     or jsd (for Jensen-Shannon)
# evaluation.type: string, can be "cv", "jackknife", "bootstrap", "built-in"
#	note: "jackknife" and "built-in" only work when a model has implicit estimation
#	of generalization error, such as random forests
# nfolds: number of folds for cross-validation or jackknifing
#	used only if evaluation.type == 'cv' or 'jackknife'
# tuning.parameters: a list of parameters for use in tuning of models
#    these are specific for each model
# 
# Outputs (as named list of):
# final.model: the tuned model, will work with generic S3 predict.mwas()
# predicted: the predicted values (based on cv or other hold-out predictions)
# evaluation.type: the type of evaluation performed
"evaluate.mwas" <- function(x, y,
		model.type=c('rf','svm','knn','glmnet')[1],
		d=NULL,
		distance.metric=c('bc','jsd')[1],
		evaluation.type=c('cv','jackknife','bootstrap','built-in')[1],
		nfolds=10,
		tuning.parameters=NULL)
){
	if(evaluation.type == 'cv'){
		evaluation <- evaluate.mwas.cv(x, y, ...)
	} else if (evaluation.type == 'bootstrap'){
		stop('Bootstrap evaluation not yet implemented\n')
	} else if (evaluation.type == 'jackknife'){
		stop('Jackknife evaluation not yet implemented\n')
	} else if (evaluation.type == 'built-in"){
		stop('Built-in evaluation not yet implemented\n')
	}
	
	evaluation$evaluation.type <- evaluation.type
	evaluation$nfolds <- nfolds
	evaluation$tuning.parameters <- tuning.parameters
	evaluation$distmat <- as.matrix(distmat)
	evaluation$distance.metric <- distance.metric
	return(evaluation)
}

# See evaluate.mwas for parameter descriptions
# output
# final.model: the tuned model, will work with generic S3 predict.mwas()
# predicted: the predicted values (based on cv or other hold-out predictions)
# probabilities: the predicted classwise probabilities
#                (based on cv or other hold-out predictions)
"evaluate.mwas.cv" <- function(x, y, nfolds=10,
			model.type,
			distmat=NULL,
			distance.metric=NULL){
	n <- nrow(x)
	n.classes <- length(unique(y))

	folds <- balanced.folds(y,nfolds)
	yhat <- y
	is.discrete <- !is.numeric(y)
	# probabilities are only used in discrete models
	prob <- matrix(0,nrow=n,ncol=n.classes)
	for(k in 1:nfolds){
		fold.ix <- folds == k
		mwas.model <- mwas.train(x[-fold.ix,,drop=F],y[-fold.ix],
										distmat=distmat,
										distance.metric=distance.metric,
										model.type=model.type,
										tuning.parameters=tuning.parameters
									  )
		prediction <- predict.mwas(mwas.model,x[fold.ix,,drop=F])
		yhat[fold.ix] <- prediction$predicted
				
		if(is.discrete){
			probabilities[fold.ix,,drop=F] <- prediction$probabilities
		}
	}
	if(is.discrete){
		evaluation <- evaluation.criteria.classification(y, yhat,
					probabilities=probabilities)
		} else {
			evaluation <- evaluation.criteria.regression(y, yhat)
		}
	}
	
	final.model <- mwas.train(x, y,
							distmat=distmat,
							distance.metric=distance.metric,
							model.type=model.type,
							tuning.parameters=tuning.parameters
					)
	evaluation$final.model <- final.model
	evaluation$predicted <- yhat
	evaluation$probabilities <- probabilities
	return(evaluation)
}

"evaluation.criteria.classification" <- function(y, yhat, probabilities=NULL){
  # Model evaluation with different criteria. 
  #
  # 1. classification accuracy
  # 2. area under the ROC (AUC)
  # 3. Matthew's correlation coefficients (MCC)
  # 4. Cohen's Kappa (not yet implemented)
  #
  # --- input: 
  #              y: true response vector
  #           yhat: predicted response vector
  #  probabilities: predicted probabilities matrix (samples x classes)
  #
  # --- output:
  #   evalobj: evaluation object - $error, $accuracy, $auc, $mcc, $kappa
  #
	y <- factor(y)
	yhat <- factor(yhat)
	is.binary <- length(levels(y)) == 2
	
    sample.num <- length(y)
    # ensure yhat has same levels as y
    levels(yhat) <- levels(y)
    
    # confusion matrix
    c.matrix <- table(y, yhat)
    
    # ROC stats
    if(!is.null(probabilities)){
	    rocobj <- roc.mwas(y, probabilities)
	}
    
    evalobj <- list()
    evalobj$error <- mean(y != yhat)
    evalobj$acc <- mean(y == yhat)
    evalobj$auc <- rocobj$auc
    
    # MCC =  Pearson's correlation of y, yhat
    # MCC and Kappa for binary classification only
    if(is.binary){
	    evalobj$mcc <- cor(y==levels(y)[1],yhat==levels(y)[1])

		# Kappa = (Pr(a)-Pr(e))/(1-Pr(e))
		# Pr(a): relative observed agreement among raters = (TP+TN)/(P+N) = acc
		# Pr(e): the hypothetical probability of chance agreement 
		#        = ((TP+FP)/(P+N)*(TP+FP)/(P+N)+(FN+TN)/(P+N)*(FP+TN)/(P+N))
		Pr.e <- (c.matrix[1,1]+c.matrix[1,2])/sum(c.matrix)*(c.matrix[1,1]*c.matrix[2,1])/sum(c.matrix)+
		  (c.matrix[2,1]*c.matrix[2,2])/sum(c.matrix)*(c.matrix[1,2]*c.matrix[2,2])/sum(c.matrix)
		evalobj$kappa <- (evalobj$acc - Pr.e)/(1 - Pr.e)
	} else {
		evalobj$mcc <- NULL
		evalobj$kappa <- NULL
	}
    return(evalobj)
}


"evaluation.criteria.regression" <- function(y, yhat){
  # Model evaluation with different criteria. 
  #
  # 1. Root mean squared error
  # 2. Pearson's correlation coefficient
  # 3. Spearman's correlation coefficient
  #
  # --- input: 
  #              y: true response vector
  #           yhat: predicted response vector
  #
  # --- output:
  #   evalobj: evaluation object - $rmse
  #
	  
}
  