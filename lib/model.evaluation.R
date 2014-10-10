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
  
