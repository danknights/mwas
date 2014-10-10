"model.evaluation.mwas" <- function(x, model, desired){
  # If the desired response is known, then evaluate the model, 
  # otherwise output the predicted labels.
  # -------------------
  # Model evaluation with different criteria. 
  # 
  # 1. classification accuracy
  # 2. area under the ROC (AUC)
  # 3. Matthew's correlation coefficients (MCC)
  # 4. Cohen's Kappa (not yet implemented)
  #
  # --- input: 
  #         x: feature vector
  #     model: trained model
  #   desired: desired output (response)
  #
  # --- output:
  #   evalobj: evaluation object - $error, $accuracy, $auc, $mcc, $kappa
  #

  predicted <- predict(model, x, decision.values = TRUE, probablity=TRUE)  # predicted output
  
  evalobj <- list()
  evalobj$prediction <- predicted
  
  if exists("desired"){ 
    # if desired response is known, then evaluate the classifier model
    # else output the predicted labels
    
    sample.num <- length(desired)

    # confusion matrix
    c.matrix <- t(sapply(levels(desired), function(level) table(predicted[desired==level])))
    
    rocobj <- roc.mwas(x, predicted, desired)
    
    evalobj$error <- sum(as.numeric(predicted != desired))/sample.num
    evalobj$acc <- (c.matrix[1,1]+c.matrix[2,2])/sum(c.matrix) # == 1 - evalobj$error 
    evalobj$auc <- rocobj$auc
    
    # MCC =  Pearson's correlation of y, yhat
    # MCC and Kappa for binary classification only
    if(is.binary){
      
    # MCC =  (TP*TN - FP*FN)/(sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN)))
    evalobj$mcc <- (c.matrix[1,1]*c.matrix[2,2])/(
      sqrt((c.matrix[1,1]+c.matrix[1,2])*
             (c.matrix[1,1]*c.matrix[2,1])*
             (c.matrix[2,2]*c.matrix[1,2])*
             (c.matrix[2,2]*c.matrix[2,1])))
    
    # Kappa = (Pr(a)-Pr(e))/(1-Pr(e))
    # Pr(a): relative observed agreement among raters = (TP+TN)/(P+N) = acc
    # Pr(e): the hypothetical probability of chance agreement 
    #        = ((TP+FP)/(P+N)*(TP+FP)/(P+N)+(FN+TN)/(P+N)*(FP+TN)/(P+N))
    Pr.e <- (c.matrix[1,1]+c.matrix[1,2])/sum(c.matrix)*(c.matrix[1,1]*c.matrix[2,1])/sum(c.matrix)+
      (c.matrix[2,1]*c.matrix[2,2])/sum(c.matrix)*(c.matrix[1,2]*c.matrix[2,2])/sum(c.matrix)
    evalobj$kappa <- (evalobj$acc - Pr.e)/(1 - Pr.e)
    }
    else {
      evalobj$mcc <- NULL
      evalobj$kappa <- NULL
    }
  }
      
  return(evalobj)
}
