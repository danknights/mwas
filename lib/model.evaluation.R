# If the desired response is known, then evaluate the model, 
# otherwise output the predicted labels.
# -------------------
# Contributors: Hu, Dan
# -------------------
# Model evaluation with different criteria. 
# 
# 1. classification accuracy
# 2. area under the ROC (AUC)
# 3. Matthew's correlation coefficients (MCC)
# 4. Cohen's Kappa
#
# --- input: 
#         x: feature vector
#     model: trained model
#   desired: desired output (response)
#
# --- output:
#   evalobj: evaluation object - $error, $accuracy, $auc, $mcc, $kappa, 
#                                $probabilities, $confusion
#
# --- 
# Last Update: 10/25/2014
#

"model.evaluation.mwas" <- function(data.set, model, desired){
  
  if (class(data.set)=="mwas"){
    x <- data.set$features
    model <- data.set$trained.model
    desired <- data.set$response
  }else x <- data.set
  
  evalobj <- list() # resutls object
  
  model.class <- tolower(class(model))
  #print(model.class)
  switch(model.class, 
         randomforest = {
           
           evalobj$probabilities <- predict(model, x, type="prob")
           evalobj$prediction <- predict(model, x, type="response")
           
         }, svm={
           predicted <- predict(model, x, decision.values = TRUE, probability=TRUE)
           
           evalobj$probabilities <- attr(predicted, "probabilities")
           evalobj$prediction <- predict(model, x)
           
         }, mlr={
           
           evalobj$prediction <- predict(model, x, type="link") # gives the linear predictors
           evalobj$probabilities <- predict(model, x, type="response") # gives the fitted probabilities
           
         }, knn={
           predicted <- predict(model, x)
           evalobj$prediction <- predicted
           evalobj$probabilities <- attr(predicted, "prob")
         }, stop('Please choose a classification method!')
  )
  
  #evalobj$prediction <- predicted
  
  is.binary <- length(levels(evalobj$prediction)) == 2
  if (!is.null(desired)){ 
    # if desired response is known, then evaluate the classifier model
    # else output the predicted labels
    
    sample.num <- length(desired)

    # confusion matrix
    c.matrix <- t(sapply(levels(desired), function(level) table(evalobj$prediction[desired==level])))

    rocobj <- roc.mwas(x, predicted = evalobj$prediction, response = desired)
    AUC <- rocobj$auc
    
    evalobj$confusion.matrix <- c.matrix
    error <- sum(as.numeric(evalobj$prediction) != as.numeric(desired))/sample.num
    accuracy <- sum(diag(c.matrix))/sum(c.matrix)  # == 1 - evalobj$error 
   
    
    # MCC =  Pearson's correlation of y, yhat
    # MCC and Kappa for binary classification only
    if(is.binary){
      
    # MCC =  (TP*TN - FP*FN)/(sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN)))
    MCC <- (c.matrix[1,1]*c.matrix[2,2])/(
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
    Kappa <- (accuracy - Pr.e)/(1 - Pr.e)
    } else {
      MCC <- NULL
      Kappa <- NULL
    }
  
    evalobj$performance = c(error, accuracy, AUC, MCC, Kappa)
    if (length(evalobj$performance)==3) { names(evalobj$performance) <- c("error", "accuracy", "AUC") 
    }else names(evalobj$performance) <- c("error", "accuracy", "AUC", " Matthews_corr_coeff","Cohens_Kappa")
  }
  return(evalobj)
}
