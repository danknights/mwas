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

#if (!require("e1071")) {
#  install.packages("e1071", dependencies = TRUE)
#  library(e1071)
#}
if (!require("kernlab")) {
  install.packages("kernlab", dependencies = TRUE)
  library(kernlab)
}
if (!require("glmnet")) {
  install.packages("glmnet", dependencies = TRUE)
  library(glmnet)
}
if (!require("randomForest")) {
  install.packages("randomForest", dependencies = TRUE)
  library(randomForest)
}
#require(e1071, quietly=TRUE, warn.conflicts=FALSE) 
#require(glmnet, quietly=TRUE, warn.conflicts=FALSE)
#require(randomForest, quietly=TRUE, warn.conflicts=FALSE)
#require(pROC, quietly=TRUE, warn.conflicts=FALSE)

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

# S3 method of predict.knn
"predict.knn" <- function(model, test){
  require(class)
  prediction <- knn(model$train, test, model$cl, k=model$k, l=model$l, prob=TRUE, use.all=TRUE)
  return(prediction)
}

bray_dist <- function(data){
  b_dist <- vegdist(data, method = "bray")
  b_dist <- as.matrix(b_dist)
  return(b_dist)
}

knn_dist <- function(data1, data2, y1, k){
  Ind1 <- dim(data1)[1]
  new_data <- rbind(data1, data2)
  dist_mat <-as.matrix(vegdist(new_data, method="bray"))
  Num <- dim(dist_mat)[1]
  y2 <- vector()
  for(id in (Ind1+1):Num){
    tt <- table(y1[match(names(sort(dist_mat[id, 1:Ind1])[1:k]), rownames(data1))])
    y2[id - Ind1] <- names(sort(tt, decreasing=T)[1])
  }
  return(y2)
}