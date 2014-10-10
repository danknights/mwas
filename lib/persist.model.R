# training function
# 
# -----
# input:
#       x : feature vector
#       y : desired response
#  is.feat: whether feature selection is required. [default: TRUE] 
#  method : classifier type
# 
#  output:
#   export training results as a file, including train model and model evaluation
# 
"train.mwas" <- function(x, y, is.feat = TRUE, 
                         method=c("RF","SVM", "knn", "MLR")[1], ...){
  if (is.feat){
    feat.set <- feature.scores.mwas(x, y, selection_threshold = 0)
    train.set <- x[,feat.set$ix]
  }
  else train.set <- x
  
  best.model <- persist.model.mwas(train.set, y, nfolds=10, classifier=method, ...)
  
  export.mwas(trained.model = best.model, feat.set = feat.set)
  return()
}


# nested cross-validation or Jackknifing
# -- input
#          x : feature vector
#          y : desired labels
#     nfolds : number of folds 
# classifier : type of classifier
#       opts : options passed from keyboards
# -- output
#  best.model : model parameters

"persist.model.mwas" <- function(x, y, nfolds=10, 
                               classifier=c("RF","SVM", "knn", "MLR")[1],...){
  # x - feature set (observation * features)
  # y - desried response
  cv.ind <- sample(dim(x)[1])
  cv.samp.num <- floor(length(cv.ind)/nfolds)
  sampl_ind <- seq(1, dim(x)[1], by=1)
  candidate.model <- list()
  candidate.rocobj <- list()
  
  # nested cross-validation
  for (cv.id in 1:nfolds){
    if(cv.id < nfolds)
      idx <- cv.ind[seq((cv.id-1)*cv.samp.num + 1, cv.id*cv.samp.num, by=1)]
    else
      idx <- cv.ind[seq((cv.id-1)*cv.samp.num + 1, length(cv.ind), by=1)]
    
    train.set <- x[idx,]
    train.labels <- y[idx]
    
    validation.set <- x[!sampl_ind %in% idx,]
    validation.labels <- y[!sampl_ind %in% idx]
    
    candidate.model[[cv.ind]] <- cross.validation(train.set, train.labels, nfolds, classifier, ...)
    candidate.rocobj[[cv.ind]] <- roc.mwas(validation.set, model = candidate.model[cv.ind], response = validation.labels)
  }
  
  ####### ISSUE: train final model on whole data set
  #best.ind <- which.max(candidate.rocobj$auc) #### find the best auc index?
  #best.model <- candidate.model[best.ind]
  best.model <- cross.validation(x, y, nfolds, classifier, ...)
  
  
  ####### ISSUE: Calculate mean and std of error/AUC, MCC, Kappa
  # using candidate.rocobj
  #
  
  
  
  # if (savefile) save(best.model, file = paste(opts$outdir,"/trained.model", collapse='', sep=''))
  # Saving file is integrated into export.mwas.R
  
  return(best.model)
}


# Jackknife 
# 1) parameters estimated from the whole sample data
# 2) each element is, in turn, dropped from the sample 
#    and the parameter of interest is estimated from this smaller sample.
# 3) the difference between the whole sample estimation and the partial estimate 
#    is computed --- called pseudo-values
# 4) The pseudo-values are used in lieu of the original values to estimate the 
#    parameter of interest and their standard deviation is used to estimate the 
#    parameter standard error.
# --------
# -- input
#          x : feature vector
#          y : desired labels
#     nfolds : number of folds 
# classifier : type of classifier
# -- output
#  model estimation : error and standard deviation

"jackknife.mwas" <- function(x, y, nfolds=10, 
                                classifier=c("RF","SVM", "knn", "MLR")[1],...){
  # x - feature set (observation * features)
  # y - desried response
  cv.ind <- sample(dim(x)[1])
  cv.samp.num <- floor(length(cv.ind)/nfolds)
  sampl_ind <- seq(1, dim(x)[1], by=1)
  candidate.model <- list()
  candidate.rocobj <- list()
  
  # jackknife function here
  #
  
  return(model.estimate)
}



