# testing function
# 
# -----
# input:
#       x : original feature vector
#    model: trained model
#       y : desired response [optional]
#  is.feat: whether feature selection is required. [default: TRUE] 
# -----
#  output:
#   export training results as a file, including train model and model evaluation
# 
"predict.mwas" <- function(x, model, y, ...){
  if (is.feat){
    feat.set <- feature.scores.mwas(x, y, selection_threshold = 0)
    train.set <- x[,feat.set$ix]
  }
  else train.set <- x
  
  best.model <- persist.model.mwas(train.set, y, nfolds=10, classifier=method, ...)
  
  export.mwas(trained.model = best.model, feat.set = feat.set)
  return()
}