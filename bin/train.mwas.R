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
  
  export.mwas(trained.model = best.model, feat = feat.set)
  return()
}