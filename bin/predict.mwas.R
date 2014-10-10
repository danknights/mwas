# testing function
# 
# -----
# input:
#       x : original feature vector
#    model: trained model with selected feature vector index
#       y : desired response [optional]
#  is.feat: whether feature selection is required. [default: TRUE] 
# -----
#  output:
#   export training results as a file, including train model and model evaluation
# 
"predict.mwas" <- function(x, model, y, ...){
  if ("feat.set" %in% model){
    test.set <- x[, model$feat.set]
  }
  else test.set <- x
  
  best.model <- model.evaluation.mwas(test.set, model$trained.model, nfolds=10, classifier=method, ...)
  
  export.mwas(trained.model = best.model, feat.set = feat.set)
  return()
}