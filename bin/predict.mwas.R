# testing function
# 
# -----
# input:
#       x : original feature vector
#    model: trained model with selected feature vector index
#       y : desired response [optional]
# -----
#  output:
#   export training results as a file, including train model and model evaluation
# 
"predict.mwas" <- function(x, model, y){
  if ("feat.set" %in% model){
    test.set <- x[, model$feat.set]
  }
  else test.set <- x
  
  if (exists(y)) pred.obj <- model.evaluation.mwas(test.set, model$trained.model, y)
  else pred.obj <- model.evaluation.mwas(test.set, model$trained.model)
  
  export.mwas(test.results = pred.obj)
  return()
}