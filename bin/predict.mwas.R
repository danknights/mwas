# testing function
# Contributors: Hu
# -----
# input:
#       x : original feature vector
#    model: trained model with selected feature vector index - class(mwas)
#       y : desired response [optional]
# -----
#  output:
#   export training results as a file, including train model and model evaluation
# -----
# Last update: 10/25/2014
#

"predict.mwas" <- function(x, model, y){
  if ("feat.set" %in% model) test.set <- x[, model$feat.set] 
  else test.set <- x
  
  if (exists(y)) pred.obj <- model.evaluation.mwas(test.set, model$trained.model, y)
  else pred.obj <- model.evaluation.mwas(test.set, model$trained.model)
  
  export.mwas(trained.model.eval = pred.obj)
}