# nested cross-validation or Jackknifing
# -- input
#          x : feature vector
#          y : desired labels
#     nfolds : number of folds 
# classifier : type of classifier
#   savefile : save model parameters
#       opts : options passed from keyboards
# -- output
#  best.model : model parameters

"test.persist.model" <- function(x, y, nfolds=10, 
                                 classifier=c("RF","SVM", "knn", "MLR")[1], 
                                 savefile = TRUE, opts, ...){
  #
  #
  #
  #
  
}