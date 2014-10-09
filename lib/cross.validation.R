"mwas.cross.validation" <- function(x, y, nfolds=10, classifier = c("RF","SVM", "knn")[1], ...){
  
  if(classifier == "RF")
    best.model <- best.randomForest(x, y, tunecontrol = tune.control(random=TRUE, sampling="cross", cross = nfolds, fix = 2/3))
  else if(classifier == "SVM")
    best.model <- best.svm(x, y, tunecontrol = tune.control(random=TRUE, sampling="cross", cross = nfolds, fix = 2/3))
  #else if(classifier == "MLR") 
  #  best.model <- best.glmnet(x, y, tunecontrol = tune.control(random=TRUE, sampling="cross", cross = nfolds, fix = 2/3))
  else if (classifier == "knn")
    best.model <- tune.knn(x, y, k =c(1,3,5,7,9,11,13,15,17,19,21,23,25), l=NULL, sampling="cross", cross = nfolds)  
  return(best.model)
}

