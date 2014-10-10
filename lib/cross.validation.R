# cross-validation 
# -- input
#          x : feature vector
#          y : desired labels
#     nfolds : number of folds 
# classifier : type of classifier
# -- output
#  best.model : model parameters

"cross.validation.mwas" <- function(x, y, nfolds=10, classifier = c("RF","SVM", "knn", "MLR")[1], ...){
  
  if(classifier == "RF")
    best.model <- tune.randomForest(x, y, tunecontrol = tune.control(random=TRUE, sampling="cross", cross = nfolds, fix = 2/3))
    # $best.performance     the error rate for the best model
    # $train.ind            list of index vectors used for splits into training and validation sets
    # $performances         error and dispersion; a data frame with all parametere combinations along with the corresponding performance results
    # $best.model           best model parameter list
    #   -- $predicted       the predicted values of the input data based on out-of-bag samples
    #   -- $err.rate        vector error rates of the prediction on the input data, the i-th element being the (OOB) error rate for all trees up to the i-th.
    #   -- $confusion       confusion matrix for the training set
    #   -- $votes           a matrix with one row for each input data point and one column for each class, giving the fraction or number of (OOB) ‘votes’ from the random forest
    #   -- $oob.times       number of times cases are ‘out-of-bag’ 
    #   -- $classes         class (category) names
    #   -- $importance      <check ?randomForest>
    #   -- $importanceSd    the “standard errors” of the permutation-based importance measure
    #   -- $localImportance <check ?randomForest>
    #   -- $ntree           number of trees grown
    #   -- $mtry            number of predictors sampled for spliting at each node
    #   -- $forest          a list that contains the entire forest
    #   -- $y               the desired labels for the training set
    #   -- $inbag, $test
  
  else if(classifier == "SVM")
    best.model <- tune.svm(x, y, tunecontrol = tune.control(random=TRUE, sampling="cross", cross = nfolds))
    # $best.model         the model trained on the complete training data using the best parameter combination.
    # $best.parameters    a 1 x k data frame, k number of parameters
    # $best.performance   best achieved performance
    # $performances       a data frame with all parametere combinations along with the corresponding performance results
    # $train.ind          list of index vectors used for splits into training and validation sets
      
  else if(classifier == "MLR") 
    best.model <- cv.glmnet(x, y, nfolds = nfolds)
    # $glmnet.fit   a fitted glmnet object for the full data
    # $cvm          mean cross-validation error    
    # $cvsd         estimate of standard error of cvm
    # $cvup         upper curve = cvm+cvsd
    # $cvlo         lower curve = cvm-cvsd
    # $nzero        number of non-zero coefficients at each lambda
    # $name         type of measure (for plotting purposes)
    # $lambda.min   value of labda that gives minimum cvm
    # $lambda.lse   largest value of lambda such that error is within 1 standard error of the minimum
    # $fit.preval, $foldid
  
  else if (classifier == "knn")
    best.model <- tune.knn(x, y, k =c(1,3,5,7,9,11,13,15,17,19,21,23,25), l=NULL, sampling="cross", cross = nfolds, prob = TRUE) 
    # $best.parameters$k          best k
    # $best.performance           classification error using best k
    # $performances$k             all candidate k's that used in cross-validation
    # $performances$error         errors for all candidates
    # $performances$dispersion    dispersions for all candidates
    # $best.model                 best model list 
  
  return(best.model)
}

