# training function
# Contributors: Hu, Dan
# -----
# input:
#       x : feature vector
#       y : desired response - should be a factor for classification problem
#  is.feat: whether feature selection is required. [default: TRUE] 
#  method : classifier type
# 
#  output:
#   export training results as a file, including train model and model evaluation
# 
# ------
#  Last update: 10/25/2014
#
"train.mwas" <- function(x, y, is.feat = TRUE, 
                         method=c("RF","SVM", "knn", "MLR")[1], valid_type = c("kfold", "jackknife")[1]){
  
  require(e1071, quietly=TRUE, warn.conflicts=FALSE) 
  require(glmnet, quietly=TRUE, warn.conflicts=FALSE)
  
  if (is.feat){
    feat.set <- feature.scores.mwas(x, y, selection_threshold = 0)
    train.set <- x[,feat.set$ix]
  }
  else train.set <- x
  
  best.model <- persist.model.mwas(train.set, y, nfolds=10, classifier=method, valid_type)
  
  #export.mwas(trained.model = best.model, feat.set = feat.set)
  return(best.model)
}


# Nested cross-validation and Jack Knife validation for model performance estimation
# And model selection (model training)

"persist.model.mwas" <- function(x, y, nfolds=10, 
                                 classifier=c("RF","SVM", "knn", "MLR")[1],
                                 valid_type=c("kfold", "jackknife")[1]){
  
  require(pROC, quietly=TRUE, warn.conflicts=FALSE)
  
  # x - feature set (observation * features)
  # y - desried response
  cv.ind <- sample(dim(x)[1])   # permutate the index
  cv.samp.num <- floor(length(cv.ind)/nfolds)
  sampl_ind <- seq(1, dim(x)[1], by=1)
  candidate.model <- list()
  candidate.rocobj <- list()
  candidate.auc <- vector("numeric", nfolds)
  candidate.error <- vector("numeric", nfolds)
  
  # nested cross-validation
  if (tolower(valid_type)=="kfold"){ # case insensitive
    for (cv.id in 1:nfolds){
      # fold index that is being hold out
      if(cv.id < nfolds)
        idx <- cv.ind[seq((cv.id-1)*cv.samp.num + 1, cv.id*cv.samp.num, by=1)]
      else # the last fold could contain less than cv.samp.num of observations
        idx <- cv.ind[seq((cv.id-1)*cv.samp.num + 1, length(cv.ind), by=1)] 
      
      validation.set <- x[idx,]
      validation.labels <- y[idx]
      
      train.set <- x[!sampl_ind %in% idx,]
      train.labels <- y[!sampl_ind %in% idx]
      
      candidate <- cross.validation.mwas(train.set, train.labels, nfolds, classifier)
      candidate.eval <- model.evaluation.mwas(validation.set, candidate$best.model, validation.labels)
      
      candidate.model <- c(candidate.model, list(candidate))
      candidate.error[cv.id] <- candidate.eval$error
      candidate.auc[cv.id] <- candidate.eval$auc
      
    }
  }else if(valid_type=="jackknife") { # jackknife function 
    jk.ind <- seq(1, dim(x)[1], by=1)
    jk.samp.num <- floor(length(cv.ind)/nfolds)
    
    candidate.model <- list()
    candidate.rocobj <- list()
    candidate.auc <- vector("numeric", nfolds)
    candidate.error <- vector("numeric", nfolds)
    
    for (jk.fold in 1:nfolds){
      # fold index that is being hold out
      if(jk.fold < nfolds)
        idx <- jk.ind[seq((jk.fold-1)*jk.samp.num + 1, jk.fold*jk.samp.num, by=1)]
      else # the last fold could contain less than jk.samp.num
        idx <- jk.ind[seq((jk.fold-1)*jk.samp.num + 1, length(jk.fold), by=1)]
      
      train.set <- x[!sampl_ind %in% idx,]
      train.labels <- y[!sampl_ind %in% idx]
      
      candidate <- cross.validation.mwas(train.set, train.labels, nfolds, classifier)
      candidate.model <- c(candidate.model, list(candidate))
      candidate.error[jk.fold] <- candidate$best.performance
      
      candidate.obj <- roc.mwas(validation.set, model=candidate$best.model, response=validation.labels)
      candidate.rocobj <- c(candidate.model, list(candidate.obj))
      candidate.auc[jk.fold] <- as.numeric(candidate.rocobj$auc)
    }
  }
  
  ####### Calculate mean and std of error - model performance estimation
  # Explains how stable the model is.
  # 
  model.perform <- list(mean.error=mean(candidate.error), std.error=sd(candidate.error),
                        mean.auc=mean(candidate.auc),std.auc=sd(candidate.auc))
  
  ####### Train final model on whole data set
  #
  best.model.obj <- cross.validation.mwas(x, y, nfolds, classifier)
  best.model <- best.model.obj$best.model
  best.model.eval <- model.evaluation.mwas(x, best.model, y)

  # if (savefile) save(best.model, file = paste(opts$outdir,"/trained.model", collapse='', sep=''))
  # Saving file is integrated into export.mwas.R
  
  export.mwas(trained.model=best.model, trained.model.eval=best.model.eval, model.perform=model.perform)
  
  return(best.model)
}

# cross-validation 
# -- input
#          x : feature vector
#          y : desired labels
#     nfolds : number of folds 
# classifier : type of classifier
# -- output
#  best.model : model parameters

"cross.validation.mwas" <- function(x, y, nfolds=10, classifier = c("RF","SVM", "knn", "MLR")[1], ...){
  
  classifier <- tolower(classifier) # case insensitive for options
  if(classifier == "rf")
    best.model <- tune.randomForest(x, y, 
                                    tunecontrol = tune.control(random=TRUE, sampling="cross", 
                                                               cross = nfolds))
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
  
  else if(classifier == "svm") {
    best.model <- tune.svm(x, y, tunecontrol = tune.control(random=TRUE, sampling="cross", cross = nfolds))
    # $best.model         the model trained on the complete training data using the best parameter combination.
    # $best.parameters    a 1 x k data frame, k number of parameters
    # $best.performance   best achieved performance
    # $performances       a data frame with all parametere combinations along with the corresponding performance results
    # $train.ind          list of index vectors used for splits into training and validation sets
    
  }else if(classifier == "mlr") {
    y <- as.numeric(y)  # for regression, response should be numeric
    best.model <- cv.glmnet(x, y, nfolds = nfolds, family="multinomial")
    # $glmnet.fit   a fitted glmnet object for the full data
    # $cvm          mean cross-validation error    
    # $cvsd         estimate of standard error of cvm
    # $cvup         upper curve = cvm+cvsd
    # $cvlo         lower curve = cvm-cvsd
    # $nzero        number of non-zero coefficients at each lambda
    # $name         type of measure (for plotting purposes)
    # $lambda.min   value of labda that gives minimum cvm
    # $lambda.1se   largest value of lambda such that error is within 1 standard error of the minimum
    # $fit.preval, $foldid
    
  }else if (classifier == "knn"){
    best.model <- tune.knn(x, y, k = seq(from=1, to=min(dim(x)[1]/nfolds*(nfolds-1), 25), by =2), 
                           l= 1, sampling="cross", cross = nfolds, prob = TRUE) 
    class(best.model$best.model) <- 'knn'
    # $best.parameters$k          best k
    # $best.performance           classification error using best k
    # $performances$k             all candidate k's that used in cross-validation
    # $performances$error         errors for all candidates
    # $performances$dispersion    dispersions for all candidates
    # $best.model                 best model list 
  }
  else stop("Undefined classification method!")
  
  return(best.model)
}