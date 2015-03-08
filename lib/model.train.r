# Training function
# Contributors: Hu, Dan
# -----
# input:
#    data : feature vector or a MWAS class object including all the parameters (the output of import.mwas)
#       y : desired response - should be a factor for classification problem
#  is.feat: whether feature selection is required. [default: TRUE] 
#  method : classifier type
#  valid_type : validation type - "kfold" or "jackknife"
# output:
#   export training results as a file, including train model and model evaluation
# 
# ------
#  Last update: 10/25/2014
#

if (!require("e1071", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("e1071", repos="http://cran.r-project.org", dependencies = TRUE)
  library("e1071", verbose=F, warn.conflicts =F)
}

if (!require("kernlab", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("kernlab", repos="http://cran.r-project.org", dependencies = TRUE)
  library("kernlab", verbose=F, warn.conflicts =F)
}
if (!require("glmnet", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("glmnet",repos="http://cran.r-project.org",  dependencies = TRUE)
  library("glmnet", verbose=F, warn.conflicts =F)
}
if (!require("randomForest", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("randomForest",repos="http://cran.r-project.org", dependencies = TRUE)
  library("randomForest", verbose=F, warn.conflicts =F)
}
if (!require("vegan", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("vegan", repos="http://cran.r-project.org", dependencies = TRUE)
  library("vegan", verbose=F, warn.conflicts =F)
}

if (!require("pROC", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("pROC", repos="http://cran.r-project.org", dependencies = TRUE)
  library("pROC", verbose=F, warn.conflicts =F)
}

#if (!require("pROC")) {
#  install.packages("pROC", dependencies = TRUE)
#  library(pROC)
#}

#require(pROC, quietly=TRUE, warn.conflicts=FALSE)

#require(e1071, quietly=TRUE, warn.conflicts=FALSE) 
#require(glmnet, quietly=TRUE, warn.conflicts=FALSE)
#require(randomForest, quietly=TRUE, warn.conflicts=FALSE)
#require(pROC, quietly=TRUE, warn.conflicts=FALSE)

"train.mwas" <- function(data.set, y=NULL, is.feat = FALSE, 
                         method=c("RF","SVM", "knn", "MLR")[1], 
                         #valid_type = c("cv", "jk")[1], 
                         nfolds = 10,
                         out.dir=NULL,
                         feat.threshold=0){
  
  if (class(data.set)=="mwas") {
    x <- data.set$features 
    y <- data.set$response
    
    is.feat <- data.set$is.feat
    if(is.null(is.feat)) is.feat=FALSE
    
    if(!is.null(data.set$feat_param)) {
      feat.threshold <- data.set$feat_param
    }else feat.threshold <- 0
    
    method <- data.set$method
    
    valid_type <- data.set$valid_type
    if(is.null(valid_type)) valid_type <- "cv"
    
    out.dir <- data.set$out.dir
  }else if(is.null(y)){
    stop("Response values are missing for the model training!")
  }else x <- data.set
  
  if (is.feat){
    if(is.null(out.dir)) out.dir <- "."
    file.out <- paste0(out.dir, "/feat_statistics")
    dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)
    feat.set <- feature.selection(x, y, selection_threshold = feat.threshold, "FDR", file.out)
    train.set <- x[,feat.set]
  }else {
    feat.set <- colnames(x)
    train.set <- x
  }
  best.model <- persist.model.mwas(train.set, y, nfolds=nfolds, classifier=method, 
                                   valid_type, is.feat = is.feat, feat.set=feat.set, 
                                   out.dir=out.dir)
  
  #export.mwas(trained.model = best.model, feat.set = feat.set)
  return(best.model)
}


# Nested cross-validation and Jack Knife validation for model performance estimation
# And model selection (model training)

"persist.model.mwas" <- function(x, y, nfolds=10, 
                                 classifier=c("RF","SVM", "knn", "MLR")[1],
                                 #valid_type=c("cv", "jk")[1], 
                                 is.feat=FALSE,
                                 feat.set=NULL, out.dir=NULL){
  
  #require(pROC, quietly=TRUE, warn.conflicts=FALSE)
  
  # x - feature set (observation * features)
  # y - desried response
  
  #cat.ind <- list()  # category indices
  #response.tab <- table(y)
  
  cv.ind <- sample(dim(x)[1])   # permutate the index
  cv.samp.num <- floor(length(cv.ind)/nfolds)
  sampl_ind <- seq(1, dim(x)[1], by=1)
  candidate.model <- list()
  candidate.rocobj <- list()
  candidate.auc <- vector("numeric", nfolds)
  candidate.error <- vector("numeric", nfolds)
  
  # nested cross-validation
  if (tolower(valid_type)=="cv"){ # case insensitive
    for (cv.id in 1:nfolds){
      # fold index that is being hold out
      if(cv.id < nfolds) {idx <- cv.ind[seq((cv.id-1)*cv.samp.num + 1, cv.id*cv.samp.num, by=1)]
      }else {# the last fold could contain less than cv.samp.num of observations
        idx <- cv.ind[seq((cv.id-1)*cv.samp.num + 1, length(cv.ind), by=1)] 
      }
      
      validation.set <- x[idx,]
      validation.labels <- y[idx]
      
      train.set <- x[!sampl_ind %in% idx,]
      train.labels <- y[!sampl_ind %in% idx]
      
      candidate <- cross.validation.mwas(train.set, train.labels, nfolds, classifier)
      candidate.eval <- model.evaluation.mwas(validation.set, candidate$best.model, validation.labels)
      
      candidate.model <- c(candidate.model, list(candidate))
      candidate.error[cv.id] <- candidate.eval$performance["error"]
      candidate.auc[cv.id] <- candidate.eval$performance["auc"]
      
    }
  }else if(valid_type=="jk") { # jackknife function 
    jk.ind <- seq(1, dim(x)[1], by=1)
    jk.samp.num <- floor(length(cv.ind)/nfolds)
    
    candidate.model <- list()
    candidate.rocobj <- list()
    candidate.auc <- vector("numeric", nfolds)
    candidate.error <- vector("numeric", nfolds)
    
    for (jk.fold in 1:nfolds){
      # fold index that is being hold out
      if(jk.fold < nfolds){
        idx <- jk.ind[seq((jk.fold-1)*jk.samp.num + 1, jk.fold*jk.samp.num, by=1)] 
      }else {# the last fold could contain less than jk.samp.num
        idx <- jk.ind[seq((jk.fold-1)*jk.samp.num + 1, length(jk.fold), by=1)]}
      
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
  class(best.model.eval) <- class(best.model)
  
  # if (savefile) save(best.model, file = paste(opts$outdir,"/trained.model", collapse='', sep=''))
  # Saving file is integrated into export.mwas.R
  if(!is.feat) feat.set <- NULL
  
  file.out <- sprintf('train_%s_model_results', classifier)
  export.mwas(trained.model=best.model, model.eval=best.model.eval, 
              trained.model.perform=model.perform, feat.set=feat.set,
              out.dir=out.dir, file.name = file.out)
  
  # class(best.model) <- "mwas"
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
  switch(classifier, rf = {
           #require(randomForest, quietly=TRUE, warn.conflicts=FALSE)
           
           #best.model <- tune.randomForest(x, y, tunecontrol = tune.control(random=TRUE, sampling="cross", cross = nfolds))
           
           best.model <- tune.randomForest(x, y, mtry=seq(from=min(round(sqrt(num_species)), round(num_species/5)), 
                                                          to=max(round(sqrt(num_species)), round(4*num_species/5)), 
                                                          by=5), 
                                         tunecontrol = tune.control(random=TRUE, sampling="cross", cross = 5))
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
         }, svm={
           #require(e1071, quietly=TRUE, warn.conflicts=FALSE) 

           best.model <- tune.svm(x, y, tunecontrol = tune.control(random=TRUE, sampling="cross", cross = nfolds), probability = TRUE)
           # $best.model         the model trained on the complete training data using the best parameter combination.
           # $best.parameters    a 1 x k data frame, k number of parameters
           # $best.performance   best achieved performance
           # $performances       a data frame with all parametere combinations along with the corresponding performance results
           # $train.ind          list of index vectors used for splits into training and validation sets
         }, mlr = {
           #require(glmnet, quietly=TRUE, warn.conflicts=FALSE)
           
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
           
         }, knn = {
           #require(e1071, quietly=TRUE, warn.conflicts=FALSE) 
           best.model <- tune.knn(x, y, k = seq(from=1, to=min(dim(x)[1]/nfolds*(nfolds-1), 15), by =2), 
                                  l= 1, sampling="cross", cross = nfolds, prob = TRUE) 
           class(best.model$best.model) <- 'knn'
           # $best.parameters$k          best k
           # $best.performance           classification error using best k
           # $performances$k             all candidate k's that used in cross-validation
           # $performances$error         errors for all candidates
           # $performances$dispersion    dispersions for all candidates
           # $best.model                 best model list 
         },
         stop("Undefined classification method!")
  )
  class(best.model) <- "mwas"
  return(best.model)
}

#"tune.ksvm" <- function()

# S3 method of predict.knn
"predict.knn" <- function(model, test){
  require(class)
  prediction <- knn(model$train, test, model$cl, k=model$k, l=model$l, prob=TRUE, use.all=TRUE)
  return(prediction)
}

"custom_dist_kernel" <- function(dataMat, method, param=1){
  switch(method,
         Matrix = {
           dist_kern <- as.kernelMatrix(exp(-dataMat/param))     
         },
         bray = {
           b_dist <- vegdist(dataMat, method = "bray")
           b_dist <- as.matrix(b_dist)
           dist_kern <- as.kernelMatrix(exp(-dataMat/param))
         },
         stop("Please specify distance type: Matrix (distance matrix) or bray (Bray-Curtis)")
  )
  return(dist_kern)
}

"custom.knn" <- function(data1, data2, y1, k){
  Ind1 <- dim(data1)[1]
  new_data <- rbind(data1, data2)
  dist_mat <-as.matrix(vegdist(new_data, method="bray"))
  Num <- dim(dist_mat)[1]
  y2 <- vector()
  for(id in (Ind1+1):Num){
    tt <- table(y1[match(names(sort(dist_mat[id, 1:Ind1])[1:k]), rownames(data1))])
    y2[id - Ind1] <- names(sort(tt, decreasing=T)[1])
  }
  return(y2)
}
