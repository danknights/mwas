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

"persist.model.mwas" <- function(x, y, nfolds=5, 
                                 classifier=c("RF","SVM", "knn", "MLR")[1],
                                 #valid_type=c("cv", "jk")[1], 
                                 is.feat=FALSE,
                                 feat.set=NULL, out.dir=NULL){
  
  #require(pROC, quietly=TRUE, warn.conflicts=FALSE)
  
  # x - feature set (observation * features)
  # y - desried response
  
  #cat.ind <- list()  # category indices
  #response.tab <- table(y)
  
  class.names <- levels(y)
  class.response <- list()
  class.length <- vector()
  cv.ind <- list()
  #xx <- matrix()
  for(id in seq_len(length(class.names))){
    class.response[[id]] <- which(y==class.names[id])
    class.length[id] <- length(class.response[[id]])
    #names(class.response[[id]]) <- class.names[id]
    cv.ind[[id]] <- sample(class.length[id])
  }
  
  #xx <- rbind()
  #cv.ind <- sample(dim(x)[1])   # permutate the index
  #cv.samp.num <- floor(length(cv.ind)/nfolds)  
  #sampl_ind <- seq(1, dim(x)[1], by=1)

  minSample.number <- min(class.length)  # even nubmers of samples for each class
  cv.samp.num <- floor(minSample.number/nfolds)
 
  # balanced feature set
  xx <- matrix()
  yy <- vector()
  for(id in seq_len(length(class.names))){
    if(id == 1){
      xx <- x[class.response[[id]][1:(cv.samp.num*nfolds)],]
    } else xx <- rbind(xx, x[class.response[[id]][1:(cv.samp.num*nfolds)],])
    yy <- c(yy, y[class.response[[id]][1:(cv.samp.num*nfolds)]])
  }
  sampl_ind <- seq(1, length(yy), by=1)
  # nested cross-validation
  #if (tolower(valid_type)=="cv"){ # case insensitive
 
  #}else if(valid_type=="jk") { # jackknife function 
  # Jackknifing for modle evaluation
  jk.ind <- seq(1, cv.samp.num*nfolds, by=1)
  #jk.samp.num <- floor(length(cv.ind)/nfolds)
  
  candidate.model <- list()
  candidate.rocobj <- list()
  candidate.auc <- vector("numeric", nfolds)
  candidate.error <- vector("numeric", nfolds)
  
  for (jk.fold in 1:nfolds){
    # fold index that is being hold out
    idx <- vector()
    for(id in seq_len(length(class.names))){
      idx <- c(idx, seq((jk.fold-1)*cv.samp.num + 1 + (id-1)*cv.samp.num*nfolds, 
                        jk.fold*cv.samp.num + (id-1)*cv.samp.num*nfolds, by=1))
    }
    
    validation.set <- xx[idx,]
    validation.labels <- yy[idx]
    
    train.set <- xx[!sampl_ind %in% idx,]
    train.labels <- yy[!sampl_ind %in% idx]
    
    candidate <- cross.validation.mwas(train.set, as.factor(train.labels), nfolds, classifier)
    candidate.model <- c(candidate.model, list(candidate))
    #candidate.error[jk.fold] <- candidate$best.performance
    
    pred.label <- predict(candidate$best.model, validation.set)
    candidate.error[jk.fold] <- length(which(pred.label != validation.labels))/length(validation.labels)
    pred.prob <- predict(candidate$best.model, validation.set, type='prob')
    
    candidate.obj <- roc.mwas(validation.set, predicted=pred[validation.labels], response=as.factor(validation.labels))
    #candidate.rocobj <- c(candidate.model, list(candidate.obj))
    candidate.auc[jk.fold] <- as.numeric(candidate.obj$auc)
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
  num_obs <- length(y)
  switch(classifier, rf = {
           #require(randomForest, quietly=TRUE, warn.conflicts=FALSE)
           
           #best.model <- tune.randomForest(x, y, tunecontrol = tune.control(random=TRUE, sampling="cross", cross = nfolds))
           
           best.model <- tune.randomForest(x, y, importance=TRUE,keep.forest=TRUE, mtry=c(round(sqrt(num_obs)/2), round(sqrt(num_obs)), round(2*sqrt(num_obs))),
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

           #best.model <- tune.ksvm(x, y, tunecontrol = tune.control(random=TRUE, sampling="cross", cross = nfolds), probability = TRUE)
           kernel <- tolower(kernel)
           sigma <- c(2^(-4:0), seq(from=2, to=2^8, by=5))
           if(kernel!='all'){
             
             best.model <- tune.ksvm(x, y, kernel, ranges=list(C=2^(-2:4), kpar=list(sigma=sigma)),
                                     probability = TRUE)
           }else{
             all_kernel <- c("linear","rbf_eu","rbf_bc", "rbf_uf")
             valid_obj <- list()
             #sigma <- c(2^(-4:0), seq(from=2, to=2^8, by=5))
             best_id <- 0
             count <- 0
             for(id in 1:length(all_kernel)){
               count <- count + 1
               valid.model <- tune.ksvm(x, y, all_kernel[id], ranges=list(C=2^(-2:4), kpar=list(sigma=sigma)),
                                       probability = TRUE)
               valid_obj <- c(valid_obj, valid.model)
               if(best_id==0){
                 best_id <- 1
               } else{
                 if(valid_obj[[best_id]]$best.performance > valid_obj$best.performance) best_id <- count
               }               
             }
             
             best.model <- valid_obj[[best_id]]
           }

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
           best.model <- tune.knn(x, y, k = seq(from=1, to=min(dim(x)[1]/nfolds*(nfolds-1), 25), by =2), 
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

"tune.ksvm" <- function(x, y, kernel, ranges=NULL, nfolds=5, probability = TRUE){
  
  kernel <- tolower(kernel)
 ###### grid search
 # Extract parameter ranges
 if(!is.null(ranges)){
   params <- names(ranges) # parameters to be tuned
   if("C" %in% params) C <- ranges$C
   if("kpar" %in% params) kpar <- ranges$kpar
 } else { # default
   C <- 1
   kpar <- NULL
 }
 
 cv_obj <- list()
 best_id <- 0
 count <- 0
 # grid search
 if(!is.null(kpar)){
   for(C_id in 1:length(C)){
     for(kpar_id in 1:length(kpar$sigma)){
       count <- count + 1
       
       switch(kernel,
              linear={
                valid_obj <- ksvm(x, y, type="C-svc", kernel="vanilladot", prob.model=probability, cross=nfolds,
                                  C=C[C_id])
              },
              rbf_eu={
                valid_obj <- ksvm(x, y, type="C-svc", kernel="rbfdot", prob.model=probability, cross=nfolds,
                                  C=C[C_id], kpar=list(sigma=kpar$sigma[kpar_id]))
              },
              rbf_bc={
                valid_obj <- ksvm(custom_dist_kernel(x, kernel, kpar$sigma[kpar_id]), y, type="C-svc",
                                  kernel="matrix", prob.model=probability, cross=nfolds, C= C[C_id])
              },
              rbf_uf={
                valid_obj <- ksvm(custom_dist_kernel(x, kernel, kpar$sigma[kpar_id]), y, type="C-svc",
                                  kernel="matrix", prob.model=probability, cross=nfolds, C= C[C_id])
              },
              stop("'kernel' should be one of 'linear','rbf_eu','rbf_bc','rbf_uf', 'all' ")
         )
       # ..@error - training error
       # ..@cross - cross validation error
       # ..@kernelf@kpar$sigma - sigma
       # ..@obj - Objective Function Value
       # ..@param$C - C
       # ..@nSV - nubmer of support vectors
       
       cv_obj <- c(cv_obj, valid_obj)
       if(best_id==0){
         best_id <- 1
       } else{
         if(cv_obj[[best_id]]@cross > valid_obj@cross) best_id <- count
       }
     }
   }
 } else {
   for(C_id in 1:length(C)){
     count <- count + 1
     
     switch(kernel,
            linear={
              valid_obj <- ksvm(x, y, type="C-svc", kernel="vanilladot", prob.model=probability, cross=nfolds,
                                C=C[C_id])
            },
            rbf_eu={
              valid_obj <- ksvm(x, y, type="C-svc", kernel="rbfdot", prob.model=probability, cross=nfolds,
                                C=C[C_id], kpar=list(sigma=kpar$sigma[kpar_id]))
            },
            rbf_bc={
              valid_obj <- ksvm(custom_dist_kernel(x, "bc", kpar$sigma[kpar_id]), y, type="C-svc",
                                kernel="matrix", prob.model=probability, cross=nfolds, C= C[C_id])
            },
            rbf_uf={
              valid_obj <- ksvm(custom_dist_kernel(x, "uf", kpar$sigma[kpar_id]), y, type="C-svc",
                                kernel="matrix", prob.model=probability, cross=nfolds, C= C[C_id])
            },
            stop("'kernel' should be one of 'linear','rbf_eu','rbf_bc','rbf_uf','all' ")
     )
     
     #valid_obj <- ksvm(x, y, type="C-svc", kernel=kernel, prob.model=probability, cross=nfolds, 
     #                   C=C[C_id])
     # ..@error - training error
     # ..@cross - cross validation error
     # ..@kernelf@kpar$sigma - sigma
     # ..@obj - Objective Function Value
     # ..@param$C - C
     # ..@nSV - nubmer of support vectors
     # ..
     
     cv_obj <- c(cv_obj, valid_obj)
     if(best_id==0){
       best_id <- 1
     } else{
       if(cv_obj[[best_id]]@cross > valid_obj@cross) best_id <- count
     }
   }
 }
 best.model <- list()
 best.model$best.model <- cv_obj[[best_id]]
 best.model$best.parameteris <- list(C=cv_obj[[best_id]]@param$C, sigma=cv_obj[[best_id]]@kernelf@kpar$sigma)
 best.model$best.performance <- cv_obj[[best_id]]@cross
 
 return(best.model)
}

# S3 method of predict.knn
"predict.knn" <- function(model, test){
  require(class)
  prediction <- knn(model$train, test, model$cl, k=model$k, l=model$l, prob=TRUE, use.all=TRUE)
  return(prediction)
}

"custom_kernel" <- function(dataMat, method, param=1){
  method <- tolower(method)
  switch(method,
         uf = { # UniFrac distance
           # dataMat must be a UniFrac distance matrix
           dist_kern <- as.kernelMatrix(exp(-dataMat/param))
         },
         bc = { # Bray-Curtis distance
           # dataMat must be a raw data table
           b_dist <- vegdist(dataMat, method = "bray")
           b_dist <- as.matrix(b_dist)
           dist_kern <- as.kernelMatrix(exp(-dataMat/param))
         },
         stop("Please specify distance type: Matrix (distance matrix) or bray (Bray-Curtis)")
  )
  return(dist_kern)
}

# "custom.knn" <- function(data1, data2, y1, k){
#   Ind1 <- dim(data1)[1]
#   new_data <- rbind(data1, data2)
#   dist_mat <-as.matrix(vegdist(new_data, method="bray"))
#   Num <- dim(dist_mat)[1]
#   y2 <- vector()
#   for(id in (Ind1+1):Num){
#     tt <- table(y1[match(names(sort(dist_mat[id, 1:Ind1])[1:k]), rownames(data1))])
#     y2[id - Ind1] <- names(sort(tt, decreasing=T)[1])
#   }
#   return(y2)
# }