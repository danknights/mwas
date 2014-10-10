"mwas.feature.scores" <- function(...) UseMethod("mwas.feature.scores") {

}

"mwas.feature.selection" <- function(x, y, feature.ids, selection_thres = 1, filename='feature_importance_scores.txt', outdir='.'){
=======
# Feature selection using randomForest
# ---
#  input:
#     x : feature vector
#     y : response
#   selection_threshold : threshold for feature selection (determines the number of feautres)
# ---
#  output:
#   list of feature vector
#   $feautres: selected feature vector
#   $ix      : selected feautre index in the feature vector
#
  result <- rf.out.of.bag(x, y)

  imp <- varimp(result$rf.model, mode=1, scale=T) # scaled is great for feature picking
  print(imp) # just to make sure the table looks as expected
  features <- imp[order(imp[:,1],decreasing=T),:]
  
  # Loop through the features in order of importance, and keep grabbing them until
  # they are no longer important (threshold > 1)
  i <- 0
  endFeatures = NULL;
  
  while (features[i,1] >= selection_threshold) { # features over threshold importance are kept
  	endFeatures <- rbind(endFeatures,features[i,:])
  	i <- i + 1
  }
  endIx <- match(rownames(endFeatures),colnames(x))
  
  return(list(features = endFeatures, ix = endIx))
}


"rf.out.of.bag" <- function(x,y, verbose=verbose, ...){
  rf.model <- randomForest(x,y,keep.inbag=TRUE,importance=TRUE,do.trace=verbose,...) # scale = TRUE
  result <- list()
  result$probabilities <- get.oob.probability.from.forest(rf.model,x)
  result$y <- y
  result$predicted <- rf.model$predicted
  result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
  result$params <- list(ntree=opts$ntree)
  result$errs <- as.numeric(result$predicted != result$y)
  result$importances <- rf.model$importance[,'MeanDecreaseAccuracy']
  result$rf.model <- rf.model
  return(result)
}

# get probability of each class using only out-of-bag predictions from RF
"get.oob.probability.from.forest" <- function(model,x){
  # get aggregated class votes for each sample using only OOB trees
  votes <- get.oob.votes.from.forest(model,x)
  # convert to probs
  probs <- sweep(votes, 1, apply(votes, 1, sum), '/')
  rownames(probs) <- rownames(x)
  colnames(probs) <- model$classes
  
  return(invisible(probs))
}

# get votes for each class using only out-of-bag predictions from RF
"get.oob.votes.from.forest" <- function(model,x){
  # get aggregated class votes for each sample using only OOB trees
  votes <- matrix(0, nrow=nrow(x), ncol=length(model$classes))
  
  rf.pred <- predict(model, x, type="vote",predict.all=T)
  for(i in 1:nrow(x)){
    # find which trees are not inbag for this sample
    outofbag <- model$inbag[i,]==0
    # get oob predictions for this sample
    votes[i,] <- table(factor(rf.pred$individual[i,][outofbag],levels=model$classes))
  }
  rownames(votes) <- rownames(x)
  colnames(votes) <- model$classes
  
  return(invisible(votes))
}
