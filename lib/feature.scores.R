"mwas.feature.scores" <- function(x, y, feature.ids, filename='feature_importance_scores.txt', outdir='.'){
  result <- rf.out.of.bag(x, y)

  save.rf.results.importances(result, feature.ids, filename, outdir)
  imp <- varimp(result$rf.model, mode=1, scale=T) # scaled is great for feature picking
  print(imp) # just to make sure the table looks as expected
  features <- imp[order(imp[:,1],decreasing=T),:]
  
  # Loop through the features in order of importance, and keep grabbing them until
  # they are no longer important (threshold > 1)
  i <- 0
  endFeatures = NULL;
  
  while (features[i,1] > 1) { # features over 1 importance are "decently" important (?)
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
