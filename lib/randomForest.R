### for Random Forests classifier
# Runs standard random forests with out-of-bag error estimation
# Return values are the same as ml.cross.validation
"rf.out.of.bag" <- function(x,y, verbose=verbose, ...){
  rf.model <- randomForest(x,y,keep.inbag=TRUE,importance=TRUE,do.trace=verbose,...)
  result <- list()
  result$probabilities <- get.oob.probability.from.forest(rf.model,x)
  result$y <- y
  result$predicted <- rf.model$predicted
  result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
  result$params <- list(ntree=opts$ntree)
  result$errs <- as.numeric(result$predicted != result$y)
  result$importances <- rf.model$importance[,'MeanDecreaseAccuracy']
  return(result)
}

# Print "feature importance scores" file
"save.rf.results.importances" <- function(result, feature.ids, filename='feature_importance_scores.txt', outdir='.'){
  filepath <- sprintf('%s/%s',outdir,filename)
  if(is.null(dim(result$importances))){
    imp <- result$importances
    imp.sd <- rep(NA,length(imp))
  } else {
    imp <- rowMeans(result$importances)
    imp.sd <- apply(result$importances, 1, sd)
  }
  output.table <- cbind(imp, imp.sd)
  rownames(output.table) <- feature.ids
  output.table <- output.table[sort(imp,dec=T,index=T)$ix,]
  colnames(output.table) <- c('Mean_decrease_in_accuracy','Standard_deviation')
  
  sink(filepath)
  cat('Feature_id\t')
  write.table(output.table,sep='\t',quote=F)
  sink(NULL)
}

