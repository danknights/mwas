# Feature selection using randomForest
# This function is adapted from randomForest module in QIIME pipeline
# Contributors: Hu, Gabe, Dan
#---
#  input:
#     x : feature vector
#     y : response
#   selection_threshold : threshold for feature selection (determines the number of feautres)
# ---
#  output:
#   list of feature vector
#   $feautres: selected feature vector
#   $ix      : selected feautre index in the feature vector
# --- 
# Last update: 10/25/2014
#

"feature.scores.mwas" <- function(x, y, selection_threshold = 1, out.dir = NULL){
  
  #require(caret, quietly=TRUE, warn.conflicts=FALSE)
  require(randomForest, quietly=TRUE, warn.conflicts=FALSE)

  rf.model <- randomForest(x,y, proximity = TRUE, importance=TRUE)
  #importances <- rf.model$importance[,'MeanDecreaseAccuracy']
  imp <- importance(rf.model, type =1, scale=T)
  importances_order <- order(imp, decreasing = T)
  
  if (is.null(out.dir)) {
    file_name = './feature_scores.txt'
  } else file_name = sprintf('%s/feature_scores.txt', out.dir)
  
  file.out <- file(file_name, 'w')
  write.table(imp, file.out, sep='\t')
  flush(file.out)
  close(file.out)
  
  ordered_feat_set <- colnames(x)[importances_order]
  ordered_feat_imp <- imp[importances_order]
  
  feat_set <- ordered_feat_set[ordered_feat_imp > selection_threshold]
  
  return(feat_set)
}

"features.mRMR.mwas" <- function(x, y, selection_threshold = 1, out.dir = NULL){
  require(mRMRe, quietly=TRUE, warn.conflicts=FALSE)
  
}
