# Feature selection using randomForest
# This function is adapted from randomForest module in QIIME pipeline
# Contributors: Hu, Dan
#---
#  input:
#         x : feature vector
#  response : response
#   selection_threshold : threshold for feature selection (determines the number of feautres)
#     method: feature selection criterion; 
#             "FDR" - based on the statistical test; 
#                     selected features have smaller false discovery rate than selection_threshold
#             "RF"  - based on random Forest feature importance values
#                     selected features have larger impartance value than selection_threshold
# ---
#  output:
#   list of feature vector
#   $feautres: selected feature vector (values)
#   $id      : selected feautre index (colnames) in the feature vector
# --- 
# Last update: 03/01/2015
#

"feature.selection" <- function(x, response, selection_threshold = 0.05, method="FDR", out.dir = NULL){
  
  feat_set <- list()
  
  if(is.null(out.dir)) out.dir <- "."
  file_name = sprintf("%s/%s_%.2f_feature_statistcs.txt", out.dir, method, selection_threshold)
  
  method <- toupper(method)
  switch(method,
         FDR = {
           feat.stats <- feature.statistics(x, response, selection_threshold, include.subset = TRUE)
           feat_set$id <- colnames(feat.stats$subset)
           feat_set$features <- x[,feat_set$id]
           
           cat(length(feat_set$id), " features are selected under the criterion FDR < ", selection_threshold)
           
           file_name <- gsub(".txt", "", file_name)
           write.statistical.test.results(feat.stats, filename=file_name)
         }, 
         RF = {
           rf.model <- randomForest(x, response, proximity = TRUE, importance=TRUE)
           #importances <- rf.model$importance[,'MeanDecreaseAccuracy']
           imp <- importance(rf.model, type =1, scale=T)
           importances_order <- order(imp, decreasing = T)
           
           # save feature importance list
           file.out <- file(file_name, 'w')
           write.table(imp, file.out, sep='\t')
           flush(file.out)
           close(file.out)
           
           ordered_feat_set <- colnames(x)[importances_order]
           ordered_feat_imp <- imp[importances_order]
           
           feat_set$id <- ordered_feat_set[ordered_feat_imp > selection_threshold]
           feat_set$features <- x[, feat_set$id]
           
           cat(length(feat_set$id), " features are selected under the criterion feat_importance > ", selection_threshold)
         }
         
  )
  return(feat_set)
}


  
