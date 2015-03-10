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
  
  num_obs <- length(response)
  
  method <- tolower(method)
  switch(method,
         fdr = {
           file_name = sprintf("%s/%s_%.2f_feature_statistcs", out.dir, method, selection_threshold)
           
           feat.stats <- feature.statistics(x, response, selection_threshold, include.subset = TRUE)
           feat_set$id <- colnames(feat.stats$subset)
           feat_set$features <- x[,feat_set$id]
           #feat_set$scores <- 
           
           cat(length(feat_set$id), " features are selected under the criterion FDR < ", selection_threshold)
           
           #file_name <- gsub(".txt", "", file_name)
           write.statistical.test.results(feat.stats, filename=file_name)
           
           ## TODO: determine the number of features for the reduced feature set
           
         }, 
         rf = {
           file_name = sprintf("%s/%s_%.2f_feature_importance_rank.txt", out.dir, method, selection_threshold)
           
           #model.rf <- tune.randomForest(x, y, importance=TRUE, mtry=seq(from=min(round(sqrt(num_species)), round(num_species/5)), to=max(round(sqrt(num_species)), round(4*num_species/5)), by=5), 
           model.rf <- tune.randomForest(x, as.factor(y), importance=TRUE, mtry=c(round(sqrt(num_obs)/2), round(sqrt(num_obs)), round(2*sqrt(num_obs))), 
                                         tunecontrol = tune.control(random=TRUE, sampling="cross", cross = 5))
           summary(model.rf)
           
           #rf.model <- randomForest(x, response, proximity = TRUE, importance=TRUE)
           #importances <- rf.model$importance[,'MeanDecreaseAccuracy']
           imp <- importance(model.rf$best.model, type =1, scale=T) 
           importances_order <- order(imp, decreasing = T)
           
           ordered_feat_set <- colnames(x)[importances_order]
           ordered_feat_imp <- imp[importances_order]
           
           ordered_feat_list <- as.matrix(ordered_feat_imp, dimnames=list(ordered_feat_set))
           ordered_feat_list <- cbind(ordered_feat_list, model.rf$best.model$importance[importances_order,length(levels(y))+1])
           colnames(ordered_feat_list) <- c("Importance_value", "MeanDecreaseAccuracy")
           
           # save feature importance list
           file.out <- file(file_name, 'w')
           #write.table(imp, file.out, sep='\t')
           #flush(file.out)
           #close(file.out)
           sink(file.out)
           cat('Features\t')
           write.table(ordered_feat_list, quote=F, sep='\t')
           sink(NULL)
           
           feat_set$id <- ordered_feat_set[ordered_feat_imp > selection_threshold]
           feat_set$features <- x[, feat_set$id]
           
           cat(length(feat_set$id), " features are selected under the criterion feat_importance > ", selection_threshold)
         }
         
  )
  return(feat_set)
}

"feature.number" <- function(feat_set){
  
}  
