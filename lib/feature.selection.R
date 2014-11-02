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

"feature.scores.mwas" <- function(x, y, selection_thres = 1){
  
  require(caret, quietly=TRUE, warn.conflicts=FALSE)
  require(randomForest, quietly=TRUE, warn.conflicts=FALSE)

  rf.model <- randomForest(x,y,keep.inbag=TRUE,importance=TRUE)
  
  features.scores <- importance(rf.model, type=1, scale=FALSE)
  
  ff <- varImp(rf.model, mode =1)
  importances <- rf.model$importance[,'MeanDecreaseAccuracy']
  varImpPlot(rf.model)
  # a matrix with nclass + 2 (for classification) or two (for regression) columns. 
  # For classification, the first nclass columns are the class-specific measures 
  # computed as mean descrease in accuracy. The nclass + 1st column is the mean 
  # descrease in accuracy over all classes. The last column is the mean decrease 
  # in Gini index. For Regression, the first column is the mean decrease in accuracy 
  # and the second the mean decrease in MSE. If importance=FALSE, the last measure is 
  # still returned as a vector.
  features.scores <- importances[order(importances, decreasing=T)]
  
  # Loop through the features in order of importance, and keep grabbing them until
  # they are no longer important (threshold > 1)
  i <- 0
  endFeatures = NULL
  
  #while (features.scores[i,1] >= selection_threshold) { # features over threshold importance are kept
  #	endFeatures <- rbind(endFeatures,features[i,:])
  #	i <- i + 1
  #}
  endIx <- match(rownames(endFeatures),colnames(x))
  
  return(list(features = endFeatures, ix = endIx))
}

