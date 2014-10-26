
"feature.selection.mwas" <- function(x, y, feature.ids, selection_thres = 1, filename='feature_importance_scores.txt', outdir='.'){

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
