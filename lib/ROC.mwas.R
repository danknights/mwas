"roc.mwas" <- function(x, predicted, response, binaryClass=TRUE, is.plot=TRUE){
  # Receiver operating characteristic - using package 'pROC' or 'ROCR'
  # ----- input:
  #         x: feature vector
  # predicted: predicted output using the trained model
  #   desired: desired output
  #  binaryClass: two-class classification or multiclass problem 
  #
  # ----- output:
  #  rocobj:  ROC object
  #           $auc              class "auc" 
  #           $sensitivities    sensitivities defining the ROC curve
  #           $specificities    specificities defining the ROC curve         
  #           $response         response vector (desired)
  #           $predictor        the predictor vector converted to numeric as used to build the ROC curve
  #
  #
  
  # predicted <- predict(model, x)
  if(binaryClass)  
    rocobj <- roc(response, as.numeric(predicted), percent=TRUE, ci=TRUE, plot=is.plot)
  else 
    rocobj <- multiclass.roc(response, as.numeric(predicted), percent=TRUE, ci=TRUE, plot=is.plot)
  
  return(rocobj)
}