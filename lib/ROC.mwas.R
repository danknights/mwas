"roc.mwas" <- function(x, predicted, response, is.plot=TRUE){
  # Receiver operating characteristic - using package 'pROC' or 'ROCR'
  # ----- input:
  #         x: feature vector
  # predicted: predicted output using the trained model
  #   desired: desired output 
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
  if (length(levels(response)) == 2)  # binary classification
    rocobj <- roc(response, as.numeric(predicted), percent=TRUE, ci=TRUE, plot=is.plot)
  else    # multi-class classification
    rocobj <- multiclass.roc(response, as.numeric(predicted), percent=TRUE, ci=TRUE, plot=is.plot)
  
  return(rocobj)
}