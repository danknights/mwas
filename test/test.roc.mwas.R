# Receiver operating characteristic - using package 'pROC' or 'ROCR'
# ----- input:
#         x: feature vector
# predicted: predicted output using the trained model
#   desired: desired response 
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
"test.roc.mwas" <- function(x, predicted, response, is.plot=TRUE){
  #
  #
  #
}