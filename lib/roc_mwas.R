# Receiver operating characteristic - using package 'pROC' or 'ROCR'
# Contributors: Hu
# 
# ----- input:
#         x: feature vector
#     model: trained model
# predicted: predicted output using the trained model (if no model is specified as input)
#   response: desired response 
#
# ----- output:
#  rocobj:  ROC object
#           $auc              class "auc" 
#           $sensitivities    sensitivities defining the ROC curve
#           $specificities    specificities defining the ROC curve         
#           $response         response vector (desired)
#           $predictor        the predictor vector converted to numeric as used to build the ROC curve
# ----- 
# Last update: 10/25/2014
#

"roc.mwas" <- function(x, model=NULL, predicted=NULL, response, is.plot=FALSE){

#  require(pROC, quietly=TRUE, warn.conflicts=FALSE)
  
  if (!is.null(model)&&is.null(predicted)) { 
    predicted <- predict(model, x) 
  }else if(is.null(model)&&is.null(predicted)){ 
    # if both model and predicted values are not provided, then informm the user and stop.
    stop("ROC fucntion needs at least one of the following values: trained model or predicted values!")
  }
  
  roc.flag <- 1
  for(id in length(levels(response))){
    roc.flag <- roc.flag * length(which(response==levels(response)[id]))
  }
  
  if(roc.flag != 0){
    if (length(levels(response)) == 2)  { # binary classification
      rocobj <- roc(response, as.numeric(predicted), percent=TRUE, ci=TRUE, plot=is.plot) 
    } else {# multi-class classification
      rocobj <- multiclass.roc(response, as.numeric(predicted), percent=TRUE, plot=is.plot)
    }
  } else {
    rocobj <- list()
    rocobj$auc <- 0
    #rocobj <- NULL
  }
  return(rocobj)
}

