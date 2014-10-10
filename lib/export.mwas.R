# Export results to files
#
# --- input:
#   trained.model : trained claissifier model and other parameters
#        feat.set : selected feature vector index
#     test.results: output from predict.mwas.R, it's an model.evaluation object
#                   if there is no desired labels are given then it only contains 
#            opts : options from keyboard
#
#------ ouput
#     export the required results to the designated directory (opts$outdir)
#
"export.mwas" <- function(trained.model, feat.set, test.results, opts, ...){
  if(exists(trained.model)) {
    # save trained model and selected feature vector
    if(exists(feat.set)){
      # save feautre scores for each feature
      # save selected feature vector with names (OTU ID etc.)
      #
      
      best.model<-list()
      best.model$trained.model <- trained.model
      best.model$features <- feat.set.ix # feature row names or index
    }
    else {
      best.model<-list()
      best.model$trained.model <- trained.model
    }
   # should save as rdt format
    save(best.model, file = paste(opts$outdir,"/trained.model", collapse='', sep=''))    
  }
  
  if (exists(test.results)){
    # save predicted labels and/or likelihood probabilities
    # save classification accuracy, if there is acc 
    # save likelihood probabilities, if there is any
    # save AUC, if there is any
    # save MCC, if there is any
    # save Cohen's Kappa, if there is any
    save.results(test.results, opts, ...)
  }  

}

#
# save results from SVM prediction object
# input: 
#    pred.obj : prediction object
#        opts : options from keyboard
# output:
#   save predicted labels, likelihood. 
#   If the desired response is given, then also output confusion matrix and AUC, MCC and Kappa
# 
"save.results.svm" <- function(pred.obj, opts){
  if (exists("desired"))
  filepath <- sprintf('%s/svm_prediction_labels_likelihood.txt', opts$outdir)
  results.table <- merge(pred.obj$predicted, attr(pred.obj$predicted, "probabilities"))
  sink(filepath)
  cat('#SampleIndex\t')
  write.table(results.table,sep='\t',quote=F)
  sink(NULL)
}

#
# save results from MLR prediction object
# input: 
#    pred.obj : prediction object
#        opts : options from keyboard
# output:
#   save predicted labels, likelihood. 
#   If the desired response is given, then also output confusion matrix and AUC, MCC and Kappa
# 
"save.results.rf" <- function(pred.obj, opts){
  filepath <- sprintf('%s/rf_prediction_labels_likelihood.txt', opts$outdir)
  results.table <- merge(pred.obj$predicted, attr(pred.obj$predicted, "probabilities"))
  sink(filepath)
  cat('#SampleIndex\t')
  write.table(results.table,sep='\t',quote=F)
  sink(NULL)
}

#
# save results from MLR prediction object
# input: 
#    pred.obj : prediction object
#        opts : options from keyboard
# output:
#   save predicted labels, likelihood. 
#   If the desired response is given, then also output confusion matrix and AUC, MCC and Kappa
# 
"save.results.mlr" <- function(pred.obj, opts){
  filepath <- sprintf('%s/MLR_prediction_labels_likelihood.txt', opts$outdir)
  results.table <- merge(pred.obj$predicted, attr(pred.obj$predicted, "probabilities"))
  sink(filepath)
  cat('#SampleIndex\t')
  write.table(results.table,sep='\t',quote=F)
  sink(NULL)
}

#
# save results from knn prediction object
# input: 
#    pred.obj : prediction object
#        opts : options from keyboard
# output:
#   save predicted labels, likelihood. 
#   If the desired response is given, then also output confusion matrix and AUC, MCC and Kappa
# 
"save.results.knn" <- function(pred.obj, opts){
  filepath <- sprintf('%s/knn_prediction_labels_likelihood.txt', opts$outdir)
  results.table <- merge(pred.obj$predicted, attr(pred.obj$predicted, "probabilities"))
  sink(filepath)
  cat('#SampleIndex\t')
  write.table(results.table,sep='\t',quote=F)
  sink(NULL)
}