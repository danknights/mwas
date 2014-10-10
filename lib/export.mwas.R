# Export results to files
#
# --- input:
#   trained.model : trained claissifier model and other parameters
#        feat.set : selected feature vector index
#     test.results: output from predict.mwas.R, it's an model.evaluation object
#    resutls.plot : parameters for plotting
#            opts : options from keyboard
#
#------ ouput
#     export the required results to the designated directory (opts$outdir)
#
"export.mwas" <- function(trained.model, feat.set, test.results, results.plot, opts){
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
    save.results.svm(test.results$prediction, opts)
  }  
  if (exists(results.plot)){
    # save plots required by users (opts$visualize)
    # 
    #
    #
    #
    #
  }
}

#
#
# 
"save.results.svm" <- function(pred.obj, opts){
  filepath <- sprintf('%s/prediction_labels_likelihood.txt', opts$outdir)
  results.table <- merge(pred.obj$predicted, attr(pred.obj$predicted, "probabilities"))
  sink(filepath)
  cat('#SampleIndex\t')
  write.table(results.table,sep='\t',quote=F)
  sink(NULL)
}