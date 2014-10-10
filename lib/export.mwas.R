# Export results to files
#
# --- input:
#   trained.model : trained claissifier model and other parameters
#        feat.set : selected feature vector index
#     test.results: output from model.evaluation
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
    save(best.model, file = paste(opts$outdir,"/trained.model", collapse='', sep=''))    
  }
  
  if (exists(test.results)){
    # save classification accuracy, if there is acc
    # save predicted labels  
    # save likelihood probabilities, if there is any
    # save AUC, if there is any
    # save MCC, if there is any
    # save Cohen's Kappa, if there is any
    
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