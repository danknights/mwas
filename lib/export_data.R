# Export results to files
#
# Contributors: Hu
# --- input:
#     trained.model : trained claissifier model and other parameters
# trained.model.eval: model evaluation on the training set by using the trained model 
#                     (error, confusion matrix, AUC, MCC, Kappa)
#     model.perform : model performance estimation (mean +/- std error/auc)
#          feat.set : selected feature vector index
#       test.results: output from predict.mwas.R, it's an model.evaluation object
#                     if there is no desired labels are given then it only contains 
#              opts : options from users
#  
# ------ ouput
#     export the required results to the designated directory (opts$outdir)
# ------ 
# Last update: 10/25/2014
#

"export.mwas" <- function(trained.model=NULL, model.eval=NULL, trained.model.perform=NULL,
                          feat.set=NULL, out.dir=NULL, file.name="predcition_results"){
  
  if(is.null(out.dir)) out.dir <- '.'
  
  if(!is.null(trained.model)) {
    # save trained model and selected feature vector
    if(!is.null(feat.set)){
      # save feautre scores for each feature
      # save selected feature vector with names (OTU ID etc.)
      #
      
      best.model<-list()
      best.model$trained.model <- trained.model
      best.model$features <- feat.set      # feature row names or index
    } else {
      best.model<-list()
      best.model$trained.model <- trained.model
    }
   # should save as rds format
   model.out.name <- sprintf('%s/trained_%s_model.rds', out.dir, class(trained.model))
   #saveRDS(best.model, file = paste(out.dir,"/trained_model.rds", collapse='', sep=''))
   saveRDS(best.model, file = model.out.name)
  }
  
  if (!is.null(model.eval)){
    # save predicted labels and/or likelihood probabilities
    # save classification accuracy, if there is acc 
    # save likelihood probabilities, if there is any
    # save AUC, if there is any
    # save MCC, if there is any
    # save Cohen's Kappa, if there is any
    file.name <- sprintf('%s/%s.xlsx', out.dir, file.name)
    save.xlsx(objects=model.eval, file.name=file.name)
  }  

  if (!is.null(trained.model.perform)){
    # save trained model evaluation
    file.name <- sprintf('%s/trained_%s_model_performance.xlsx', out.dir, class(trained.model))
    save.xlsx(objects=trained.model.perform, file.name=file.name)
  }  
   
  if (!is.null(feat.set)){
    feat_file_name <- sprintf('%s/selected_feature_sets', out.dir)
    file.out <- file(feat_file_name, 'w')
    write.table(feat.set, file.out, sep='\t')
    flush(file.out)
    close(file.out)
  }
}

# save data tabel in a Excel workbook.
# The original code is by Rob Kabacoff 
# from: http://www.r-bloggers.com/quickly-export-multiple-r-objects-to-an-excel-workbook/
# 
"save.xlsx" <- function (objects, file.name){
  require(xlsx, quietly = TRUE, warn.conflicts=FALSE)
  
  #objects <- list(...)
  #fargs <- as.list(match.call(expand.dots = TRUE))
  #objnames <- as.character(fargs)[-c(1, 2)]
  objnames <- names(objects)
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (i == 1) {
      write.xlsx(objects[[i]], file.name, sheetName = objnames[i])
    }else write.xlsx(objects[[i]], file.name, sheetName = objnames[i], append = TRUE)
  }
  print(paste("Workbook", file.name, "has", nobjects, "worksheets."))
}


# saves list of results from differentiation.test to file (or prints)
"write.differentiation.test.results" <- function(results, filename='differentiated.features.txt'){
  if(!is.null(filename)){
    scipen.save <- options('scipen')
    options(scipen=20)
    hits <- cbind(results$pvalues, results$qvalues)
    hits <- cbind(hits, results$classwise.means)
    colnames(hits)[1:2] <- c('pvalue','qvalue')
    hits <- hits[!is.na(hits[,1]),,drop=F]
    hits <- hits[order(hits[,1]),]
    sink(filename)
    cat('Feature\t')
    write.table(hits,quote=F,sep='\t')
    sink(NULL)
    options(scipen=scipen.save)
  }
}

