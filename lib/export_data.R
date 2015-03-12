# Export results to files
#
# Contributors: Hu, Dan
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

# require(xlsx, quietly = TRUE, warn.conflicts=FALSE)
# if (!require("xlsx", quietly=TRUE, warn.conflicts = FALSE)) {
#  install.packages("xlsx", dependencies = TRUE)
#  library(xlsx, verbose=F, warn.conflicts =F)
# }

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
    file.name <- sprintf('%s/%s.txt', out.dir, file.name)
    
    # save as .txt format
    scipen.save <- options('scipen') 
    options(scipen=20)                         # avoid exponential notation
    sink(file.name)
    cat('The Jackknife model evaluation: \n')
    cat('\t Error rate: ', model.eval$mean.error, " +/- ", model.eval$std.error, "\n")
    cat('\t AUC from ROC: ', model.eval$mean.auc, " +/- ", model.eval$std.auc, "\n")
    sink(NULL)
    options(scipen=scipen.save)
    #save.xlsx(objects=model.eval, file.name=file.name)
  }  

  if (!is.null(trained.model.perform)){
    # save trained model evaluation
    file.name1 <- sprintf('%s/trained_%s_model_perform_labels.txt', out.dir, class(trained.model))
    #save.xlsx(objects=trained.model.perform, file.name=file.name)
    pred.table <- cbind(trained.model.perform$prediction, trained.model.perform$probabilities)
    rownames(pred.table) <- trained.model.perform$rownames
    #colnames(pred.table) <- c("labels", trained.model.perform$colnames)
    sink(file.name1)
    write.table(pred.table,quote=F,sep='\t')
    sink(NULL)
    
    file.name2 <- sprintf('%s/trained_%s_model_performance.txt', out.dir, class(trained.model))
    
    perform_table <- list()
    perform_table$confusion.matrix <- trained.model.perform$confusion.matrix
    perform_table$performance <- trained.model.perform$performance
    
    sink(file.name2)
    lapply(perform_table, print, file=file.name2, append=TRUE)
    sink(NULL)
  }  
   
  if (!is.null(feat.set)){
    feat_file_name <- sprintf('%s/selected_feature_sets.txt', out.dir)
    file.out <- file(feat_file_name, 'w')
    #names(feat.set) <- "Features"
    #cat("ID\tFeatures - ")
    write.table(feat.set, file.out, sep='\t')
    flush(file.out)
    close(file.out)
  }
}

# save data tabel in a Excel workbook.
# The original code is by Rob Kabacoff 
# from: http://www.r-bloggers.com/quickly-export-multiple-r-objects-to-an-excel-workbook/
# 
# "save.xlsx" <- function (objects, file.name){
# 
#  # objects <- list(...)
#   fargs <- as.list(match.call(expand.dots = TRUE))
#   objnames <- as.character(fargs)[-c(1, 2)]
#   objnames <- names(objects)
#   nobjects <- length(objects)
#   for (i in 1:nobjects) {
#     if (i == 1) {
#       write.xlsx(objects[[i]], file.name, sheetName = objnames[i])
#     }else write.xlsx(objects[[i]], file.name, sheetName = objnames[i], append = TRUE)
#   }
#   print(paste("Workbook", file.name, "has", nobjects, "worksheets."))
# }


# saves list of results from feature.statistics to file (or prints)
"write.statistical.test.results" <- function(results, out.dir=NULL, filename='feature_statistics'){
  
  if(is.null(out.dir)) out.dir <- '.'
  
  # save as .txt format
  scipen.save <- options('scipen') 
  options(scipen=20)                         # avoid exponential notation
  hits <- cbind(results$pvalues, results$qvalues)
  hits <- cbind(hits, results$classwise.means)
  colnames(hits)[1:2] <- c('pvalue','qvalue')
  hits <- hits[!is.na(hits[,1]),,drop=F]     # remove all NA values
  hits <- hits[order(hits[,1]),]
  filename1 <- sprintf("%s/%s.txt", out.dir, filename)
  sink(filename1)
  cat('Features\t')
  write.table(hits,quote=F,sep='\t')
  sink(NULL)
  options(scipen=scipen.save)
  
  # save as .xlsx format
  #filename2 <- sprintf("%s.xlsx", filename)
  #hits <- cbind(rownames(hits), hits)
  #colnames(hits)[1] <- "Features"
  #write.xlsx(hits, filename2, row.names=FALSE, sheetName = results)
  
}

