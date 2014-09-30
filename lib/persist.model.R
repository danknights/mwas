"mwas.persist.model" <- function(x, y, nfolds=10, 
                               classifier=c("RF","SVM", "knn", )[1], 
                               savefile = TRUE, opts, ...){
  # x - feature set (observation * features)
  # y - desried response
  cv.ind <- sample(dim(x)[1])
  cv.samp.num <- floor(length(cv.ind)/nfolds)
  sampl_ind <- seq(1, dim(x)[1], by=1)
  candidate.model <- list()
  candidate.rocobj <- list()
  for (cv.id in 1:nfolds){
    if(cv.id < nfolds)
      idx <- cv.ind[seq((cv.id-1)*cv.samp.num + 1, cv.id*cv.samp.num, by=1)]
    else
      idx <- cv.ind[seq((cv.id-1)*cv.samp.num + 1, length(cv.ind), by=1)]
    
    train.set <- x[idx,]
    train.labels <- y[idx]
    
    validation.set <- x[!sampl_ind %in% idx,]
    validation.labels <- y[!sampl_ind %in% idx]
    
    candidate.model[[cv.ind]] <- ml.cross.validation(train.set, train.labels, nfolds, classifier, ...)
    candidate.rocobj[[cv.ind]] <- ml.roc(validation.set, candidate.model[cv.ind], validation.labels)
  }
  best.ind <- which.max(candidate.rocobj$auc) #### find the best auc index?
  best.model <- candidate.model[best.ind]
  
  if (savefile) save(best.model, file = paste(opts$outdir,"/trained.model", collapse='', sep=''))
  return(best.model)
}

