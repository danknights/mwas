# Function: import parameters 
#
# Contributors: Hu
# -----
# input:
#       opts : options from the user
# 
# output:
#   param.list : a table of parameters that needed in each corresponding fucntion
# 
# ------
#  Last update: 10/27/2014
#

"import.train.params" <- function(opts){
  mapping <-  load.qiime.mapping.file(opts$map_fp)   # mapping file
  
  if (grep(".biom$",opts$OTU_fp)) {
    biom_table <- read_biom(opts$OTU_fp)         # OTU table - biom format
    otus <- t(as.matrix(biom_data(biom_table)))      # OTU table - classic format
  } else {
    tryCatch(otus <- read.delim(opts$OTU_fp, sep='\t',
                                comment='',head=T,row.names=1,check.names=F),error = function(err) 
                                  print("Couldn't parse OTU table. If BIOM format, use .biom extension"))
  }
  
  feat.Data <- otus # feature data for training
  
  response <- droplevels(factor(data.list$map[[opts$category]])) # desired labels 
  
  #print(dim(feat.Data))
  #print(response)
  param.list <- list(features=feat.Data, response=response, is.feat=opts$feat, method=opts$method, 
                     c.params=opts$param, valid_type=opts$validType)
  class(param.list) <- "mwas"
  
  return(param.list)
}

"import.predict.params" <- function(opts){
  if(!is.null(opts$map_fp) {
    mapping <-  load.qiime.mapping.file(opts$map_fp)         # mapping file
    response <- droplevels(factor(mapping[[opts$category]])) # desired labels 
  } else{
    mapping <- NULL
    response <- NULL
  }
 
  if (grep(".biom$",opts$OTU_fp)) {
    biom_table <- read_biom(opts$OTU_fp)           # OTU table - biom format
    otus <- t(as.matrix(biom_data(biom_table)))      # OTU table - classic format
  } else {
    tryCatch(otus <- read.delim(opts$OTU_fp, sep='\t',
                                comment='',head=T,row.names=1,check.names=F),error = function(err) 
                                  print("Couldn't parse OTU table. If BIOM format, use .biom extension"))
  }
  
  best.model <- readRDS(opts$param_fp)
  
  if("feat.set" %in% best.model) feat.Data <- otus[, best.model$feat.set]
  else feat.Data <- otus
  
  
  param.list <- list(features=feat.Data, trained.model=best.model$trained.model, response=response)
  class(param.list) <- "mwas"
  
  return(param.list)
}

"import.plot.params" <- function(opts){
  mapping <-  load.qiime.mapping.file(opts$map_fp)         # mapping file
  dist.matrix <- load.qiime.distance.matrix(opts$param_fp) # distance matrix
  
  if (grep(".biom$",opts$OTU_fp)) {
    biom_table <- read_biom(opts$OTU_fp)         # OTU table - biom format
    otus <- t(as.matrix(biom_data(biom_table)))      # OTU table - classic format
  } else {
    tryCatch(otus <- read.delim(opts$OTU_fp, sep='\t',
                                comment='',head=T,row.names=1,check.names=F),error = function(err) 
                                  print("Couldn't parse OTU table. If BIOM format, use .biom extension"))
  }
  
  feat.Data <- otus # feature data for training
  
  response <- droplevels(factor(data.list$map[[opts$category]])) # desired labels 
  
  #print(dim(feat.Data))
  #print(response)
  param.list <- list(features=feat.Data, response=response, dist.matrix=opts$feat, method=opts$method, 
                     c.params=opts$param, valid_type=opts$validType)
  class(param.list) <- "mwas"
  
  return(param.list)
}

"import.stats.params" <- function(opts){
  
  if (grep(".biom$",opts$OTU_fp)) {
    biom_table <- read_biom(opts$OTU_fp)         # OTU table - biom format
    otus <- t(as.matrix(biom_data(biom_table)))      # OTU table - classic format
  } else {
    tryCatch(otus <- read.delim(opts$OTU_fp, sep='\t',
                                comment='',head=T,row.names=1,check.names=F),error = function(err) 
                                  print("Couldn't parse OTU table. If BIOM format, use .biom extension"))
  }
  
  feat.Data <- otus # feature data for training
  
  param.list <- list(features=feat.Data, response=response, is.feat=opts$feat, method=opts$method, 
                     c.params=opts$param, valid_type=opts$validType)
  class(param.list) <- "mwas"
  
  return(param.list)
}