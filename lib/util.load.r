# This file is adapted from QIIME I/O module for R 
# (https://github.com/qiime/qiime.io/blob/master/R/I-methods.r)
# --------------
# Contributors: Dan, Gabe, PJ, Hu
# --------------
# Functions:
# reads a QIIME otu/metadata/taxon/distance table.
# Support legacy formats, where
# the header may or may not start with '#', 
# and comment lines can be anywhere in the file.
# return value is a matrix unless as.data.frame is TRUE
# --------------
# Last Update: 10/25/2014
#

if (!require("biom")) {
  install.packages("biom", dependencies = TRUE)
  library(biom)
}

# Read both BIOM format and classic OTU table
"read.qiime.table.mwas" <- function(filepath, as.data.frame=FALSE){
  #f <- file(filepath,'r')
  #if (grep(".biom",f)) {  # it's not a logical value, length is zero. Not working! 
  #if (file_ext(filepath) == "biom"){ # Need additonal pacakage: Tools
  
  if (tolower(tail(strsplit(filepath, "[.]")[[1]],n = 1)) == "biom"){ 
    # extract the file extension, if it's a BIOM format, then use BIOM package
    # require(biom, quietly=TRUE, warn.conflicts=FALSE)
    
    biom_table <- read_biom(filepath)          	   # OTU table - biom format
    datatable <- as.matrix(biom_data(biom_table))  # OTU table - classic format
    metadata <- observation_metadata(biom_table) 
    if(!is.null(metadata)){
      lineages <- apply(metadata, 1, function(x) paste(x, collapse=';')) # lineage names
    } else lineages <- NULL
  } else {
    # otherwise, the file could be a classic OTU table/mapping file etc.
    tryCatch(datatable <- read.qiime.classic.table(filepath, as.data.frame=as.data.frame), error = function(err) 
      print("If BIOM format, use .biom extension"))
    lineages <- NULL
  }
  
  if(!as.data.frame) datatable <- as.matrix(datatable)
  return(list(otus=datatable, lineages=lineages))
}

"read.qiime.classic.table" <- function(filepath, as.data.frame=FALSE){
  header.index <- get.header.index(filepath)
  # read the header
  f <- file(filepath,'r')
  header <- scan(filepath, what='character', sep='\t',comment='',skip=header.index-1,quote='"',
                 nlines=1,quiet=TRUE)
  close(f)
  # read the rest of the table
  datatable <- read.table(filepath,sep='\t',skip=header.index, comment='#',quote='"',
                          head=F,row.names=1,check=FALSE,strip.white=TRUE)
  
  # set column names using header
  colnames(datatable) <- header[-1]
  
  if(!as.data.frame) datatable <- as.matrix(datatable)
  return(datatable)
}


"load.qiime.mapping.file" <- function(filepath){
    return(read.qiime.classic.table(filepath, as.data.frame=TRUE))
}

"load.qiime.otu.table" <- function(filepath,include.lineages=FALSE){
    otus.obj <- read.qiime.table.mwas(filepath, as.data.frame=TRUE)
    otus <- otus.obj$otus
    lineages <-  otus.obj$lineages # from biom_read
    # drop "Consensus Lineage" column if present
    if(otu.table.has.metadata(colnames(otus))){
        C <- ncol(otus)
        lineages <- as.character(otus[,C])
        otus <- otus[,-C]
    }
    #} else {
    #    lineages <- NULL
    #}
    
    otus <- as.matrix(t(otus)) # need to check
    
    if(include.lineages){
        return(list(otus=otus,lineages=lineages))
    } else {
        return(otus=otus)
    }
}

# TRUE if last column is "Consensus Lineage" or "OTU Metadata"
"otu.table.has.metadata" <- function(headers){
    C <- length(headers)
    has.metadata <- grepl('consensus[ ]lineage|otu[ ]*metadata|taxonomy|kegg_pathways',
                          headers[C], ignore.case=TRUE)
    return(has.metadata)
}

# returns the index of the header line
# note: lines after the header may be comments with '#'
# read.table should therefore be called with (skip=header.index, comment='#')
"get.header.index" <- function(filepath){
    ncolumns.per.line <- NULL
    
    # read lines until the first line without a '#'
    # for each line, obtain the number of tab-delimited columns
    linecount <- 0
    start.character <- '#'
    while(start.character == '#'){
        linecount <- linecount + 1
        f <- file(filepath,'r') # open file in read mode
        line <- scan(f,what='character',skip=linecount-1,nlines=1, sep='\t',
                     quote='"', quiet=TRUE)
        close(f)
        # ncolumns is the number of entries in this line
        # not including trailing empties
        ncolumns <- max(which(sapply(line,nchar) > 0))
        ncolumns.per.line <- c(ncolumns.per.line, ncolumns)
        start.character <- substring(line[1],1,1)
    }
    
    # first non-comment line gives the number of columns
    C <- ncolumns.per.line[linecount]
    if(linecount == 1){
        # if there are no comment lines, then the first line is the header
        header.index <- 1
    } else {
        if(any(ncolumns.per.line[-linecount] == C)){
            # if there is a comment line with the correct number of columns,
            # it is the header
            header.index <- max(which(ncolumns.per.line[-linecount] == C))
        } else if(any(ncolumns.per.line[-linecount] == C+1)) {
            # if there is no comment line with the exacet number of columns, 
            # but one more column than the data tabel, then this line is considered
            # as the header. 
            # The last column is "Description" column, and all the data don't 
            # have any values for this entry column. 
           header.index <- max(which(ncolumns.per.line[-linecount] == C+1))
        } else {
            # if there is no comment line with the correct number of columns,
            # the first non-comment line is the header
            header.index <- linecount
        }
    }

    return(header.index)
}

# The same as load.qiime.otu.table(taxon_table_path)
#"load.qiime.taxon.table" <- function(filepath){
#    taxa <- as.matrix(t(read.table(filepath,sep='\t',head=T,row.names=1,check=FALSE,quote='"')))
#    return(taxa)
#}

"load.qiime.distance.matrix" <- function(filepath){
    d <- as.matrix(read.table(filepath,sep='\t',head=T,row.names=1,check=FALSE,quote='"'))
    return(d)
}

"load.qiime.pcoa.file" <- function(filepath){
	d <- as.matrix(read.table(filepath,sep='\t',head=T,row.names=1,check=FALSE,quote='"'))
	return(d)
}

"load.qiime.feature.stats" <- function(filepath){
	d < - as.matrix(read.table(filepath,sep='\t',head=T,row.names=1,check=FALSE,quote='"'))
	return(d)
}

# ensure map, data table, etc., contain the same samples in the same order
"remove.nonoverlapping.samples" <- function(map=NULL,otus=NULL,taxa=NULL,distmat=NULL){
    IDs <- NULL
    objects <- list(map=map,otus=otus,taxa=taxa,distmat=distmat)

    # find overlapping samples in all tables
    for(obj in objects){
        if(!is.null(obj)) {
            if(is.null(IDs)){
                IDs <- rownames(obj)
            } else {
                IDs <- intersect(rownames(obj), IDs)
            }
        }
    }
    
    # drop non-overlapping samples 
    for(i in 1:length(objects)){
        if(!is.null(objects[[i]])) {
            objects[[i]] <- objects[[i]][IDs,,drop=F]
            # for mapping file, drop any empty levels from factors that might
            # have occurred due to dropped samples
            if(i == 1) objects[[i]] <- droplevels(objects[[i]])
            # for distance matrix, get subset of columns too
            if(i == 4) objects[[i]] <- objects[[i]][,IDs]
        }
    }
    
    return(objects)
}

# returns vector of cluster ids for clusters with internal
# complete-linkage correlation of min.cor
"cluster.by.correlation" <- function(x, min.cor=.5){
  #     library('fastcluster')
  cc <- cor(x,use='pairwise.complete.obs',method='pear')
  if(ncol(x) == 379) browser()
  cc <- as.dist(1-cc)
  hc <- hclust(cc)
  res <- cutree(hc,h=1-min.cor)
  names(res) <- colnames(x)
  return(res)
}

# returns vector of cluster ids for clusters with internal
# complete-linkage correlation of min.cor
#
# by default, chooses cluster reps as highest-variance member
# if select.rep.fcn=mean
"collapse.by.correlation" <- function(x, min.cor=.5, select.rep.fcn=c('var','mean','lowest.mean',
                                                                      'longest.name', 'shortest.name')[2],
                                      verbose=FALSE){
  if(verbose) cat('Clustering',ncol(x),'features...')
  gr <- cluster.by.correlation(x, min.cor=min.cor)
  if(verbose) cat('getting means...')
  if(select.rep.fcn == 'mean'){
    v <- apply(x,2,function(xx) mean(xx,na.rm=TRUE))
  } else if(select.rep.fcn == 'lowest.mean'){
    v <- apply(x,2,function(xx) -mean(xx,na.rm=TRUE))
  } else if(select.rep.fcn == 'longest.name'){
    v <- nchar(colnames(x))
  } else if(select.rep.fcn == 'shortest.name'){
    v <- -nchar(colnames(x))
  } else {
    v <- apply(x,2,function(xx) var(xx,use='complete.obs'))
  }
  if(verbose) cat('choosing reps...')
  reps <- sapply(split(1:ncol(x),gr),function(xx) xx[which.max(v[xx])])
  if(verbose)
    cat(sprintf('collapsed from %d to %d.\n',ncol(x), length(reps)))
  return(list(reps=reps, groups=gr))
}


"inspect.env" <- function(){
  objs <- ls(1)
  objs <- objs[!is.null(objs)]
  res <- split(objs, sapply(objs, function(xx) class(get(xx))))
  for(i in seq_along(res)){
    cat('\n')
    cat(names(res)[i], ':\n',sep='')
    print(res[[i]])
  }
  invisible(res)
}


"shorten.taxonomy" <- function(ids,delim=';',num.levels=1, must.include.level=7){
  ids <- gsub('[kpcofgs]__$','Other',ids)
  ids <- gsub('[kpcofgs]__','',ids)
  newids <- ids
  ids <- strsplit(ids,delim)
  for(i in seq_along(ids)){
    n <- length(ids[[i]])
    j <- n
    while(ids[[i]][j] == 'Other' || ids[[i]][j] == '') j <- j - 1
    start.level <- min(must.include.level,j-num.levels+1)
    start.level <- max(1,start.level)
    newids[i] <- paste(ids[[i]][start.level:j],collapse=' ')
    if(j < n) newids[i] <- paste('Uncl. ',newids[i],sep='')
  }
  
  # add indices to duplicate names
  counts <- table(newids)
  if(max(counts) > 1){
    for(name in names(counts)){
      if(counts[name] > 1){
        ix <- which(newids == name)
        for(i in seq_along(ix)){
          newids[ix[i]] <- paste(name,i,sep=' ')
        }
      }
    }
  }
  
  return(newids)
}

# Get balanced folds where each fold has close to overall class ratio
#"balanced.folds" <- function(y, nfolds=10){
#  folds = rep(0, length(y))
#  y <- as.factor(y)
#  classes = levels(y)
#  # size of each class
#  Nk = table(y)
#  # -1 or nfolds = len(y) means leave-one-out
#  if (nfolds == -1 || nfolds == length(y)){
#    invisible(1:length(y))
#  } else{
#    # Can't have more folds than there are items per class
#    nfolds = min(nfolds, max(Nk))
#    # Assign folds evenly within each class, then shuffle within each class
#    for (k in 1:length(classes)){
#      ixs <- which(y==classes[k])
#      folds_k <- rep(1:nfolds, ceiling(length(ixs) / nfolds))
#      folds_k <- folds_k[1:length(ixs)]
#      folds_k <- sample(folds_k)
#      folds[ixs] = folds_k
#    }
#    invisible(folds)
#  }
#}

#"is.outlier" <- function(x,range=1.5){
#   # normalize first
#   x <- (x - mean(x)) / sd(x)
#   sx <- summary(x)
#   q1 <- sx[2]
#   q3 <- sx[5]
#   iqr <- q3 - q1
#   ret <- x > q3  + range * iqr | x < q1 - range * iqr
#   if(any(is.na(ret))) ret <- rep(FALSE,length(x))
#   if(all(ret)) ret <- !ret # don't call all outliers
#   return(ret)
# }

#"factor.to.numeric"<-function(x)
# {
#   x <- factor(x)
#   levels(x) <- 1:length(levels(x))
#   x <- as.numeric(x)
#   return(x)
# }

PARAMETERS <<- NULL
# used for parsing the params file. read into global var to avoid rereading too many times.
"parse.params" <- function(functionname, parameter)
{
	if(length(PARAMETERS)==0)
		PARAMETERS <- read.table("../config/params",sep=" ",row=1)
	ret <- levels(PARAMETERS[paste(functionname,parameter,sep=":"),])
	return(ret)
}

