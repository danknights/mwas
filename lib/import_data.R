# Function: import parameters for different functional modules
# -----------------
# Contributor: Hu
# -----------------
# input:
#       opts : options from the user
#              parameters are slightly different for the four modes:
#              train -- input_fp
#                       category
#                       map_fp
#                       is_feat (is.feat_select)
#                       feat_param (parameters for feature selection)
#                       method
#                       method_param
#                       suppress_relative_abundance_conversion (optional)
#                       min_prevalence
#                       transform_type
#                       collapse_table
#                       validType (validation type)
#                       nfolds (# of folds in cross validation)
#                       outdir
#            predict -- input_fp
#                       map_fp (optional)
#                       category (optional, must be given with map_fp)
#                       method (trained model object file)
#                       suppress_relative_abundance_conversion (optional)
#                       min_prevalence
#                       transform_type
#                       collapse_table
#                       outdir
#               plot -- input_fp
#                       map_fp
#                       category
#                       method (plot type)
#                       feat_stats (optional)
#                       pcoa_fp (optional)
#                       distance_fp (optional)
#                       suppress_relative_abundance_conversion (optional)
#                       min_prevalence
#                       transform_type
#                       collapse_table
#                       which_taxa
#                       filter_kegg (optional, for heatmap)
#                       outdir
#         statistics -- 
#
# --------
# output:
#   param.list : a table of parameters that needed in each corresponding fucntion
# 
# -------------
#  Last update: 11/12/2014
#

require(biom, quietly=TRUE, warn.conflicts=FALSE)

"import.train.params" <- function(opts){

  mapping <-  load.qiime.mapping.file(opts$map_fp)   # mapping file
  
  otu <- load.qiime.otu.table(opts$input_fp, include.lineages=FALSE)
  
  # preprocessing
  if(opts$suppress_relative_abundance_conversion) {
    is.relative.conversion = FALSE
  } else {
    if (sum(rowSums(otu)) != dim(otu)[1]) is.relative.conversion = TRUE
    else is.relative.conversion = FALSE
  }
  preporcessed.obj <- preprocess.mwas(input.data = otu, 
                                      map = mapping, 
                                      min_prevalence = opts$min_prevalence,
                                      transform_type = opts$transform_type,
                                      is.collapse = opts$collapse_table,
                                      is.relative.conversion=is.relative.conversion
                                      )
  feat.Data <- preporcessed.obj$otu
  mapping <- preporcessed.obj$map
  
  response <- droplevels(factor(mapping[[opts$category]])) # desired labels 
  
  param.list <- list(features=feat.Data, 
                     response=response, 
                     is.feat=opts$is_feat, 
                     feat.param = opts$feat_param,
                     method=opts$method, 
                     c.params=opts$method_param, 
                     valid_type=opts$validType, 
                     out.dir=opts$outdir, 
                     nfolds=opts$nfolds)
  # c.params is parameter sets for the classifier
  class(param.list) <- "mwas"
  
  return(param.list)
}

"import.predict.params" <- function(opts){
  
  otu <- load.qiime.otu.table(opts$input_fp, include.lineages=FALSE)
  
  if(!is.null(opts$map_fp)) {
    mapping <-  load.qiime.mapping.file(opts$map_fp)         # mapping file
    
    response <- droplevels(factor(mapping[[opts$category]])) # desired labels 
  } else{
    mapping <- NULL
    response <- NULL
  }

  # preprocessing
  if(opts$suppress_relative_abundance_conversion) {
    is.relative.conversion = FALSE
  } else {
    if (sum(rowSums(otu)) != dim(otu)[1]) is.relative.conversion = TRUE
    else is.relative.conversion = FALSE
  }
  
  preporcessed.obj <- preprocess.mwas(input.data = otu, 
                                      map = mapping, 
                                      min_prevalence = opts$min_prevalence,
                                      transform_type = opts$transform_type,
                                      is.collapse = opts$collapse_table,
                                      is.relative.conversion=is.relative.conversion
                                      )
  otus <- preporcessed.obj$otu
  mapping <- preporcessed.obj$map
  
  best.model <- readRDS(opts$method)
  
  #if("feat.set" %in% best.model) {
  if(!is.null(best.model$feat.set)) {
      feat.Data <- otus[, best.model$feat.set]
  } else feat.Data <- otus
  
  param.list <- list(features=feat.Data, 
                     trained.model=best.model$trained.model, 
                     response=response, 
                     out.dir=opts$outdir)
  class(param.list) <- "mwas"
  
  return(param.list)
}

"import.plot.params" <- function(opts){
  
  require('RColorBrewer', quietly=TRUE, warn.conflicts=FALSE)
  require('vegan', quietly=TRUE, warn.conflicts=FALSE)
  
  # load OTU table/taxon table in either BIOM or txt format 
  otu_table <- load.qiime.otu.table(opts$input_fp, include.lineages=TRUE)  # OTU table - feature data for training
  otu <- otu_table$otus
  kegg <- otu_table$lineages
  if (!is.null(kegg)&&is.null(names(kegg))) { # if the lineage doesn't have OTU IDs
    kegg <- setNames(otu_table$lineages, colnames(otu))
  }

  # load differentiated feature table
  if(!is.null(opts$feat_stats_fp)){ 
    feat_stats <- read.table(opts$feat_stats_fp,sep='\t',head=T,row=1,check=F,quote='"',comment='')  
  }else feat_stats <- NULL 
  
  # load mapping file
  if(!is.null(opts$map_fp)){
    m <-  load.qiime.mapping.file(opts$map_fp)         # mapping file
  }else m <- NULL
  
  if(is.null(opts$pcoa_fp)){
    if(is.null(opts$distance_fp)){
      d <- as.matrix(vegdist(otu))
    } else {
      d <- load.qiime.distance.matrix(opts$distance_fp)
    }
    
    pc <- cmdscale(d,k=5)
  } else {
    pc <- read.table(opts$pcoa_fp,sep='\t',row=1,head=T)
    #d <- load.qiime.pcoa.file(opts$distance_fp)
    if(rownames(pc)[nrow(pc)] == '% variation explained'){
      pc <- pc[1:(nrow(pc)-2),1:min(5,ncol(pc))]
    }
    if(mean(rownames(otu) %in% rownames(pc)) < 1){
      stop('Taxon table row names do not match PC file row names')
    }
    pc <- pc[rownames(otu),]
  }
  
  # whether supress converting to relative abundance
  if(opts$suppress_relative_abundance_conversion) {
    is.relative.conversion = FALSE
  } else {
    if (sum(rowSums(otu)) != dim(otu)[1]) is.relative.conversion = TRUE
    else is.relative.conversion = FALSE
  }
  
  preporcessed.obj <- preprocess.mwas(input.data = otu, 
                                      map = m, 
                                      distMat = d,
                                      kegg = kegg, 
                                      min_prevalence = opts$min_prevalence,
                                      transform_type = opts$transform_type,
                                      is.filter.kegg = opts$filter_kegg,
                                      is.collapse = opts$collapse_table,
                                      is.relative.conversion=is.relative.conversion
                                      )
  
  otu=preporcessed.obj$otu
  kegg_pathways=preporcessed.obj$kegg_pathways 
  m = preporcessed.obj$map
  distMat = preporcessed.obj$distMat
  
  response <- droplevels(factor(m[[opts$category]])) # desired labels 
  names(response) <- rownames(m)
  #print(dim(x))
  #otus <- x
  #processed.obj <- preprocess.mwas(input.data=x, map=m, min_prevalence=opts$min_prevalence, transform_type=opts$transform_type)
  #x <- processed.obj$otu
  #kegg_pathways <- processed.obj$kegg_pathways
  #kegg_pathways <- NULL
  
  # check that taxon.names are in taxon table
  if(is.null(opts$which_taxa)){
    taxon.names <- colnames(otu)[rev(order(colMeans(otu)))]
    taxon.names <- taxon.names[1:min(opts$nplot, length(taxon.names))]
  } else {
    taxon.names <- strsplit(opts$which_taxa,',')[[1]]
    if(!all(taxon.names %in% colnames(otu))){
    #if(!all(sapply(taxon.names, function(xx) ifelse(length(grep(xx, colnames(x), value=F))>0, T, F)))){
      stop(paste('The following taxa are not present in the taxon table:',
                 paste(taxon.names[!(taxon.names %in% colnames(otu))],collapse=', '),
                 '\n'))
    }
  }
  
  if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)
  
  param.list <- list(otu=otu, 
                     m=m, 
                     out.dir=opts$outdir,
                     pc=pc, 
                     kegg_pathways = kegg_pathways,
                     is.shorten.taxa = opts$shorten_taxa, 
                     taxon.names = taxon.names,
                     category = opts$category,
                     response = response,
                     is.multiple_axes = opts$multiple_axes,
                     alpha = opts$alpha,
                     plot.type = opts$method,
                     feat_stats = feat_stats)
  class(param.list) <- "mwas"
  
  return(param.list)
}

"import.stats.params" <- function(opts){

  mapping <-  load.qiime.mapping.file(opts$map_fp)   # mapping file
  
  feat.Data <- load.qiime.otu.table(opts$input_fp)  # OTU table - feature data for training
  
  response <- droplevels(factor(mapping[, opts$category])) # desired labels 
  
  param.list <- list(features=feat.Data, response=response, alpha=opts$alpha, is.parametric=opts$parametric, 
                     include.subset=opts$subset)
  class(param.list) <- "mwas"
  
  return(param.list)
}