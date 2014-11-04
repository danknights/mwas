# Function: import parameters for different functional modules
# -----------------
# Contributors: Hu
# -----------------
# input:
#       opts : options from the user
# --------
# output:
#   param.list : a table of parameters that needed in each corresponding fucntion
# 
# -------------
#  Last update: 10/27/2014
#

require(biom, quietly=TRUE, warn.conflicts=FALSE)

"import.train.params" <- function(opts){

  mapping <-  load.qiime.mapping.file(opts$map_fp)   # mapping file
  
  otus <- load.qiime.otu.table(opts$OTU_fp)  # OTU table - feature data for training
  
  samp.obj <- remove.nonoverlapping.samples(map=mapping,otus=otus)
  
  mapping <- samp.obj$map
  feat.Data <- samp.obj$otus
  
  response <- droplevels(factor(mapping[[opts$category]])) # desired labels 
  
  #print(dim(feat.Data))
  #print(response)
  param.list <- list(features=feat.Data, response=response, is.feat=opts$feat, method=opts$method, 
                     c.params=opts$param, valid_type=opts$validType, out.dir=opts$outdir)
  # c.params is parameter sets for the classifier
  class(param.list) <- "mwas"
  
  return(param.list)
}

"import.predict.params" <- function(opts){
  
  otus <- load.qiime.otu.table(opts$OTU_fp)  # OTU table
  
  if(!is.null(opts$map_fp)) {
    mapping <-  load.qiime.mapping.file(opts$map_fp)         # mapping file
    
    samp.obj <- remove.nonoverlapping.samples(map=mapping,otus=otus)
    
    mapping <- samp.obj$map
    otus <- samp.obj$otus
    
    response <- droplevels(factor(mapping[[opts$category]])) # desired labels 
  } else{
    mapping <- NULL
    response <- NULL
  }
  
  best.model <- readRDS(opts$param_fp)
  
  if("feat.set" %in% best.model) feat.Data <- otus[, best.model$feat.set]
  else feat.Data <- otus
  
  param.list <- list(features=feat.Data, trained.model=best.model$trained.model, response=response, 
                     out.dir=opts$outdir)
  class(param.list) <- "mwas"
  
  return(param.list)
}

"import.plot.params" <- function(opts){
  
  require('RColorBrewer', quietly=TRUE, warn.conflicts=FALSE)
  require('vegan', quietly=TRUE, warn.conflicts=FALSE)
  
  otu_table <- load.qiime.otu.table(opts$OTU_fp, include.lineages=TRUE)  # OTU table - feature data for training
  x <- otu_table$otus
 
  # differentiated feature table
  if(!is.null(opts$feat_stats)){ 
    feat_stats <- read.table(opts$feat_stats,sep='\t',head=T,row=1,check=F,quote='"',comment='')  
  }else feat_stats <- NULL 
  
  if (!is.null(otu_table$lineages)) kegg <- setNames(otu_table$lineages, rownames(x))
  
  if(!is.null(opts$map_fp)){
    m <-  load.qiime.mapping.file(opts$map_fp)         # mapping file
    # check rownames in mapping file matrix
    #missing.taxa.samples <- setdiff(rownames(x), rownames(m))
    #missing.map.samples <- setdiff(rownames(m), rownames(x))
    #if(length(missing.taxa.samples) > 0){
    
    #stop(sprintf('\n\nError: one or more sample names from taxonomy table (%s, ...) not present in metadata table (%s, ...).',
    #               paste(sort(missing.taxa.samples)[1:5],collapse=', '),
    #               paste(sort(missing.map.samples)[1:5],collapse=', ')))
    #}
    #x <- x[intersect(rownames(x),rownames(m)),,drop=F]
    #m <- droplevels(m[rownames(x),,drop=F])
    #m[rownames(x),,drop=F]
  }else m <- NULL
  samp.obj <- remove.nonoverlapping.samples(map=m,otus=x)
  
  m <- samp.obj$map
  x <- samp.obj$otus
  
  response <- droplevels(factor(m[[opts$category]])) # desired labels 
  #print(dim(x))
  #otus <- x
  #processed.obj <- preprocess.mwas(input.data=x, map=m, min_prevalence=opts$min_prevalence, transform_type=opts$transform_type)
  #x <- processed.obj$otu
  #kegg_pathways <- processed.obj$kegg_pathways
  kegg_pathways <- NULL
  # check that taxon.names are in taxon table
  if(is.null(opts$which_taxa)){
    taxon.names <- colnames(x)[rev(order(colMeans(x)))]
    taxon.names <- taxon.names[1:min(opts$nplot, length(taxon.names))]
  } else {
    taxon.names <- strsplit(opts$which_taxa,',')[[1]]
    #if(!all(taxon.names %in% colnames(x))){
    if(!all(sapply(taxon.names, function(xx) ifelse(length(grep(xx, colnames(x), value=F))>0, T, F)))){
      stop(paste('The following taxa are not present in the taxon table:',
                 paste(taxon.names[!(taxon.names %in% colnames(x))],collapse=', '),
                 '\n'))
    }
  }
  
  if(opts$shorten_taxa){
    colnames(x) <- shorten.taxonomy(colnames(x))
    taxon.names <- shorten.taxonomy(taxon.names)
  }
  
  if(is.null(opts$pcoa_fp)){
    if(is.null(opts$distance_fp)){
      d <- vegdist(x)
    } else {
      d <- read.table(opts$distance_fp,sep='\t',head=T,row=1,check=F)
      # check rownames in distance matrix
      missing.taxa.samples <- union(setdiff(rownames(x), rownames(d)), setdiff(rownames(x), colnames(d)))
      missing.distance.samples <- union(setdiff(rownames(d), rownames(x)), setdiff(colnames(d), rownames(x)))
      if(length(missing.taxa.samples) > 0){
        stop(sprintf('\n\nError: one or more sample names from taxonomy table (%s, ...) not present in distance table (%s, ...).',
                     paste(sort(missing.taxa.samples)[1:5],collapse=', '),
                     paste(sort(missing.distance.samples)[1:5],collapse=', ')))
      }
      d <- d[rownames(x),rownames(x)]
      d <- as.dist(d)
    }
    
    pc <- cmdscale(d,k=5)
  } else {
    pc <- read.table(opts$pcoa_fp,sep='\t',row=1,head=T)
    if(rownames(pc)[nrow(pc)] == '% variation explained'){
      pc <- pc[1:(nrow(pc)-2),1:min(5,ncol(pc))]
    }
    if(mean(rownames(x) %in% rownames(pc)) < 1){
      stop('Taxon table row names do not match PC file row names')
    }
    pc <- pc[rownames(x),]
  }
  
  if(is.null(opts$category)) {
    outdir <- sprintf('%s/gradients.pdf',opts$outdir)
    is.gradient = FALSE
  } else {
  #  if(!is.element(opts$column,colnames(m))) stop(paste(opts$column,'not in mapping file.'))
   # fp <- sprintf('%s/pcoa.pdf',opts$outdir)
  #  is.gradient = TRUE
  }
  
  param.list <- list(x=x, pc=pc, outdir=opts$outdir, m=m, 
                     kegg_pathways = kegg_pathways,
                     taxon.names = taxon.names, 
                     category = opts$category,
                     response = response,
                     is.multiple_axes = opts$multiple_axes,
                     category_order = opts$category_order,
                     is.sort_by_abundance = opts$sort_by_abundance, 
                     num_taxa = opts$nplot,
                     alpha = opts$alpha,
                     x_axis_label = opts$x_axis_label,
                     out.dir = opts$outdir,
                     plot.type = opts$plottype,
                     min_prevalence = opts$min_prevalence, 
                     transform_type = opts$transform_type,
                     suppress_relative_abundance_conversion = opts$suppress_relative_abundance_conversion,
                     #kegg = kegg,
                     feat_stats = feat_stats)
  class(param.list) <- "mwas"
  
  return(param.list)
}

"import.stats.params" <- function(opts){

  mapping <-  load.qiime.mapping.file(opts$map_fp)   # mapping file
  
  feat.Data <- load.qiime.otu.table(opts$OTU_fp)  # OTU table - feature data for training
  
  response <- droplevels(factor(mapping[, opts$category])) # desired labels 
  
  param.list <- list(features=feat.Data, response=response, alpha=opts$alpha, is.parametric=opts$parametric, 
                     include.subset=opts$subset)
  class(param.list) <- "mwas"
  
  return(param.list)
}