
# Visualization function: including heatmaps, beeswarm, violin, boxplots, etc
# 
# ------------
# Contributors: Dan, PJ, Emmanuel, Gabe, Hu
# ------------
# Input:
#
# -------
# Output:
#
# -------
# Last update: 11/04/2014
#
# Need to fix heatmap plots

#require(beeswarm, quietly=TRUE, warn.conflicts=FALSE)

"plot.mwas" <- function(data, ...){ 
  
  options <- list(...)

  if (class(data) == "mwas"){
    x <- data$otu 
    pc <- data$pc 
    out.dir <- data$out.dir
    m <- data$m 
    response <- data$response
    taxon.names <- data$taxon.names
    is.shorten.taxa <- data$is.shorten.taxa
    category <- data$category
    is.multiple_axes <- data$is.multiple_axes
    alpha <- data$alpha
    plot.type <- data$plot.type
    feat_stats <- data$feat_stats
    nplot <- opts$nplot
  } else{
    x <- data 
    pc <- options$pc 
    out.dir <- options$out.dir
    m <- options$m 
    response <- options$response
    taxon.names <- options$taxon.names
    is.shorten.taxa <- options$is.shorten.taxa
    category <- options$category
    is.multiple_axes <- options$is.multiple_axes
    alpha <- options$alpha
    plot.type <- options$plot.type
    feat_stats <- options$feat_stats
  }
  #print(plot.type)
  switch(plot.type,
         beeswarm = {
           if (!is.null(alpha)) {
             # if alpha is given, then use the differentiated test feature table
             # otherwise use the whole taxa in the original table.
             
             diff.obj <- diff.plot.parameters(x, response,  alpha, feat_stats, nplot, is.shorten.taxa)
             run.beeswarm(diff.obj$new_taxon_table, 
                          diff.obj$response, 
                          diff.obj$filename,
                          out.dir = out.dir)
           }else  { # plot all the taxa that are provided in the table
             #print("Beeswarm plot...")
             run.beeswarm(x = x, response=response, filename="beeswarm-plot.pdf", out.dir = out.dir)
           }
         },
         gradient = {
           #processed.data <- preprocess.mwas(data)
           plot.gradients(x=x, 
                          pc=pc, out.dir=out.dir, m=m,
                          taxon.names=taxon.names, 
                          category=category,
                          is.multiple_axes=is.multiple_axes,
                          is.shorten.taxa = is.shorten.taxa)
         },
         heatmap = { # need to fix
           heatmap.mwas(x, map, diff.features, cluster.var=c("Sex", "Treatment"), 
                        color.var=names(color.list), color.list, 
                        kegg_pathways=kegg_pathways, 
                        heatmap.title=heatmap.title, 
                        outputfile=fp)
         }, 
         scatterplot = {
           if (!is.null(alpha)) {
             # if alpha is given, then use the differentiated test feature table
             # otherwise use the whole taxa in the original table.
             
             diff.obj <- diff.plot.parameters(x, response,  alpha, feat_stats, is.shorten.taxa)
             run.2d.scatterplot(diff.obj$new_taxon_table, 
                                diff.obj$response, 
                                diff.obj$filename,
                                out.dir = out.dir)             
           }else  {
             print("2D scatter plot...")
             run.2d.scatterplot(x=x, response=response, out.dir = out.dir)
           }
         },
         stop("Please assign the correct plot type!(Optioins: beeswarm, graidents, heatmap, scatterplot.")
    )
}

# This function is replaced by "run.beeswarm" function
#"plot.beeswarm.dt" <- function(cols, x_axis_label, hit.ix, env, outdir) {
#  require(beeswarm, quietly=TRUE, warn.conflicts=FALSE)
#  
#  cols <-  c(brewer.pal(9,'Set1'),brewer.pal(9,'Pastel1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[-6]
#  cols[1:2] <- cols[2:1]
#  cols <- sprintf('%sbb',cols)
#  if (!length(hit.ix)) hit.ix = seq(1,dim(x)[2])
#  #print(x[,1])
#  #print(env)
#  #pdf(filename.i,width=8.5,height=11)
#  #par(oma=rep(c(0,0),2),mar=c(12,8,8,3), cex.axis=.65, cex.lab=.65, cex=2)
#  for(i in hit.ix){
#    taxon <- x[,i]
#    taxon.name <- colnames(x)[i]
#    filename.i <- sprintf('%s/beeswarm-%s.pdf',outdir, taxon.name)
#    pdf(filename.i,width=8.5,height=11)
#    par(oma=rep(c(0,0),2),mar=c(12,8,8,3), cex.axis=.65, cex.lab=.65, cex=2)
#   cex.lab <- max(.45, 1 - strwidth(taxon.name, units='figure') / 4)
#    cex.lab <- .45
#   cex.axis <- .45
#    beeswarm(taxon ~ env,corral='random',las=2,
#             col='#000000bb',
#             bg=cols,
#             pch=21,ylab=taxon.name,xlab=x_axis_label,cex.lab=cex.lab, cex.axis=cex.axis)
#    bxplot(taxon ~ env,add=TRUE)
#    dev.off()
#  }
#}


# Beeswarm of several conditions included in the mapping file
#
# Contributors: Emmanuel, Hu
# ------
# input:
#        x:  OTU/taxon table
# response:  response lable
# ------
# output:
#     save the plot as PDF file
# ------
# Last update: 11/03/2014
#
"run.beeswarm" <- function(x, response, filename="beeswarm-plot.pdf", out.dir){
  
  x2beeswarmfile <- cbind(response,x)
  beeswarmfile <- x2beeswarmfile[order(response,decreasing=TRUE),]
  beeswarmfiletitle <- colnames(beeswarmfile)
  beeswarmfile2 <-beeswarmfile
  beeswarmfile3 <- beeswarmfile2[, -1]
  
  cols <-  c(brewer.pal(9,'Set1'),brewer.pal(9,'Pastel1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[-6]
  cols[1:2] <- cols[2:1]
  cols <- sprintf('%sbb',cols)
  
  if(!is.null(out.dir)) {
    file.out <- sprintf('%s/%s', out.dir, filename)
    pdf(file.out,width=6,height=7)
    #print(file.out)
  }
  
  for (i in 1:dim(beeswarmfile3)[2]){
    #print(i)
    par(mar =c(13,4,3,2))
    beeswarm(beeswarmfile3[,i] ~ response, data = beeswarmfile2,
             pch = 21,
             corral='random',las = 2,
             col='#000000bb',
             bg=cols,
             main = beeswarmfiletitle[i+1],
             #labels= "",
             labels=levels(response),
             xlab="", ylab ="")
    bxplot(beeswarmfile3[,i] ~ response,add=TRUE)
    #text(x = 1:3, y=-0.05, labels=levels(response), srt = 20, pos=1, xpd=TRUE)
  }
  if(!is.null(out.dir)) dev.off()
  return("Please check the PDF file.")
}

# Scatterplot of 2 conditions included in the mapping file
# and perform spearman correlation for each taxon
#
# Contributors: Emmanuel, Hu
# -------
# Input:
#    m:  mapping file
#   x: otu table
# -------
# Output:
#   save plot as a PDF file
# -------
# Last update: 11/04/2014
#

"run.2d.scatterplot" <- function(x, response, filename="scatter-plot.pdf", out.dir){

  level.names = levels(response) 
  level.num = length(level.names)
  for (id in 1:(level.num-1)){
    for (jd in (id+1):level.num){
      
      if(!is.null(out.dir)) {
        file.out <- sprintf('%s/categories-%d-%d-%s', out.dir, id, jd, filename)
        pdf(file.out,width=4,height=4)
        #print(file.out)
      }
      print(id)
      print(jd)
      x1 <- as.data.frame(x[names(response[response==level.names[id]]),])
      x2 <- as.data.frame(x[names(response[response==level.names[jd]]),])
      
      if (dim(x1)[1] > dim(x2)[1]) {
        rand.id <- sample(dim(x1)[1], size = dim(x2)[1])
        xx1 <- x1[rand.id, ]
        xx2 <- x2
      } else{
        rand.id <- sample(dim(x2)[1], size = dim(x1)[1])
        xx1 <- x1
        xx2 <- x2[rand.id, ]
      }
      
      for (i in 1:dim(x1)[2]){
        if ((max(x1[,i])<0.1) & (max(x2[,i])< 0.1)){
          
          plot(xx1[,i],xx2[,i],xlim=c(0,0.1),ylim=c(0,0.1), col=c("firebrick2"), pch= 20, 
               xlab=level.names[id],ylab=level.names[jd],main = colnames(xx1)[i], cex.main=1.5)
          abline(0, 1)
          
        } else if ((max(x1[,i])<0.3) & (max(x2[,i])< 0.3)){
          
          plot(xx1[,i],xx2[,i],xlim=c(0,0.4),ylim=c(0,0.4), col=c("firebrick2"), pch= 20, 
               xlab=level.names[id],ylab=level.names[jd],main= colnames(xx1)[i],cex.main=1.5)
          abline(0, 1)
          
        } else {
          plot(xx1[,i],xx2[,i],xlim=c(0,0.6),ylim=c(0,0.9), col=c("firebrick2"), pch= 20, 
               xlab=level.names[id],ylab=level.names[jd],main= colnames(xx1)[i],cex.main=1.5)
          abline(0, 1)
        }
        legend("topleft", sprintf('Spearman Coeff.: %.4f', cor(xx1[,i], xx2[,i], method="spearman")),
               pt.cex=0.1, adj=0.05, cex=0.85)
      }
      if(!is.null(out.dir)) dev.off()
    }
  }
 
}

"plot.gradients" <- function(x, pc, out.dir, m=NULL, taxon.names=NULL, category=NULL,
                             is.multiple_axes=FALSE, is.shorten.taxa = TRUE){
  
  if (is.null(category)){ 
    # if the category is not specified, then plot gradient
    # else show the metadata.
    filename <- "gradient-plot.pdf"
  } else filename <- sprintf("%s-plot.pdf", category)
  
  if(!is.null(out.dir)) {
    file.out <- sprintf('%s/%s', out.dir, filename)
  }else file.out <- filename
  
  if(is.multiple_axes){
    pdf(file.out,width=11,height=3.75)
    par(mfrow=c(1,3))
    combs <- combn(1:3,2)
  } else {
    pdf(file.out,width=6,height=5)
    combs <- matrix(1:2,ncol=1)
  }
  
  #  if(is.null(opts$category)) {
  #    outdir <- sprintf('%s/gradients.pdf',opts$outdir)
  #    is.gradient = FALSE
  #  } else {
  #  if(!is.element(opts$column,colnames(m))) stop(paste(opts$column,'not in mapping file.'))
  # fp <- sprintf('%s/pcoa.pdf',opts$outdir)
  #  is.gradient = TRUE
  #  }
  
  if(is.null(category)){
    for(i in seq_along(taxon.names)){
      for(j in 1:ncol(combs)){
        show.gradients(x[,taxon.names[i]], pc[,combs[,j]], incl.legend=TRUE,pt.alpha='CC',
                       axis.labels=sprintf('PC%d',combs[,j]),
                       title.text=sprintf('%s - PC%d v PC%d',taxon.names[i],combs[1,j],combs[2,j]))
      }
    }
  } else {
    for(j in 1:ncol(combs)){
      show.metadata(m[,category], pc[,combs[,j]], incl.legend=TRUE,pt.alpha='CC',
                    axis.labels=sprintf('PC%d',combs[,j]),
                    title.text=sprintf('%s - PC%d v PC%d',category,combs[1,j],combs[2,j]))
    }
  }
  if(!is.null(out.dir)) dev.off()
  return("Please check the PDF file.")
}

"diff.plot.parameters" <- function(x, response,  alpha, feat_stats=NULL, nplot=50, is.shorten.taxa=TRUE){
  
  if (!is.null(feat_stats)){
    # if the feature statistic table is provided, then load from file
    # else calculate q-values
    ft.qvalue <- subset(feat_stats, qvalues <= alpha) # ft - differentiated feature vector
    ft.qvalue <- t(ft.qvalue)
    ft.qvalue <- ft.qvalue[,order(ft.qvalue["pvalues", ], decreasing=F)[1:min(nplot,dim(ft.qvalue)[2])]]
    
    keep_bugs <- colnames(x)[colnames(x) %in% colnames(ft.qvalue)]
  } else { 
    
    # if the feature stats table is not given, 
    diff.table <- differentiation.test(x, response, alpha = alpha)
    hit.ix <- order(diff.table$pvalues, decreasing=F)[1:min(nplot,length(diff.table$pvalues))]
    ft.qvalue <- diff.table$qvalues[hit.ix]
    
    keep_bugs <- colnames(x)[colnames(x) %in% names(ft.qvalue)]
  }
  
  if (length(keep_bugs)==0) stop("Please input a taxon table (not OTU table); or there is no overlapped taxa information between these two tables.")
  new_taxon_table <- x[, keep_bugs]
  
  # shorten taxonomy name if specified
  if(is.shorten.taxa) {
    if (!is.null(feat_stats)) colnames(ft.qvalue) <- shorten.taxonomy(colnames(ft.qvalue))
    else names(ft.qvalue) <- shorten.taxonomy(names(ft.qvalue))
    
    colnames(new_taxon_table) <- shorten.taxonomy(colnames(new_taxon_table))
  }
  filename <- sprintf("beeswarm-plot-alpha-%.3f.pdf", alpha)
  
  return(list(new_taxon_table=new_taxon_table, response=response, filename=filename))
}

#QQ plot
# http://GettingGeneticsDone.blogspot.com/
# See http://gettinggeneticsdone.blogspot.com/p/copyright.html

# Define the function
#"qqplot.pvals"   <- function(pvector, main=NULL, ...) {
#  o = -log10(sort(pvector,decreasing=F))
#  e = -log10( 1:length(o)/length(o) )
#  plot(e,o,pch=19,cex=1, main=main, ...,
#       xlab=expression(Expected~~-log[10](italic(p))),
#       ylab=expression(Observed~~-log[10](italic(p))),
#       xlim=c(0,max(e)), ylim=c(0,max(o)))
#  lines(e,e,col="red")
#}