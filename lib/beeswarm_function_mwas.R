#beeswarm of several conditions included in the mapping file (here Treatment: preCT, postCT)
#require m:  mapping file 
# x: otu table 

source('~/DKpostdoc/src/mwas/lib/util.load.r')


#do beeswarm 
"run.beeswarm" <- function(m,x,filename='beeswarm_pre_postCT.pdf'){
  x <- x[order(row.names(x)),]
  m <- m[order(row.names(m)),]
  x2beeswarmfile <- cbind(m,x)
  x2beeswarmfile <- x2beeswarmfile[order(x2beeswarmfile$Treatment,decreasing = FALSE),]
  beeswarmfile <- x2beeswarmfile
  attributes(beeswarmfile)
  beeswarmfiletitle <- names(beeswarmfile)
  beeswarmfiletitle <- shorten.taxonomy(beeswarmfiletitle)
  colnames(beeswarmfile) <- c(beeswarmfiletitle)
  category <- beeswarmfile$Treatment
  beeswarmfile2 <-beeswarmfile[,-c(1,2,4,5)]
  beeswarmfile3 <- asin(sqrt(beeswarmfile2[-1]))
 
  require(beeswarm)
  
  if(!is.null(filename)) pdf(filename,width=4,height=4)
  for (i in 1:length(beeswarmfile3)){
    beeswarm(beeswarmfile3[,i] ~ category, data = beeswarmfile2,
             pch = 16, 
             col = rainbow(8),
             labels = c("postCT", "preCT"),
             main = colnames(beeswarmfile3)[i],
             xlab = "", ylab ="")
    
  }
  if(!is.null(filename)) dev.off()
}  









