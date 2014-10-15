#beeswarm of several conditions included in the mapping file (here Treatment: preCT, postCT)
#require m:  mapping file 
# x: otu table 

x <- x[order(row.names(x)),]
m <- m[order(row.names(m)),]
x2beeswarmfile <- cbind(m,x)

x2beeswarmfile <- x2beeswarmfile[order(x2beeswarmfile$Treatment),]
beeswarmfile <- x2beeswarmfile

attributes(beeswarmfile)
beeswarmfiletitle <- names(beeswarmfile)

"shorten.taxonomy" <- function(beeswarmfiletitle,delim=';'){
  beeswarmfiletitle <- gsub('[kpcofgs]__','',beeswarmfiletitle)
  newbeeswarmfiletitle <- beeswarmfiletitle
  beeswarmfiletitle <- strsplit(beeswarmfiletitle,delim)
  for(k in seq_along(beeswarmfiletitle)){
    n <- length(beeswarmfiletitle[[k]])
    j <- n
    while(beeswarmfiletitle[[k]][j] == 'Other' || beeswarmfiletitle[[k]][j] == '') j <- j - 1
    newbeeswarmfiletitle[k] <- beeswarmfiletitle[[k]][j]
  }
  return(newbeeswarmfiletitle)
}

beeswarmfiletitle <- shorten.taxonomy(beeswarmfiletitle)
colnames(beeswarmfile) <- c(beeswarmfiletitle)
attributes(beeswarmfile)
beeswarmfiletitle <- names(beeswarmfile)

CT <- beeswarmfile$Treatment

beeswarmfile2 <-beeswarmfile[,-c(1,2,4,5)]
beeswarmfile3 <- asin(sqrt(beeswarmfile2[-1]))
attach(beeswarmfile3)
attributes(beeswarmfile3)
tax <- names(beeswarmfile3)

#do beeswarm 
"run.beeswarm" <- function(beeswarmfile2, filename='beeswarm_pre_postCT.pdf'){
  require(beeswarm)
  
  if(!is.null(filename)) pdf(filename,width=4,height=4)
  for (i in beeswarmfile3){
    
    main.title = tax[i]
    
    beeswarm(i ~ CT, data = beeswarmfile2, 
             pch = 16, 
             col = rainbow(8),
             labels = c("preCT", "postCT"),
             main = main.title,
             xlab = "", ylab ="")
    
  }
  if(!is.null(filename)) dev.off()
}  









