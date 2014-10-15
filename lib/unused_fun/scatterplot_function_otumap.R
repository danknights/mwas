#scatterplot of 2 conditions included in the mapping file (here Treatment preCT, postCT)


#require m:  mapping file 
# x: otu table 

#order following names of the patients in mapping file and otu table
x <- x[order(row.names(x)),]
m <- m[order(row.names(m)),]

x2scatterplot <- cbind(m,x)

x2scatterplot<- x2[order(x2$Treatment),]
scatterplotfile <- x2scatterplot

scatterplotfile2 <-scatterplotfile[,-c(1,2,3,4,5)]
attributes(scatterplotfile2)
scatterplotitle2 <- names(scatterplotfile2)


"shorten.taxonomy" <- function(scatterplotitle2,delim=';'){
  scatterplotitle2 <- gsub('[kpcofgs]__','',scatterplotitle2)
  newscatterplotitle2 <- scatterplotitle2
  scatterplotitle2 <- strsplit(scatterplotitle2,delim)
  for(k in seq_along(scatterplotitle2)){
    n <- length(scatterplotitle2[[k]])
    j <- n
    while(scatterplotitle2[[k]][j] == 'Other' || scatterplotitle2[[k]][j] == '') j <- j - 1
    newscatterplotitle2[k] <- scatterplotitle2[[k]][j]
  }
  return(newscatterplotitle2)
}

scatterplotitle2 <- shorten.taxonomy(scatterplotitle2)
colnames(scatterplotfile2) <- c(scatterplotitle2)


after <- asin(sqrt(scatterplotfile2[1:15,]))
before <- asin(sqrt(scatterplotfile2[16:30,]))
attach(scatterplotfile2)

attributes(scatterplotfile2)
scatterplotitle <- names(scatterplotfile2)


"run.scatterplot" <- function(scatterplotfile2, filename='scatterplot_pre_postCT.pdf'){
  if(!is.null(filename)) pdf(filename,width=4,height=4)
  for (i in scatterplotitle[1:175]){
    title.text = i
    if ((max(before[,i])<0.1) & (max(after[,i])< 0.1)){
      plot(before[,i],after[,i],xlim=c(0,0.1),ylim=c(0,0.1), col=c("firebrick2"), pch= 20 , xlab='before chemotherapy',ylab='afterchemotherapy',main= title.text,cex.main=1.5)
      abline(0, 1)
    } else if ((max(before[,i])<0.3) & (max(after[,i])< 0.3)){
      plot(before[,i],after[,i],xlim=c(0,0.4),ylim=c(0,0.4), col=c("firebrick2"), pch= 20 , xlab='before chemotherapy',ylab='afterchemotherapy',main= title.text,cex.main=1.5)
      abline(0, 1)
    } else {
      plot(before[,i],after[,i],xlim=c(0,0.6),ylim=c(0,0.9), col=c("firebrick2"), pch= 20 , xlab='before chemotherapy',ylab='afterchemotherapy',main= title.text,cex.main=1.5)
      abline(0, 1)
    }
    
  }
  if(!is.null(filename)) dev.off()
}  


