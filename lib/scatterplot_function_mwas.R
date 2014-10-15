# scatterplot of 2 conditions included in the mapping file (here Treatment preCT, postCT)
# and perform spearman correlation for each taxon
# require m:  mapping file 
# x: otu table 
# require shorten.taxonomy function
  
source('~/DKpostdoc/src/mwas/lib/util.r')

"run.scatterplot" <- function(m,x, filename='scatterplot_pre_postCT2.pdf'){
  x <- x[order(row.names(x)),]
  m <- m[order(row.names(m)),]
  x2scatterplot <- cbind(m,x)
  x2scatterplot<- x2scatterplot[order(x2scatterplot$Treatment),]
  scatterplotfile <- x2scatterplot
  scatterplotfile2 <-scatterplotfile[,-c(1,2,3,4,5)]
  attributes(scatterplotfile2)
  scatterplotitle2 <- names(scatterplotfile2)
  scatterplotitle2 <- shorten.taxonomy(scatterplotitle2)
  colnames(scatterplotfile2) <- c(scatterplotitle2)
  after <- asin(sqrt(scatterplotfile2[1:15,]))
  before <- asin(sqrt(scatterplotfile2[16:30,]))
  attach(scatterplotfile2)
  attributes(scatterplotfile2)
  scatterplotitle <- names(scatterplotfile2)
  
  if(!is.null(filename)) pdf(filename,width=4,height=4)
  for (i in scatterplotitle[1:16]){
    title.text = i
    if ((max(before[,i])<0.1) & (max(after[,i])< 0.1)){
      plot(before[,i],after[,i],xlim=c(0,0.1),ylim=c(0,0.1), col=c("firebrick2"), pch= 20 , xlab='before chemotherapy',ylab='afterchemotherapy',main= title.text,cex.main=1.5)
      #abine(0,1) #separe en 2 le graph
      abline(before[,i],after[,i])
    } else if ((max(before[,i])<0.3) & (max(after[,i])< 0.3)){
      plot(before[,i],after[,i],xlim=c(0,0.4),ylim=c(0,0.4), col=c("firebrick2"), pch= 20 , xlab='before chemotherapy',ylab='afterchemotherapy',main= title.text,cex.main=1.5)
      abline(before[,i],after[,i])
    } else {
      plot(before[,i],after[,i],xlim=c(0,0.6),ylim=c(0,0.9), col=c("firebrick2"), pch= 20 , xlab='before chemotherapy',ylab='afterchemotherapy',main= title.text,cex.main=1.5)
      abline(before[,i],after[,i])
    }
  }
  if(!is.null(filename)) dev.off()
}  








# perform an unique scatterplot of all the bugs
# and calculate spearman correlation
# require m:  mapping file 
# x: otu table

"run.scatterplot.all.data" <- function(scatterplotfile2, filename='scatterplot_pre_postCT_all_data.pdf'){
  if(!is.null(filename)) pdf(filename,width=4,height=4)
  scatterplotfileAV <- scatterplotfile2[-c(1:15),]
  AV <- stack(scatterplotfileAV)
  scatterplotfileAP <- scatterplotfile2[-c(16:30),]
  AP <- stack(scatterplotfileAP)
  AVAP <- cbind(AV,AP)
  AVAP <- AVAP[,-c(2,4)]
  colnames(AVAP) <- c('before' , 'after')
  plot(AVAP$before, AVAP$after, main="Scatterplot pre post CT ALL", 
       xlab="preCT ", ylab="postCT ", col=c("firebrick2"), pch= 20)
  res <- cor.test(a$before, a$after, method = 'spearman')
  res$estimate 
  res$p.value
  if(!is.null(filename)) dev.off()
  print(paste("rho:",res$estimate, "and p-value:",res$p.value))
}
