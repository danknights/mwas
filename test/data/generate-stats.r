m <- read.table('gg-map-adults.txt',sep='\t',head=T,row=1,check=F,comment='')
dim(m)
x <- read.table('taxa/merged-taxa.txt',sep='\t',head=T,row=1,check=F)
dim(x)
x <- t(as.matrix(read.table('taxa/merged-taxa.txt',sep='\t',head=T,row=1,check=F)))
x <- x[rownames(m),]
dim(m)
dim(x)
m$COUNTRY
pvals <- apply(x,2,function(xx) kruskal.test(m$COUNTRY ~ xx)$p.value)
pvals <- apply(x,2,function(xx) kruskal.test(xx ~ m$COUNTRY)$p.value)
qvals <- p.adjust(pvals,'fdr')
sort(qvals)[1:10]
classwise.means <- apply(x,2, function(xx) sapply(split(xx, m$COUNTRY),mean))
dim(classwise.means)
mean(x[m$COUNTRY == levels(m$COUNTRY)[1],1])
classwise.means[,1]
classwise.means <- t(apply(x,2, function(xx) sapply(split(xx, m$COUNTRY),mean)))
res <- data.frame(pvalues=pvals, qvalues=qvals, classwise.means)
colnames(res)
sink('stats/taxon-stats-table.txt'); cat('Taxon\t'); write.table(res,sep='\t',quote=F); sink(NULL)
