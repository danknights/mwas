# Example - using the whole pipeline
#
#
source('~/Documents/R/mwas_git/lib/util.load.r')
library(biom)
setwd("~/Documents/R/mwas_git")

mapping <- load.qiime.mapping.file('test/data/gg-map-adults.txt')
dim(mapping)
#x <- read.table('test/data/GG_100nt_even10k-adults-biom',sep='\t',head=T,row=1,check=F)
biom_table <- load.qiime.otu.table('test/data/GG_100nt_even10k-adults-s20.biom')
otu <- t(as.matrix(biom_table$shape))
dim(otu)
otu <- otu[rownames(mapping),]

data.list <- remove.nonoverlapping.samples(map = mapping, otus = otu)
feat.Data <- data.list$otus # feature data for training


m$COUNTRY
pvals <- apply(otu,2,function(xx) kruskal.test(m$COUNTRY ~ xx)$p.value)
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
