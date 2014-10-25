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
dim(biom_table)
otus <- biom_table[rownames(mapping),]

feat.Data <- otus
response <- as.factor(mapping[,"COUNTRY"])

source('~/Documents/R/mwas_git/lib/train.r')

best.model <- train.mwas(feat.Data, response, is.feat = FALSE)