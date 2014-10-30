# Example - using the whole pipeline
#
#
source('~/Documents/R/mwas_git/lib/util.load.r')
library(biom)
setwd("~/Documents/R/mwas_git")

opts <- list(make_option(c("-w", "--mode"),type='character',
                    help="Function mode [required]"),
make_option(c("-i","--OTU_fp"), type="character",
            help="BIOM format or classic format of OTU table or other matrix [requried]"),
make_option(c("-m","--map_fp"), type="character",
            help="Mapping file  [required]."),
make_option(c("-c","--category"), type="character",
            help="Column name in the mapping file [requried]"),
make_option(c("-t", "--method"),type='character',     #default="RF",
            help="Classifier type [required for model training]"),
make_option(c("-o", "--outdir"),type='character',default=".",
            help="Output directory [default: %default]")
)



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