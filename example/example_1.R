# Example - using the whole pipeline
#
#
source('~/Documents/R/mwas_git/lib/util.load.r')
source('~/Documents/R/mwas_git/lib/stats.r')
library(biom)
setwd("~/Documents/R/mwas_git")

opts <- list()
opts$mode <- "train"
opts$input_fp <- "test/data/GG_100nt_even10k-adults-s20.biom"
opts$map_fp <- "test/data/gg-map-adults.txt"
opts$category <- "COUNTRY"
opts$outdir <- "example/"

################### train
mapping <- load.qiime.mapping.file('test/data/gg-map-adults.txt')
dim(mapping)
#x <- read.table('test/data/GG_100nt_even10k-adults-biom',sep='\t',head=T,row=1,check=F)
biom_table <- load.qiime.otu.table('test/data/GG_100nt_even10k-adults-s20.biom')
dim(biom_table)
otus <- biom_table[rownames(mapping),]

feat.Data <- otus
response <- as.factor(mapping[,"COUNTRY"])

diff <- differentiation.test(feat.Data, response)
write.differentiation.test.results(diff,  filename='example/differentiated.features.txt')

diff.qvalue <- subset(diff, diff$qvalues < 0.05)

ft<- read.table('example/differentiated.features.txt',sep='\t',head=T,row=1,check=F,quote='"',comment='')
ft.qvalue <-  subset(ft, qvalue < 0.05)
source('~/Documents/R/mwas_git/lib/model.train.r')

best.model <- train.mwas(feat.Data, response, is.feat = FALSE, method="svm")
################### train
opts <- list()
opts$mode <- "train"
opts$input_fp <- "test/data/GG_100nt_even10k-adults-s20.biom"
opts$map_fp <- "test/data/gg-map-adults.txt"
opts$category <- "COUNTRY"
opts$outdir <- "example/"
opts$method <- "svm"
opts$feat <- FALSE

case.mode <- tolower(opts$mode)

mwas.obj <- import.train.params(opts)
model.obj <- train.mwas(mwas.obj)
print("Training a model is finished!")

################# predict
pred_opts <- list()
pred_opts$input_fp <- "test/data/GG_100nt_even10k-adults-s20.biom"
pred_opts$map_fp <- "test/data/gg-map-adults.txt"
pred_opts$category <- "COUNTRY"
pred_opts$outdir <- "example/"
pred_opts$param_fp <- "example/trained_model.rds"

pred.obj <- import.predict.params(pred_opts)
results <- model.evaluation.mwas(pred.obj)
export.mwas(model.eval=results, out.dir=pred_opts$outdir, file.name="prediction")

################# feature selection - RF
opts <- list()
opts$mode <- "train"
opts$input_fp <- "test/data/GG_100nt_even10k-adults-s20.biom"
opts$map_fp <- "test/data/gg-map-adults.txt"
opts$category <- "COUNTRY"
opts$outdir <- "example/"
opts$method <- "svm"
opts$feat <- TRUE

case.mode <- tolower(opts$mode)

mwas.obj <- import.train.params(opts)
model.obj <- train.mwas(mwas.obj)
print("Training a model is finished!")

################# plot test
opts <- list()
opts$mode <- "plot"
opts$input_fp <- "test/data/taxa/merged-taxa.txt"
opts$map_fp <- "test/data/gg-map-adults.txt"
opts$category <- "COUNTRY"
opts$outdir <- "example/"
opts$feat_stats <- "test/data/stats/taxon-stats-table.txt"
opts$min_prevalence <- 0.1
opts$plot_type <- "beeswarm"
opts$shorten_taxa <- FALSE
opts$alpha <- 0.05
opts$shorten_taxa <- TRUE
opts$filter_kegg <- FALSE
opts$collapse_table <- FALSE
opts$suppress_relative_abundance_conversion <- FALSE

mwas.obj <- import.plot.params(opts)
plot(mwas.obj)
###### test plot gradients
opts <- list()
opts$mode <- "plot"
opts$input_fp <- "test/data/taxa/merged-taxa.txt"
opts$map_fp <- "test/data/gg-map-adults.txt"
opts$category <- "COUNTRY"
opts$outdir <- "example/"
opts$feat_stats <- "test/data/stats/taxon-stats-table.txt"
opts$min_prevalence <- 0.1
opts$plot_type <- "scatterplot"
opts$shorten_taxa <- FALSE
opts$alpha <- 0.05
opts$shorten_taxa <- TRUE
opts$which_taxa <- 'k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides'

mwas.obj <- import.plot.params(opts)
plot(mwas.obj)
