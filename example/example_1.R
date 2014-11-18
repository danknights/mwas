# Example - using the whole pipeline
#
#


#library(biom)
setwd("~/Documents/R/mwas_git")
# Source files
#file.sources = list.files(c("C:/folder1", "C:/folder2"),
file.sources = list.files("lib", pattern="*.R$",
                          #file.sources = list.files(paste(Sys.getenv('MWAS_DIR'),'/lib',sep=''), pattern="*.R$",
                          full.names=TRUE, ignore.case=TRUE)
invisible(sapply(file.sources, source, .GlobalEnv))

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
opts$outdir <- "example/train_model"
opts$method <- "svm"
opts$feat <- FALSE
opts$suppress_relative_abundance_conversion <- FALSE
opts$collapse_table <- FALSE
opts$min_prevalence <- 0.01
opts$transform_type <- "none"
opts$validType <- "cv"
opts$nfolds <- 10

if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)

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
opts$feat_stats_fp <- "test/data/stats/taxon-stats-table.txt"
opts$min_prevalence <- 0.1
opts$method <- "beeswarm"
opts$shorten_taxa <- FALSE
opts$alpha <- 0.001
opts$shorten_taxa <- TRUE
opts$filter_kegg <- FALSE
opts$collapse_table <- FALSE
opts$suppress_relative_abundance_conversion <- FALSE
opts$nplot <- 20

mwas.obj <- import.plot.params(opts)
plot(mwas.obj)

###### test plot gradients
opts <- list()
opts$mode <- "plot"
opts$input_fp <- "test/data/taxa/merged-taxa.txt"
opts$map_fp <- "test/data/gg-map-adults.txt"
opts$category <- "COUNTRY"
opts$outdir <- "example/gradients"
opts$feat_stats <- "test/data/stats/taxon-stats-table.txt"
opts$min_prevalence <- 0.1
opts$method <- "gradients"
opts$shorten_taxa <- FALSE
opts$alpha <- 0.05
opts$shorten_taxa <- TRUE
opts$which_taxa <- 'k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides'
opts$filter_kegg <- FALSE
opts$collapse_table <- FALSE
opts$suppress_relative_abundance_conversion <- FALSE
opts$multiple_axes <- FALSE

mwas.obj <- import.plot.params(opts)
plot(mwas.obj)
