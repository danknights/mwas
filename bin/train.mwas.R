#!/usr/bin/env Rscript
## main function calling all testing, model training, model testing and visualization
# ------  
#  input: 
#  1. OTU table - BIOM or classic format (required)
#  2. map file   (required)
#  3. category name(s) (required)
#  4. mode - "train", "test", "plot", "statistics"
#  5. method - classifier method (optional); if omitted then need to specify a trained model file 
#  6. classifier parameters- kernel function only for SVM classification method
#  7. number of folds in cross-validation [default: 10]
#  8. output directory: save classification model, model evaluation results, visualization file etc.
#  9. plot type - heatmap, beeswarm, violin, gradient
#  10. feature selection option - TRUE or FALSE
#  11. feature selection parameter - threshold in RF
#  12. statistical options
#

require(biom, quietly=TRUE, warn.conflicts=FALSE)
require(optparse, quietly=TRUE, warn.conflicts=FALSE) 
require(e1071, quietly=TRUE, warn.conflicts=FALSE) 
require(kernlab, quietly=TRUE, warn.conflicts=FALSE)
require(randomForest, quietly=TRUE, warn.conflicts=FALSE)
require(glmnet, quietly=TRUE, warn.conflicts=FALSE)
require(pROC, quietly=TRUE, warn.conflicts=FALSE)

source(paste(Sys.getenv('MWAS_DIR'),'/lib/gradients.r',sep=''))
source(paste(Sys.getenv('MWAS_DIR'),'/lib/util.r',sep=''))
source(paste(Sys.getenv('MWAS_DIR'),'/lib/util.load.r',sep=''))
source(paste(Sys.getenv('MWAS_DIR'),'/lib/stats.r',sep=''))

####################### Parse INPUT options #####

# make option list and parse command line
option_list <- list(
  make_option(c("-i","--data_table"), type="character",
              help="OTU table (BIOM or classic format), taxon table, other feature table [requried]."),
  make_option(c("-m","--map_fp"), type="character",
              help="Mapping file  [required]."),
  make_option(c("-c","--category"), type="character",
              help="Column name in the mapping file for classification/regression. If column is numeric, assumes regression [required]"),
  make_option(c("-t", "--method"),type='character', default="RF",
              help="Classifier type: RF, SVM, KNN [default %default]"),
  make_option(c("-p", "--param"),type='character',default="SVM-kernel:radial,RF-ntree:1000",
              help="Comma-separated list of classifier parameters, e.g. SVM kernel type [default: %default]"),
  make_option(c("-f", "--nfolds"),type='numeric',default=10,
              help="Number of folds in cross-validation used for modeling evaluation (-1 means leave-one-out) [default: %default]"),
  make_option(c("-s", "--feat"),type='logical',default=TRUE,
              help="Option for feature selection [default: %default]"),
  make_option(c("-b", "--feat_param"),type='numeric',default=0,
              help="Parameter for feature selection [default: %default]"),
  make_option(c("-o", "--outdir"),type='character',default=".",
              help="Output directory [default: %default]")
)
opts <- parse_args(OptionParser(option_list=option_list),
                   args=commandArgs(trailing=TRUE))

# create output directory if needed
if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)

######################## Load data #######
mapping <-  load.qiime.mapping.file(opts$map_fp)   # mapping file


if (grep(".biom$",opts$data_table)) {
	biom_table <- read_biom(opts$data_table)         # OTU table - biom format
	otus <- t(as.matrix(biom_data(biom_table)))        # OTU table - classic format
} else {
	trycatch(otus <- read.delim(opts$OTU_table_fp, sep='\t',
	comment='',head=T,row.names=1,check.names=F),error = function(err) 
		print("Couldn't parse OTU table. If BIOM format, use .biom extension"))
}

data.list <- remove.nonoverlapping.samples(map = mapping, otus = otus)
feat.Data <- data.list$otus # feature data for training

response <- droplevels(factor(data.list$map[[opts$category]])) # desired labels 

print(dim(feat.Data))
print(response)

# model.obj <- cross.validation.mwas(feat.Data, desired.labels, nfold=opts$fold, classifier=opts$method, savefile=TRUE, opts)
# results <- model.evaluation.mwas(testing.set, opts$model)

### save results ###