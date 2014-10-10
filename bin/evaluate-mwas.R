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

source("../lib/I-methods.r") # for OTUs, mapping files operation
#source("../lib/load_library.r")

#source("Documents/FMT_stability/dist_median.r")
#source('~/Documents/FMT_stability/hmp_fmt/lib/health_index_methods.R')

# convert arguments to vector
#myPackages <- c("biom", "optparse", "e1071", "kernlab","randomForest", "glmnet", "pROC", "ROCR")
#load.library(myPackages)
require(biom)
require(optparse) 
require(e1071) 
require(kernlab)
require(randomForest)
require(glmnet)
require(pROC)
#require(ROCR)

####################### Parse INPUT options #####

# make option list and parse command line
option_list <- list(
  make_option(c("-i","--OTU_table_fp"), type="character",
              help="BIOM format or classic format of OTU table [requried]."),
  make_option(c("-m","--map_fp"), type="character",
              help="Mapping file  [required]."),
  make_option(c("-c","--category"), type="character",
              help="Column name in the mapping file [requried]"),
  make_option(c("-t", "--method"),type='character',     #default="RF",
              help="Classifier type [required for model training]"),
  make_option(c("-d", "--mode"),type='character',
              help="Function mode [required]"),
  make_option(c("-k", "--param"),type='character',default="radial",
              help="Classifier parameter, e.g. SVM kernel type [default: %default]"),
  make_option(c("-f", "--fold"),type='numeric',default=10,
              help="Number of folds in cross-validation [default: %default]"),
  make_option(c("-o", "--outdir"),type='character',default=".",
              help="Output directory [default: %default]")
  make_option(c("-v", "--plottype"),type='character', default="heatmap",
              help="Plot type [default: %default]")
  make_option(c("-s", "--feat"),type='logical',default=TRUE,
              help="Option for feature selection [default: %default]")
  make_option(c("-b", "--feat_param"),type='numeric',default=0,
              help="Parameter for feature selection [default: %default]")
  make_option(c("-a", "--statistcs"),type='character',default="linear",
              help="Statistical testing options [default: %default]")
)
opts <- parse_args(OptionParser(option_list=option_list),
                   args=commandArgs(trailing=TRUE))

# create output directory if needed
if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)

######################## Load data #######
mapping <-  load.qiime.mapping.file(opts$map_fp)   # mapping file

if (grep(".biom",opts$map_fp)) {
	biom_table <- read_biom(opts$OTU_table_fp)         # OTU table - biom format
	otus <- t(as.matrix(biom_data(biom_table)))        # OTU table - classic format
}
else {
	trycatch(otus <- read.delim(opts$OTU_table_fp, sep='\t',
	comment='',head=T,row.names=1,check.names=F),error = function(err) 
		print("Couldn't parse OTU table. If BIOM format, use .biom extension"))
}
data.list <- remove.nonoverlapping.samples(map = mapping, otus = otus)
# 
feat.Data <- data.list$otus                        # feature data for training
desired.response <- factor(data.list$map[[opts$categories]]) # desired lables 
(categories)
desired.labels <- desired.response
levels(desired.labels) <- 0:length(levels(desired.labels))-1


if (exists(opts$method, mode="character")) { # training mode
  training.set <- feat.Data
  model.obj <- cross.validation.mwas(feat.Data, desired.labels, nfold=opts$fold, classifier=opts$method, savefile=TRUE, opts)
  ## need feature selection and testing.set
  
}
else{ # tesing mode
  testing.set <- feat.Data 
  ## need feature selection process
  model.obj <- load(opts$model)
}

results <- model.evaluation.mwas(testing.set, opts$model)

### save results ## 
