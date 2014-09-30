source("../lib/I-methods.r") # for OTUs, mapping files operation
source("../lib/load_library.r")

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
require(ROCR)

####################### Parse INPUT options #####

# make option list and parse command line
option_list <- list(
  make_option(c("-i","--OTU_table_fp"), type="character",
              help="BIOM format of OTU table [requried]."),
  make_option(c("-m","--map_fp"), type="character",
              help="Mapping file  [required]."),
  make_option(c("-c","--category"), type="character",
              help="Column name in the mapping file [requried]"),
  make_option(c("-m", "--method"),type='character',default="RF",
              help="Classifier type [default: %default]"),
  make_option(c("-k", "--kernel"),type='character',default="radial",
              help="SVM kernel type [default: %default]"),
  make_option(c("-f", "--fold"),type='numeric',default=10,
              help="Number of folds in cross-validation [default: %default]"),
  make_option(c("-o", "--outdir"),type='character',default=".",
              help="Output directory [default: %default]")
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
	trycatch(otus <- read.delim('map-subset-imputed.txt', sep='\t',
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

#Training.Data <- cbind(feat.Data, desired.labels)      # combine feature and labels together
# ix <- which(!is.na(desired.labels))                  # indices of non-NA's

#print(desired.labels)
cv.ind <- sample(dim(feat.data)[1])


training.set <- feat.Data
model.obj <- ml.cross.validation(feat.Data, desired.labels, nfold=opts$fold, classifier=opts$method, savefile=TRUE, opts)
