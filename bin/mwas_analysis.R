#!/usr/bin/env Rscript
## Command line function calling all testing, model training, model testing and visualization
# 
# Contributors: Hu
# ------  
#  input: 
#   <general options>
#   1. w - mode - "train", "predict", "plot", "statistics"
#   2. i - OTU table - BIOM or classic format (required)
#   3. m - map file   (required)
#   4. c - category name(s) (required)
#   5. o - output directory: save classification model, model evaluation results, visualization file etc.
#   <classifier related>
#   6. t - method - classifier method (optional); if omitted then need to specify a trained model file 
#   7. k - classifier parameters- kernel function only for SVM classification method
#   8. e - validation type (k-fold cross-validation - cv; jackknifing - jk)
#   9. f - number of folds in cross-validation [default: 10]
#  10. s - feature selection option - TRUE or FALSE
#  11. b - feature selection parameter - threshold in RF
#  <plot> (options are capital)
#  12. T - plot type - heatmap, beeswarm, violin, gradient
#  13.
#  <statistical testing>
#  13. a - statistical options
#
# ----
# Last update: 10/25/2014
#

source(paste(Sys.getenv('MWAS_DIR'),'/lib/I-methods.r',sep='')) # for OTUs, mapping files operation
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
  make_option(c("-w", "--mode"),type='character',
              help="Function mode [required]"),
  make_option(c("-i","--OTU_fp"), type="character",
              help="BIOM format or classic format of OTU table or other matrix [requried]"),
  make_option(c("-m","--map_fp"), type="character",
              help="Mapping file  [required]."),
  make_option(c("-c","--category"), type="character",
              help="Column name in the mapping file [requried]"),
  make_option(c("-t", "--method"),type='character',     #default="RF",
              help="Classifier type [required for model training]"),
  make_option(c("-k", "--param_fp"),type='character',default="radial",
              help="For training: classifier parameter, e.g. SVM kernel type [default: %default]; 
              For predicting: the directory of the trained model.
              For plotting: distance matrix file."),
  make_option(c("-e", "--validType"),type='character',default="cv",
              help="Validation type (k-fold cross-validation [cv] or Jackknifing [jk]) [default: %default]"),
  make_option(c("-f", "--fold"),type='numeric',default=10,
              help="Number of folds in cross-validation [default: %default]"),
  make_option(c("-o", "--outdir"),type='character',default=".",
              help="Output directory [default: %default]")
  make_option(c("-v", "--plot_type"),type='character', default="heatmap",
              help="Plot type [default: %default]")
  make_option(c("-s", "--feat"),action="store_true", default=TRUE
              help="Flag for feature selection [default: %default]")
  make_option(c("-b", "--featparam"),type='numeric',default=0,
              help="Parameter for feature selection [default: %default]")
  make_option(c("-a", "--statistcs"),type='character',default="linear",
              help="Statistical testing options [default: %default]")
)
opts <- parse_args(OptionParser(option_list=option_list),
                   args=commandArgs(trailing=TRUE))

# create output directory if needed
if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)

######################## Mode selection   #######
# Four available modes:
#  1. train
#  2. predict
#  3. plot
#  4. statistics
case.mode <- tolower(opts$mode) # case insensitive

switch(case.mode, 
       train = {
         ########################    Load data     #######
         mwas.obj <- import.train.params(opts, type="train")
         train.mwas(mwas.obj)
         print("Training a model is finished!")
       }, 
       predict = {
         mwas.obj <- import.predict.params(opts)
         results <- model.evaluation.mwas(mwas.obj)
         export.mwas(results)
       },
       plot = {
         
       },
       statistics = {
         
       },
       stop("Please specify a function mode: train, predict, plot, statistics.")
       )

######################## 
########################    Load data     #######
table.data <- import.mwas(opts)
