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

# Source files
#file.sources = list.files(c("C:/folder1", "C:/folder2"),
#file.sources = list.files("lib", pattern="*.R$",
file.sources = list.files(paste(Sys.getenv('MWAS_DIR'),'/lib',sep=''), pattern="*.R$",
                          full.names=TRUE, ignore.case=TRUE)
invisible(sapply(file.sources, source, .GlobalEnv))

####################### Parse INPUT options #####
require(optparse, quietly=TRUE, warn.conflicts=FALSE)
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
  make_option(c("-o", "--outdir"),type='character',default=".",
              help="Output directory [default: %default]"),
  # classification parameters
  make_option(c("-t", "--method"),type='character',     #default="RF",
              help="Classifier type [required for model training]"),
  make_option(c("-k", "--param_fp"),type='character',default="radial",
              help="For training: classifier parameter, e.g. SVM kernel type [default: %default]; 
              For predicting: the directory of the trained model."),
  make_option(c("-e", "--validType"),type='character',default="cv",
              help="Validation type (k-fold cross-validation [cv] or Jackknifing [jk]) [default: %default]"),
  make_option(c("-f", "--fold"),type='numeric',default=10,
              help="Number of folds in cross-validation [default: %default]"),
  make_option(c("-s", "--feat"),action="store_true", default=FALSE,
              help="Flag for feature selection [default: %default]"),
  make_option(c("-b", "--feat_param"),type='numeric',default=0,
              help="Parameter for feature selection [default: %default]"),
  # Statistical test parameters
  make_option(c("-a", "--statistcs"),type='character',default="linear",
              help="Statistical testing options [default: %default]"),
  # plot parameters
  make_option(c("-V", "--plot_type"),type='character', default="heatmap",
              help="Plot type [default: %default]"),
  make_option(c("-D","--distance_fp"), type="character",default=NULL,
              help="QIIME-formatted distance table file (optional). If omitted, the script uses Bray-Curtis distance."),
  make_option(c("-P","--pcoa_fp"), type="character",default=NULL,
              help="QIIME-formatted pcoa table file (optional). If omitted, the script uses Bray-Curtis distance. If included, takes priority over --distance_fp."),
  make_option(c("-W", "--which_taxa"), type="character", default=NULL,
              help="Comma-separated list of taxa to plot [default: plot top --nplot taxa]"),
  make_option(c("-S", "--shorten_taxa"),action='store_true',default=FALSE,
              help="Shorten taxonomy names to lowest defined level. [default: %default]"),
  make_option(c("-M", "--multiple_axes"),action='store_true',default=FALSE,
              help="Show PC1 v PC2, PC1 v PC3, PC2 v PC3 in 3 separate plots. [default: %default]"),
  make_option(c("-N", "--nplot"), type="numeric", default=10,
              help="Number of taxa to plot (in order of decreasing mean). Ignored if --which_taxa exists [default: %default]"),
  make_option(c("-T", "--transform_type"), type="character", default="norm_asin_sqrt",
              help="Relative abundance transform type (none, asin_sqrt, or norm_asin_sqrt) [default: norm_asin_sqrt]"),
  make_option(c("-X", "--x_axis_label"), type="character", default="",
              help="Label for x axis [default: blank]"),
  make_option(c("-C", "--suppress_relative_abundance_conversion"),action='store_true',default=FALSE,
              help="Do not convert input to relative abundances (assumes already relative abundances). [default: %default]"),
  make_option(c("-R", "--sort_by_abundance"),action='store_true',default=FALSE,
              help="Sort resulting plots by decreasing relative abundance (instead of significance) [default: %default]"),
  make_option(c("-A","--alpha"), type='numeric', default=.05,
              help='Maximum false discovery rate to report. Ignored if --which_taxa exists or --nplot exists [default: %default]'),
  make_option(c("-O","--category_order"), type="character", default=NULL,
              help="Optional ordering of categories (comma-separated) [default alphabetical].")
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
case_mode <- tolower(opts$mode) # case insensitive

switch(case_mode, 
       train = { mwas.obj <- import.train.params(opts)
                 train.mwas(mwas.obj)
                 print("Training is finished!")
       }, 
       predict = { mwas.obj <- import.predict.params(opts)
                   results <- model.evaluation.mwas(mwas.obj)
                   export.mwas(model.eval=results, out.dir=opts$outdir, file.name="prediction_results")
         #print("Prediction is finished!")
       }, 
       plot = { mwas.obj <- import.plot.params(opts)
                plot(mwas.obj)
         #print("Visualization")
       }, 
       statistics = { mwas.obj <- import.stats.params(opts)
                      model.statistical.test.mwas(mwas.obj)
         #print("statistics")
       }, 
       stop("Please specify a function mode: train, predict, plot, statistics.")
)
