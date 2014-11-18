#!/usr/bin/env Rscript
## Command line function calling all testing, model training, model testing and visualization
# 
# Contributor: Hu
# ------  
#  input: 
#   <general options>
#   1. w - mode - "train", "predict", "plot", "statistics"
#   2. i - OTU table - BIOM or classic format (required)
#   3. m - map file   (required)
#   4. c - category name(s) (required)
#   5. o - output directory: save classification model, model evaluation results, visualization file etc. [default: "."]
#   6. M - method: classifier type for "training", trained classifier model file name for "predict", 
#                  plot type for "plot", statistical test type for "statistics"
#   <preprocessing options>
#   7. t - transformation type [default: none]
#   8. r - suppress conversion of relative abundance from absolute abunance [default: FALSE]
#   9. p - minimum prevalence of OTU to keep [default: 0.01 %]
#  10. b - whether collapse the table by correlation [default: FALSE]
#   <classifier related>
#  11. C - classifier parameters- kernel function only for SVM classification method
#  12. v - validation type (k-fold cross-validation - cv; jackknifing - jk)
#  13. n - number of folds in cross-validation [default: 10]
#  14. f - Is feature selection option [default: FALSE]
#  15. s - feature selection parameter [default: 0]
#  <plot> (options are capital)
#  16. F - feature statistics table file from differentiated test
#  17. P - PCoA matrix file 
#  18. D - distance matrix file
#  19. T - which taxa to plot
#  20. S - whether shorten taxon names [default: FALSE]
#  21. A - False discovery rate, alpha [default: NULL]
#  22. X - whhether plot pair-wise PC plots [default: FALSE]
#  23. K - For heatmap plot, whether filter kegg list [default: FALSE]
#  <statistical testing>
#  24. a - statistical options
#
# ----
# Last update: 11/12/2014
#

# Source files
#file.sources = list.files(c("C:/folder1", "C:/folder2"),
file.sources = list.files("lib", pattern="*.R$",
#file.sources = list.files(paste(Sys.getenv('MWAS_DIR'),'/lib',sep=''), pattern="*.R$",
                          full.names=TRUE, ignore.case=TRUE)
invisible(sapply(file.sources, source, .GlobalEnv))

####################### Parse INPUT options #####
require(optparse, quietly=TRUE, warn.conflicts=FALSE)
# make option list and parse command line
option_list <- list(
  make_option(c("-w", "--mode"),type='character',
              help="Function mode [required]"),
  # preprocessing parameters
  make_option(c("-t", "--transform_type"), type="character", default="none",
              help="Relative abundance transform type (none, asin_sqrt, or norm_asin_sqrt) [default: norm_asin_sqrt]"),
  make_option(c("-r", "--suppress_relative_abundance_conversion"),action='store_true',default=FALSE,
              help="Do not convert input to relative abundances (assumes already relative abundances). [default: %default]"),
  make_option(c("-p","--min_prevalence"), type="numeric",default=0.0,
              help="Minimum fraction of samples in which taxon must be present to be included [default: %default]."),
  make_option(c("-b","--collapse_table"), type="logical", action='store_true', default=FALSE,
              help="Collapse otu table by correlation."), 
  # universal paramters
  make_option(c("-i","--input_fp"), type="character",
              help="BIOM format or classic format of OTU table or other matrix [requried]"),
  make_option(c("-m","--map_fp"), type="character",
              help="Mapping file  [required]."),
  make_option(c("-c","--category"), type="character",
              help="Column name in the mapping file [requried]"),
  make_option(c("-o", "--outdir"),type='character',default=".",
              help="Output directory [default: %default]"),
  make_option(c("-M", "--method"),type='character',  
              help="Classifier type [required for model training] or trained model file [required for predicting] 
              or plot type [required for plotting]."),
  # classification parameters
  make_option(c("-C", "--method_param"),type='character',default="radial",
              help="For training: classifier parameter, e.g. SVM kernel type [default: %default]."),
  make_option(c("-v", "--validType"),type='character',default="cv",
              help="Validation type (k-fold cross-validation [cv] or Jackknifing [jk]) [default: %default]"),
  make_option(c("-n", "--nfolds"),type='numeric',default=10,
              help="Number of folds in cross-validation [default: %default]"),
  make_option(c("-f", "--is_feat"),type="logical", action="store_true", default=FALSE,
              help="Flag for feature selection [default: %default]"),
  make_option(c("-s", "--feat_param"),type='numeric',default=0.0,
              help="Parameter for feature selection [default: %default]. 
              For plotting options, it specifies the directory of feature statistics file."),
  # plot parameters
  make_option(c("-F","--feat_stats_fp"), type="character",default=NULL,
              help="Feature statistics table file from differentiated test."),
  make_option(c("-D","--distance_fp"), type="character",default=NULL,
              help="QIIME-formatted distance table file (optional). If omitted, the script uses Bray-Curtis distance."),
  make_option(c("-P","--pcoa_fp"), type="character",default=NULL,
              help="QIIME-formatted pcoa table file (optional). If omitted, the script uses Bray-Curtis distance. If included, takes priority over --distance_fp."),
  make_option(c("-T", "--which_taxa"), type="character", default=NULL,
              help="Comma-separated list of taxa to plot [default: plot top --nplot taxa]"),
  make_option(c("-S", "--shorten_taxa"),type="logical", action='store_true',default=FALSE,
              help="Shorten taxonomy names to lowest defined level. [default: %default]"),
  make_option(c("-X", "--multiple_axes"),type="logical", action='store_true',default=FALSE,
              help="Show PC1 v PC2, PC1 v PC3, PC2 v PC3 in 3 separate plots. [default: %default]"),
  make_option(c("-A","--alpha"), type='numeric', default=NULL,
              help='Maximum false discovery rate to report. Ignored if --which_taxa exists exists [default: %default]'),
  make_option(c("-K","--filter_kegg"), type="logical", action='store_true', default=FALSE,
              help='Whether filter kegg list'),
  make_option(c("-N", "--nplot"), type="numeric", default=50,
              help="Number of taxa to plot (in order of decreasing significance). Ignored if --which_taxa exists [default: %default]"),
  # Statistical test parameters
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

case_mode <- tolower(opts$mode) # case insensitive

switch(case_mode, 
       train = { mwas.obj <- import.train.params(opts)
                 train.mwas(mwas.obj)
                 #print("Training is finished!")
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
