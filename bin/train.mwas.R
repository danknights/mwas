#!/usr/bin/env Rscript
## Main function calling all testing, model training, model testing and visualization
# Contributors: Hu, Dan
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
# -------
#  Last Update: 10/25/2014
#

######################## Load data #######
mapping <-  load.qiime.mapping.file(opts$map_fp)   # mapping file

if (grep(".biom$",opts$data_table)) {
	biom_table <- read_biom(opts$data_table)         # OTU table - biom format
	otus <- t(as.matrix(biom_data(biom_table)))      # OTU table - classic format
} else {
	trycatch(otus <- read.delim(opts$OTU_table_fp, sep='\t',
	comment='',head=T,row.names=1,check.names=F),error = function(err) 
		print("Couldn't parse OTU table. If BIOM format, use .biom extension"))
}

feat.Data <- otus # feature data for training

response <- droplevels(factor(data.list$map[[opts$category]])) # desired labels 

print(dim(feat.Data))
print(response)

######################## Train and evaluation model #######
model.obj <- cross.validation.mwas(feat.Data, response, nfold=opts$nfolds, classifier=opts$method, savefile=TRUE, opts)
# results <- model.evaluation.mwas(testing.set, opts$model)

### save results ###