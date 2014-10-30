# Training wrapper function
# Contributors: Hu, Dan
# ------  
#  input: 
#     opts : options from the user

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