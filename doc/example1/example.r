# this is an example file for a hypothetical data set
source(paste(Sys.getenv('MWAS_DIR'),'/src/lib/util.r',sep=''))
source(paste(Sys.getenv('MWAS_DIR'),'/src/lib/wrap.edgeR.r',sep=''))

library('edgeR')

# load data into data list
d <- load.experiment(map.fp='/path/to/map',
		otu.fp='/path/to/otus-non-rarefied',
		taxa.fp=c('/path/to/taxa-file-1','/path/to/taxa-file-2',...),
		beta.fp=c('/path/to/bdiv-file-1','/path/to/bdiv-file-2',...),
		alpha.fp=c('/path/to/adiv-table')
	)


# d is a data list containing map, taxa, otus, etc.
with(d, {

	# example usage of edgeR with covariates
	ix <- map$Treatment %in% c('Treatment1', 'Treatment2')
	dge.res <- exact.test.edgeR.covariates(x=taxa[[6]][ix,],
										y=map$Treatment[ix],
										covariates=map[ix,c('Gender','Age')])
	tt <- topTags(dge.res)
	print(tt)
})