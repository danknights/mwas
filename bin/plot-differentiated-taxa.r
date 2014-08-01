#!/usr/bin/env Rscript
# Usage: run non-parametric test for differentiation according to metadata:
#
# run with Rscript make-beeswarm-plots.r -i taxontable -m map.txt -c column_name -o outdir
# 
# Usage: plot specific taxa (must match row name exactly from taxontable)
# run with Rscript make-beeswarm-plots.r -i taxontable -w 'Bacteroides,Prevotella' -o outdir

# REQUIRED GLOBAL VARIABLES: PLEASE EDIT
source(paste(Sys.getenv('MWAS_DIR'),'/lib/gradients.r',sep=''))
source(paste(Sys.getenv('MWAS_DIR'),'/lib/util.r',sep=''))
source(paste(Sys.getenv('MWAS_DIR'),'/lib/stats.r',sep=''))

library('RColorBrewer')
library('optparse')
library('vegan')
library('beeswarm')

# make option list and parse command line
option_list <- list(
    make_option(c("-i","--input_fp"), type="character",
        help="QIIME-formatted input taxon table (tab-delimited text, not biom) [required]."),
    make_option(c("-m","--map_fp"), type="character",default=NULL,
        help="QIIME-formatted mapping file (optional). If provided, only samples in both taxon table and mapping file will be plotted."),
    make_option(c("-c","--column"), type="character",default=NULL,
        help="Name of metadata column to color plot by (optional). If included, does not plot gradients."),
    make_option(c("-M","--min_prevalence"), type="numeric",default=.1,
        help="Minimum fraction of samples in which taxon must be present to be included [default: %default]."),
    make_option(c("-w", "--which_taxa"), type="character", default=NULL,
        help="Comma-separated list of taxa to plot [default: plot top --nplot taxa]"),
    make_option(c("-s", "--shorten_taxa"),action='store_true',default=FALSE,
        help="Shorten taxonomy names to lowest defined level. [default: %default]"),
    make_option(c("-n", "--nplot"), type="numeric", default=NULL,
        help="Number of taxa to plot (in order of decreasing significance). Ignored if --which_taxa exists [default: %default]"),
	make_option(c("-a","--alpha"), type='numeric', default=.05,
		help='Maximum false discovery rate to report. Ignored if --which_taxa exists or --nplot exists [default: %default]'),
    make_option(c("-o", "--outdir"), type="character", default='.',
        help="Output directory [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), 
    args=commandArgs(trailing=TRUE))

# create output directory if needed
if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)

# LOAD DATA
x <- t(read.table(opts$input_fp,sep='\t',head=T,row=1,check=F))
if(!is.null(opts$map_fp)){
	m <- read.table(opts$map_fp,sep='\t',head=T,row=1,check=F,comment='')
	x <- x[intersect(rownames(x),rownames(m)),,drop=F]
	m <- droplevels(m[rownames(x),,drop=F])
}

# remove rare features
x <- x[,colMeans(x > 0) >= opts$min_prevalence]

# asin sqrt
x <- asin(sqrt(x))

# check that taxon.names are in taxon table
if(is.null(opts$which_taxa)){
	taxon.names <- colnames(x)[rev(order(colMeans(x)))]
	taxon.names <- taxon.names[1:min(opts$nplot, length(taxon.names))]
} else {
	taxon.names <- strsplit(opts$which_taxa,',')[[1]]

	if(!all(taxon.names %in% colnames(x))){
		stop(paste('The following taxa are not present in the taxon table:',
					paste(taxon.names[!(taxon.names %in% colnames(x))],collapse=', '),
					'\n'))
	}
} 
if(!is.element(opts$column,colnames(m))) stop(paste(opts$column,'not in mapping file.'))

# identify differentiated features
diff.tests <- differentiation.test(x, m[,opts$column], parametric=FALSE)

if(is.null(opts$which_taxa)){
	if(is.null(opts$nplot)){
		hit.ix <- which(diff.tests$qvalues <= opts$alpha)
	} else {
		hit.ix <- order(diff.tests$pvalues)[1:min(opts$nplot,length(diff.tests$pvalues))]
	}
} else {
	hit.ix <- match(opts$which_taxa, colnames(x))
}

if(length(hit.ix) == 0){
	stop('No hits found.\n')
} else {
	cat("There were",length(hit.ix),"taxa significant at FDR of",opts$alpha,'\n')
}

# save hits
write.differentiation.test.results(diff.tests, filename=sprintf('%s/test_results.txt',opts$outdir))

if(opts$shorten_taxa){
	colnames(x) <- shorten.taxonomy(colnames(x))
	taxon.names <- shorten.taxonomy(taxon.names)
}
		
# plots
env <- m[,opts$column]
fp <- sprintf('%s/beeswarms.pdf',opts$outdir)
pdf(fp,width=4,height=4)
par(oma=c(4,0,0,0),mar=c(4,4,.5,.5))
cols <- y.colors <- c(brewer.pal(9,'Set1'),brewer.pal(9,'Pastel1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[-6]
cols <- sprintf('%sbb',cols)
for(i in hit.ix){
	taxon <- x[,i]
	taxon.name <- colnames(x)[i]
	beeswarm(taxon ~ env,corral='random',las=2,
		col='#000000bb',
		bg=cols,
		pch=21,ylab=taxon.name,xlab='')
}
dev.off()
