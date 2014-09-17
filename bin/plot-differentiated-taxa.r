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
    make_option(c("-t", "--transform_type"), type="character", default="norm-asin-sqrt",
        help="Relative abundance transform type (none, asin-sqrt, or norm-asin-sqrt) [default: norm-asin-sqrt]"),
    make_option(c("-X", "--x_axis_label"), type="character", default="",
        help="Label for x axis [default: blank]"),
    make_option(c("-s", "--shorten_taxa"),action='store_true',default=FALSE,
        help="Shorten taxonomy names to lowest defined level. [default: %default]"),
    make_option(c("-r", "--sort_by_abundance"),action='store_true',default=FALSE,
        help="Sort resulting plots by decreasing relative abundance (instead of significance) [default: %default]"),
    make_option(c("-n", "--nplot"), type="numeric", default=NULL,
        help="Number of taxa to plot (in order of decreasing significance). Ignored if --which_taxa exists [default: %default]"),
	make_option(c("-a","--alpha"), type='numeric', default=.05,
		help='Maximum false discovery rate to report. Ignored if --which_taxa exists or --nplot exists [default: %default]'),
    make_option(c("-O","--category_order"), type="character", default=NULL,
        help="Optional ordering of categories (comma-separated) [default alphabetical]."),	
    make_option(c("-o", "--outdir"), type="character", default='.',
        help="Output directory [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), 
    args=commandArgs(trailing=TRUE))

# create output directory if needed
if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)

# LOAD DATA
x <- t(read.table(opts$input_fp,sep='\t',head=T,row=1,check=F,quote='"'))
if(!is.null(opts$map_fp)){
	m <- read.table(opts$map_fp,sep='\t',head=T,row=1,check=F,comment='',quote='"')
	x <- x[intersect(rownames(x),rownames(m)),,drop=F]
	m <- droplevels(m[rownames(x),,drop=F])
}

# remove rare features
x <- x[,colMeans(x > 0) >= opts$min_prevalence]

# data transform
if(opts$transform_type == 'asin-sqrt'){
	x <- asin(sqrt(x))
} else if(opts$transform_type == 'norm-asin-sqrt'){
	x <- asin(sqrt(x))/asin(sqrt(1))
} else if(opts$transform_type != 'none'){
	stop(paste('Unrecognized data transform type:',opts$transform_type))
}

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
diff.tests <- differentiation.test(x, m[,opts$column], alpha=2, parametric=FALSE)

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

# sort by relative abundance if requested
if(opts$sort_by_abundance){
	hit.ix <- hit.ix[order(colMeans(x)[hit.ix],decreasing=TRUE)]
}

if(opts$shorten_taxa){
	colnames(x) <- shorten.taxonomy(colnames(x))
	taxon.names <- shorten.taxonomy(taxon.names)
}
		
# plots
env <- m[,opts$column]

# fix order if requested
if(!is.null(opts$category_order)){
	level.order <- strsplit(opts$category_order,',')[[1]]
	expected.levels <- sort(unique(as.character(env)))
	if(!identical(sort(unique(level.order)),expected.levels)){
		stop(paste("--category_order list does not contain the same categories",
				   "as the provided mapping file:",
				   paste(sort(unique(env)),collapse=', '),'\n',sep=' '))
	}
	env <- factor(env,levels=level.order)
}

cols <-  c(brewer.pal(9,'Set1'),brewer.pal(9,'Pastel1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[-6]
cols[1:2] <- cols[2:1]
cols <- sprintf('%sbb',cols)
for(i in hit.ix){
	taxon <- x[,i]
	taxon.name <- colnames(x)[i]
	filename.i <- sprintf('%s/beeswarm-%s.pdf',opts$outdir,taxon.name)
	pdf(filename.i,width=8.5,height=11)
	par(oma=rep(c(0,0),2),mar=c(12,8,8,3), cex.axis=.65, cex.lab=.65, cex=2)
	cex.lab <- max(.45, 1 - strwidth(taxon.name, units='figure') / 4)
	cex.lab <- .45
	cex.axis <- .45
	beeswarm(taxon ~ env,corral='random',las=2,
		col='#000000bb',
		bg=cols,
		pch=21,ylab=taxon.name,xlab=opts$x_axis_label,cex.lab=cex.lab, cex.axis=cex.axis)
	dev.off()
}
