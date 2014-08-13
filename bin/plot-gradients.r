#!/usr/bin/env Rscript
# run with Rscript plot-gradients.r -i taxontable -w 'Bacteroides,Prevotella' -o outdir
# 

# REQUIRED GLOBAL VARIABLES: PLEASE EDIT
source(paste(Sys.getenv('MWAS_DIR'),'/lib/gradients.r',sep=''))
source(paste(Sys.getenv('MWAS_DIR'),'/lib/util.r',sep=''))

library('RColorBrewer')
library('optparse')
library('vegan')

# make option list and parse command line
option_list <- list(
    make_option(c("-i","--input_fp"), type="character",
        help="QIIME-formatted input taxon table (tab-delimited text, not biom) [required]."),
    make_option(c("-m","--map_fp"), type="character",default=NULL,
        help="QIIME-formatted mapping file (optional). If provided, only samples in both taxon table and mapping file will be plotted."),
    make_option(c("-c","--column"), type="character",default=NULL,
        help="Name of metadata column to color plot by (optional). If included, does not plot gradients."),
    make_option(c("-d","--distance_fp"), type="character",default=NULL,
        help="QIIME-formatted distance table file (optional). If omitted, the script uses Bray-Curtis distance."),
    make_option(c("-p","--pcoa_fp"), type="character",default=NULL,
        help="QIIME-formatted pcoa table file (optional). If omitted, the script uses Bray-Curtis distance. If included, takes priority over --distance_fp."),
    make_option(c("-w", "--which_taxa"), type="character", default=NULL,
        help="Comma-separated list of taxa to plot [default: plot top --nplot taxa]"),
    make_option(c("-s", "--shorten_taxa"),action='store_true',default=FALSE,
        help="Shorten taxonomy names to lowest defined level. [default: %default]"),
    make_option(c("-x", "--multiple_axes"),action='store_true',default=FALSE,
        help="Show PC1 v PC2, PC1 v PC3, PC2 v PC3 in 3 separate plots. [default: %default]"),
    make_option(c("-n", "--nplot"), type="numeric", default=10,
        help="Number of taxa to plot (in order of decreasing mean). Ignored if --which_taxa exists [default: %default]"),
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

if(opts$shorten_taxa){
	colnames(x) <- shorten.taxonomy(colnames(x))
	taxon.names <- shorten.taxonomy(taxon.names)
}

if(is.null(opts$pcoa_fp)){
	if(is.null(opts$distance_fp)){
		d <- vegdist(x)
	} else {
		d <- read.table(opts$distance_fp,sep='\t',head=T,row=1,check=F)
		d <- d[rownames(x),rownames(x)]
		d <- as.dist(d)
	}
	pc <- cmdscale(d,k=5)
} else {
	pc <- read.table(opts$pcoa_fp,sep='\t',row=1,head=T)
	if(rownames(pc)[nrow(pc)] == '% variation explained'){
		pc <- pc[1:(nrow(pc)-2),1:min(5,ncol(pc))]
	}
	if(mean(rownames(x) %in% rownames(pc)) < 1){
		stop('Taxon table row names do not match PC file row names')
	}
	pc <- pc[rownames(x),]
}

# plots
if(is.null(opts$column)) {
	fp <- sprintf('%s/gradients.pdf',opts$outdir)
} else {
	if(!is.element(opts$column,colnames(m))) stop(paste(opts$column,'not in mapping file.'))
	fp <- sprintf('%s/pcoa.pdf',opts$outdir)
}
if(opts$multiple_axes){
	pdf(fp,width=11,height=3.75)
	par(mfrow=c(1,3))
	combs <- combn(1:3,2)
} else {
	pdf(fp,width=6,height=5)
	combs <- matrix(1:2,ncol=1)
}

if(is.null(opts$column)){
	for(i in seq_along(taxon.names)){
		for(j in 1:ncol(combs)){
			show.gradients(x[,taxon.names[i]], pc[,combs[,j]], incl.legend=TRUE,pt.alpha='CC',
				axis.labels=sprintf('PC%d',combs[,j]),
				title.text=sprintf('%s - PC%d v PC%d',taxon.names[i],combs[1,j],combs[2,j]))
		}
	}
} else {
	for(j in 1:ncol(combs)){
		show.metadata(m[,opts$column], pc[,combs[,j]], incl.legend=TRUE,pt.alpha='CC',
			axis.labels=sprintf('PC%d',combs[,j]),
			title.text=sprintf('%s - PC%d v PC%d',opts$column,combs[1,j],combs[2,j]))
	}
}
dev.off()


