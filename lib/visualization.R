"plot.beeswarm.dt" <- function(cols, x_axis_label, hit.ix, env, outdir) {
	require(beeswarm)
	cols <-  c(brewer.pal(9,'Set1'),brewer.pal(9,'Pastel1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[-6]
	cols[1:2] <- cols[2:1]
	cols <- sprintf('%sbb',cols)
	if (!length(hit.ix)) hit.ix = seq(1,dim(x)[2])
	print(x[,1])
	print(env)
	for(i in hit.ix){
		taxon <- x[,i]
		taxon.name <- colnames(x)[i]
		filename.i <- sprintf('%s/beeswarm-%s.pdf',outdir,taxon.name)
		pdf(filename.i,width=8.5,height=11)
		par(oma=rep(c(0,0),2),mar=c(12,8,8,3), cex.axis=.65, cex.lab=.65, cex=2)
		cex.lab <- max(.45, 1 - strwidth(taxon.name, units='figure') / 4)
		cex.lab <- .45
		cex.axis <- .45
		beeswarm(taxon ~ env,corral='random',las=2,
			col='#000000bb',
			bg=cols,
			pch=21,ylab=taxon.name,xlab=x_axis_label,cex.lab=cex.lab, cex.axis=cex.axis)
		bxplot(taxon ~ env,add=TRUE)
		dev.off()
	}
}

# Scatterplot of 2 conditions included in the mapping file (here Treatment preCT, postCT)
# and perform spearman correlation for each taxon
#
# Contributors: Emmanuel
# -------
# Input:
#    m:  mapping file
# x: otu table
# -------
# Output:
#   save plot as a PDF file
# -------
# Last update: 10/25/2014
#

source('~/DKpostdoc/src/mwas/lib/util.r')

"run.scatterplot" <- function(m,x, filename=NULL, ...){
    x <- x[order(row.names(x)),]
    m <- m[order(row.names(m)),]
    x2scatterplot <- cbind(m,x)
    x2scatterplot<- x2scatterplot[order(x2scatterplot$Treatment),]
    scatterplotfile <- x2scatterplot
    scatterplotfile2 <-scatterplotfile[,-c(1,2,3,4,5)]
    attributes(scatterplotfile2)
    scatterplotitle2 <- names(scatterplotfile2)
    scatterplotitle2 <- shorten.taxonomy(scatterplotitle2)
    colnames(scatterplotfile2) <- c(scatterplotitle2)
    after <- asin(sqrt(scatterplotfile2[1:15,]))
    before <- asin(sqrt(scatterplotfile2[16:30,]))
    attach(scatterplotfile2)
    attributes(scatterplotfile2)
    scatterplotitle <- names(scatterplotfile2)
    
    if(!is.null(filename)) pdf(filename,width=4,height=4)
    for (i in scatterplotitle[1:16]){
        title.text = i
        if ((max(before[,i])<0.1) & (max(after[,i])< 0.1)){
            plot(before[,i],after[,i],xlim=c(0,0.1),ylim=c(0,0.1), col=c("firebrick2"), pch= 20 , xlab='before chemotherapy',ylab='afterchemotherapy',main= title.text,cex.main=1.5)
            #abine(0,1) #separe en 2 le graph
            abline(before[,i],after[,i])
        } else if ((max(before[,i])<0.3) & (max(after[,i])< 0.3)){
            plot(before[,i],after[,i],xlim=c(0,0.4),ylim=c(0,0.4), col=c("firebrick2"), pch= 20 , xlab='before chemotherapy',ylab='afterchemotherapy',main= title.text,cex.main=1.5)
            abline(before[,i],after[,i])
        } else {
            plot(before[,i],after[,i],xlim=c(0,0.6),ylim=c(0,0.9), col=c("firebrick2"), pch= 20 , xlab='before chemotherapy',ylab='afterchemotherapy',main= title.text,cex.main=1.5)
            abline(before[,i],after[,i])
        }
    }
    if(!is.null(filename)) dev.off()
}

# perform an unique scatterplot of all the bugs
# and calculate Spearman correlation
# require m:  mapping file
# x: otu table

"run.scatterplot.all.data" <- function(scatterplotfile2, filename='scatterplot_pre_postCT_all_data.pdf'){
    if(!is.null(filename)) pdf(filename,width=4,height=4)
    scatterplotfileAV <- scatterplotfile2[-c(1:15),]
    AV <- stack(scatterplotfileAV)
    scatterplotfileAP <- scatterplotfile2[-c(16:30),]
    AP <- stack(scatterplotfileAP)
    AVAP <- cbind(AV,AP)
    AVAP <- AVAP[,-c(2,4)]
    colnames(AVAP) <- c('before' , 'after')
    plot(AVAP$before, AVAP$after, main="Scatterplot pre post CT ALL",
    xlab="preCT ", ylab="postCT ", col=c("firebrick2"), pch= 20)
    res <- cor.test(a$before, a$after, method = 'spearman')
    res$estimate
    res$p.value
    if(!is.null(filename)) dev.off()
    print(paste("rho:",res$estimate, "and p-value:",res$p.value))
}

# Beeswarm of several conditions included in the mapping file
#
# Contributors: Emmanuel
# ------
# input:
#        x:  otu table
#      map:  mapping file
# ------
# output:
#     save the plot as PDF file
# ------
# Last update: 10/25/2014
#

#source('~/DKpostdoc/src/mwas/lib/util.load.r')

"run.beeswarm" <- function(x, map, filename='beeswarm_plot.pdf'){
    x <- x[order(row.names(x)),]
    m <- map[order(row.names(m)),]
    
    x2beeswarmfile <- cbind(m,x)
    x2beeswarmfile <- x2beeswarmfile[order(x2beeswarmfile$Treatment,decreasing = FALSE),]
    beeswarmfile <- x2beeswarmfile
    attributes(beeswarmfile)
    beeswarmfiletitle <- names(beeswarmfile)
    beeswarmfiletitle <- shorten.taxonomy(beeswarmfiletitle)
    colnames(beeswarmfile) <- c(beeswarmfiletitle)
    category <- beeswarmfile$Treatment
    beeswarmfile2 <-beeswarmfile[,-c(1,2,4,5)]
    beeswarmfile3 <- asin(sqrt(beeswarmfile2[-1]))
    
    require(beeswarm)
    
    if(!is.null(filename)) pdf(filename,width=4,height=4)
    for (i in 1:length(beeswarmfile3)){
        beeswarm(beeswarmfile3[,i] ~ category, data = beeswarmfile2,
        pch = 16,
        col = rainbow(8),
        labels = c("postCT", "preCT"),
        main = colnames(beeswarmfile3)[i],
        xlab = "", ylab ="")
        
    }
    if(!is.null(filename)) dev.off()
}

#!/usr/bin/env Rscript
# run with Rscript plot-gradients.r -i taxontable -w 'Bacteroides,Prevotella' -o outdir

# REQUIRED GLOBAL VARIABLES: PLEASE EDIT
source(paste(Sys.getenv('MWAS_DIR'),'/lib/gradients.r',sep=''))
source(paste(Sys.getenv('MWAS_DIR'),'/lib/util.r',sep=''))

require('RColorBrewer')
require('optparse')
require('vegan')

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
    # check rownames in mapping file matrix
    missing.taxa.samples <- setdiff(rownames(x), rownames(m))
    missing.map.samples <- setdiff(rownames(m), rownames(x))
    if(length(missing.taxa.samples) > 0){
        stop(sprintf('\n\nError: one or more sample names from taxonomy table (%s, ...) not present in metadata table (%s, ...).',
        paste(sort(missing.taxa.samples)[1:5],collapse=', '),
        paste(sort(missing.map.samples)[1:5],collapse=', ')))
    }
    
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
        # check rownames in distance matrix
        missing.taxa.samples <- union(setdiff(rownames(x), rownames(d)), setdiff(rownames(x), colnames(d)))
        missing.distance.samples <- union(setdiff(rownames(d), rownames(x)), setdiff(colnames(d), rownames(x)))
        if(length(missing.taxa.samples) > 0){
            stop(sprintf('\n\nError: one or more sample names from taxonomy table (%s, ...) not present in distance table (%s, ...).',
            paste(sort(missing.taxa.samples)[1:5],collapse=', '),
            paste(sort(missing.distance.samples)[1:5],collapse=', ')))
        }
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
source(paste(Sys.getenv('MWAS_DIR'),'/lib/visualization.R',sep=''))

require('RColorBrewer')
require('optparse')
require('vegan')

# make option list and parse command line
option_list <- list(
make_option(c("-i","--input_fp"), type="character",
help="QIIME-formatted input taxon, OTU, or feature table (see --input_type) (tab-delimited text, not biom) [required]."),
make_option(c("-I","--input_format"), type="character", default='Taxon table',
help="Format of classic QIIME-formatted input table: \"Taxon table\" or \"OTU table\"  [default: %default]."),
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
make_option(c("-C", "--suppress_relative_abundance_conversion"),action='store_true',default=FALSE,
help="Do not convert input to relative abundances (assumes already relative abundances). [default: %default]"),
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
if(opts$input_format == "Taxon table"){
    x <- t(read.table(opts$input_fp,sep='\t',head=T,row=1,check=F,quote='"'))
} else if(opts$input_format == "OTU table") {
    x <- t(read.table(opts$input_fp,sep='\t',head=T,row=1,check=F,quote='"',comment='',skip=1))
} else {
    stop("Unknown input format\n")
}
if(!is.null(opts$map_fp)){
    m <- read.table(opts$map_fp,sep='\t',head=T,row=1,check=F,comment='',quote='"')
    # check rownames in mapping file matrix
    missing.taxa.samples <- setdiff(rownames(x), rownames(m))
    missing.map.samples <- setdiff(rownames(m), rownames(x))
    if(length(missing.taxa.samples) > 0){
        stop(sprintf('\n\nError: one or more sample names from taxonomy table (%s, ...) not present in metadata table (%s, ...).',
        paste(sort(missing.taxa.samples)[1:5],collapse=', '),
        paste(sort(missing.map.samples)[1:5],collapse=', ')))
    }
    x <- x[intersect(rownames(x),rownames(m)),,drop=F]
    m <- droplevels(m[rownames(x),,drop=F])
}

# remove rare features (do by minimum prevalence or average prevalence)
print(dim(x))
if (dim(x)[2] > 1) x <- x[,colMeans(x > 0) >= opts$min_prevalence, drop = FALSE]

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

diff.tests = NULL;
# identify differentiated features
print(dim(m))
if (dim(m)[2] > 1)
try (diff.tests <- differentiation.test(x, m[,opts$column], alpha=2, parametric=FALSE),silent=T)


#if(is.null(opts$which_taxa)){
if (!is.null(diff.tests)) {
    if(is.null(opts$nplot)){
        hit.ix <- which(diff.tests$qvalues <= opts$alpha)
    } else {
        hit.ix <- order(diff.tests$pvalues)[1:min(opts$nplot,length(diff.tests$pvalues))]
    }
} else {
    hit.ix <- match(opts$which_taxa, colnames(x))
}

if(length(hit.ix) == 0){
    print('No hits found.')
} else {
    cat("There were",length(hit.ix),"taxa significant at FDR of",opts$alpha,'\n')
}

# save hits
if (length(hit.ix)) write.differentiation.test.results(diff.tests, filename=sprintf('%s/test_results.txt',opts$outdir))

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

plot.beeswarm.dt(cols, opts$x_axis_label, hit.ix, env, opts$outdir)

