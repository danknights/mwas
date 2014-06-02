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
        help="Input taxon table [required]."),
    make_option(c("-m","--map_fp"), type="character",default=NULL,
        help="Mapping file (optional). If provided, only samples in both taxon table and mapping file will be plotted."),
    make_option(c("-w", "--which_taxa"), type="character", default=NULL,
        help="Comma-separated list of taxa to plot [default: plot top --nplot taxa]"),
    make_option(c("-s", "--shorten_taxa"),action='store_true',default=FALSE,
        help="Shorten taxonomy names to lowest defined level. [default: %default]"),
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
x <- t(read.table(opts$input_fp,sep='\t',head=T,row=1,check=F))
if(!is.null(opts$map_fp)){
	m <- read.table(opts$map_fp,sep='\t',head=T,row=1,check=F,comment='')
	x <- x[intersect(rownames(x),rownames(m)),]
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

pc <- cmdscale(vegdist(x),k=5)


# plots
fp <- sprintf('%s/gradients.pdf',opts$outdir)
pdf(fp,width=11,height=3.75)
par(mfrow=c(1,3))
for(i in seq_along(taxon.names)){
# 	fp <- sprintf('%s/gradient-%s.pdf',opts$outdir,taxon.names[i])
# 	pdf(fp,width=11,height=3.75)
# 	par(mfrow=c(1,3))
	combs <- combn(1:3,2)
	for(j in 1:ncol(combs)){
		show.gradients(x[,taxon.names[i]], pc[,combs[,j]], incl.legend=TRUE,pt.alpha='AA',
			axis.labels=sprintf('PC%d',combs[,j]),
			title.text=sprintf('%s - PC%d v PC%d',taxon.names[i],combs[1,j],combs[2,j]))
	}
}
dev.off()


