"visualize.mwas" <- function(data, params){
  # Visualization
  # - beeswarm
  # - 
  #
  #
  
}

"plot.beeswarm.dt" <- function(cols, x_axis_label, hit.ix, env, outdir) {
	require(beeswarm)
	cols <-  c(brewer.pal(9,'Set1'),brewer.pal(9,'Pastel1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[-6]
	cols[1:2] <- cols[2:1]
	cols <- sprintf('%sbb',cols)
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
		dev.off()
	}
}
