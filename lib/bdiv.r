# regress out any additional distance between samples
# that is explained by the given features
# regression is performed in order of features provided,
# not sure if this makes a difference.
#
#
# d is distance matrix/object, m is metadata table
"residual.distance" <- function(d,m){
	d <- as.matrix(d)
	if(is.null(dim(m))) m <- data.frame(m=m)
	if(mean(rownames(d) %in% rownames(m)) < 1) stop('rownames in d that are not in m.')
	m <- droplevels(m[rownames(d),,drop=F])
	
	# get actual distances and model matrix of between-group indicators as vectors
	v <- dist.as.named.vector(d)
	mm <- between.group.indicator.matrix(m)
	fit <- lm(v ~ mm)
	vr <- resid(fit)
	dr <- as.matrix(d)
	dr[upper.tri(dr)] <- vr
	dr[lower.tri(dr)] <- vr
	rownames(dr) <- rownames(m)
	colnames(dr) <- rownames(m)
	return(dr)
}

"between.group.indicator.matrix" <- function(m){
	# a binary model matrix with a column for every type of pair of groups
	# entry of a given column is 1 if that distance in v
	# is a distance between the two groups represented by that column
	mm <- NULL
	
	# loop through columns of m
	# initialize a zero-distance matrix for precomputation only
	d.0 <- matrix(0,nrow=nrow(m),ncol=nrow(m))
	rownames(d.0) <- rownames(m); colnames(d.0) <- rownames(m)
	for(k in 1:ncol(m)){
		if(is.numeric(m[,k])){
			# real=valued case: convert real-valued column to distance
			d.ij <- abs(outer(m[,k],m[,k],'-'))
			mm <- cbind(mm,as.numeric(as.dist(d.ij)))
			colnames(mm)[ncol(mm)] <- colnames(m)[k]
		} else {
			# discrete case: loop through all pairwise combinations of factor levels
			# discrete case
			group.names <- levels(droplevels(as.factor(m[,k])))
			for(i in 1:(length(group.names)-1)){
				for(j in (i+1):length(group.names)){
					if(i != j){
						d.ij <- d.0
						i.ix <- m[,k] == group.names[i]
						j.ix <- m[,k] == group.names[j]

						# contains 1 if comparison in v.ids is a comparison between
						# groups i and j
						d.ij[i.ix,j.ix] <- 1
						d.ij[j.ix,i.ix] <- 1
						mm <- cbind(mm,as.numeric(as.dist(d.ij)))
						colnames(mm)[ncol(mm)] <- paste(colnames(m)[k],'_',group.names[i],'_v_',group.names[j],sep='')
					}
				}
			}
		}

	}
	return(mm)
}

"dist.as.named.vector" <- function(d){
	# convert distances to vector
	v <- as.numeric(as.dist(d))
	# get names of samples being compared in v
	v.ids <- t(combn(rownames(m),2))
	names(v) <- paste(v.ids[,1],v.ids[,2],sep='_v_')
	return(v)
}