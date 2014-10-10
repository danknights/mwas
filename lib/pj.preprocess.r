preprocess(otu, minOTUInSamples=.001)
{
	source("util.r") 
	otu <- sweep(otu, 1, rowSums(otu), '/')
	otu <- otu[, colMeans(otu) > minOTUInSamples, drop=FALSE]
	otu <- asin(sqrt(otu))
	ret <- collapse.by.correlation(otu, .95)
	otu <- otu[, ret$reps]

	return(otu)
}