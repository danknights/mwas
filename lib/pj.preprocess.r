preprocess(otu, minOTUInSamples=.001)
{
	source("util.r") 
	source("I-methods.r")
		
	otu <- (remove.nonoverlapping.samples(map = map, otus = otus))$otus

	otu <- sweep(otu, 1, rowSums(otu), '/')
	otu <- otu[, colMeans(otu) > minOTUInSamples, drop=FALSE]
	otu <- asin(sqrt(otu))
	ret <- collapse.by.correlation(otu, .95)
	otu <- otu[, ret$reps]

	return(otu)
}