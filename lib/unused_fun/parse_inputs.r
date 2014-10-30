# Universal input reader and parser. Reads in all types of input files supported by 
# MWAS, does sanity and error checking, intersection, and empty line dropping.
# Automatically detects correct input format for QIIME/BIOM-type inputs. 
# Returns an object with the processed objects. Note that feature stats table
# is ordered/subset based on the result of the final feature list, not vice versa.

"load.inputs" <- function(featureTable=NULL, mapFile=NULL, labelCols=NULL, distMat=NULL, PCoAfile=NULL, featureStats=NULL, dropEmpty=T) {
	# Read in all specified files
	fTab = ifelse(is.null(featureTable), NULL, read.qiime.table(featureTable));
	map = ifelse(is.null(mapFile), NULL, load.qiime.mapping.file(mapFile));
	dist = ifelse(is.null(distMat), NULL, load.qiime.distance.matrix(distMat));
	pcoa = ifelse(is.null(PCoAfile), NULL, load.qiime.pcoa.file(PCoAfile)); # new
	fS = ifelse(is.null(featureStats), NULL, load.qiime.feature.stats(featureStats)); 
	isct = unique(c(rownames(fTab), rownames(map), rownames(dist), rownames(pcoa)));
	if (is.null(fTab) && is.null(map) && is.null(dist) && is.null(pcoa)) # no input
		stop("No input file path is valid.");
	
	# Process sample name intersection
	if (!is.null(fTab)) isct = intersect(isct, rownames(fTab));
	if (!is.null(map)) isct = intersect(isct, rownames(map));
	if (!is.null(dist)) isct = intersect(isct, rownames(dist));
	if (!is.null(pcoa)) isct = intersect(isct, rownames(pcoa));
	if (length(isct) < 1) stop("Intersection empty; check sample names across files.");
	
	if (!is.null(fS) && !is.null(fTab)) # order/subset fS based on fTab's features
		fS = fS[match(intersect(rownames(fS),colnames(fTab)),rownames(fS)),];
	
	# Re-order all inputs to match the list of samples in common (subsetting)
	if (!is.null(fTab)) fTab = fTab[match(isct, rownames(fTab)),];
	if (!is.null(map)) map = map[match(isct, rownames(map)),];
	if (!is.null(dist)) dist = dist[match(isct, rownames(dist)),];
	if (!is.null(pcoa)) pcoa = pcoa[match(isct, rownames(pcoa)),];
	
	return (list(fTab=fTab, map=map, dist=dist, pcoa=pcoa, fS=fS, isct=isct));

}