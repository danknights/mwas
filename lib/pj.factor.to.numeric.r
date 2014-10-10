factor.to.numeric<-function(x)
{
	x <- factor(x)
   	levels(x) <- 1:length(levels(x))
   	x <- as.numeric(x)
   	return(x)
}
