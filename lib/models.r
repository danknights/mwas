# Generic tuning method for 
"tune.mwas" <- function(...) UseMethod("tune.mwas")
"predict.mwas" <- function(...) UseMethod("predict.mwas")


"mwas.train" <- function(x, y,
						distmat=distmat,
						distance.metric=distance.metric,
						model.type=model.type,
						tuning.parameters=tuning.parameters
				){
				
				
}


"predict.mwas.rf" <- function(){

}


"predict.mwas.svm" <- function(){

}

"predict.mwas.knn" <- function(){

}

"predict.mwas.glmnet" <- function(){

}

"predict.mwas.null" <- function(){

}

