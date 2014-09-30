"mwas.feature.scores" <- function(x, y, feature.ids, filename='feature_importance_scores.txt', outdir='.'){
  result <- rf.out.of.bag(x, y)
  save.rf.results.importances(result, feature.ids, filename, outdir)
  
}