#' RandomSubset
#'
#' Creates a subset of a Seurat object
#'
#' @param object A Seurat object to find subset
#' @param by.factor The factor to determine subset size (1 - ∞)
#' @param group.by The name of the metadata column used for determining sample probabilities.
#' @export
#' @examples
#' obj = RandomSubset(obj, by.factor = 2, group.by = "subclusters")

RandomSubset = function(object, by.factor, group.by = NULL) {
  stopifnot(class(object) == "Seurat")
  stopifnot(by.factor > 1)
  
  # Calculate subsample probabilities based on group.by
  # Deliberate separate if statements
  if (is.null(group.by)) {
    prob = rep(1, ncol(object))
  }
  if (!is.null(group.by) & group.by == "ident") {
    object$ident = Idents(object)
  }
  if (!is.null(group.by)) {
    group.by = factor(object@meta.data[,group.by])
    prob = mapvalues(
      group.by, 
      levels(group.by), 
      table(group.by)**(-1))
    prob = as.numeric(prob) / length(levels(group.by))
  }
  
  # Subsample and return
  n = ncol(object)
  n_sub = ceiling(n / by.factor)
  return(object[,sample(n, n_sub, prob = prob)])
 
}
