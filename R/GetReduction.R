#' GetReduction
#'
#' Extract the reduction from a Seurat object
#'
#' @param object A Seurat object from which to extract the reduction 
#' @param reduction The name of the reduction to extract
#' @param dims The dimension(s) of the reduction to extract. Defaults to the first two
#' @return A matrix w/ reduction coordinates of the specified dims, per cell.
#' @export
#' @examples
#' GetReduction(seuratobj, "umap", 1)
#' GetReduction(seuratobj, "pca", 1:2)

GetReduction <- function(object, reduction, dims){
  # Input validation
  stopifnot(class(object) == "Seurat")
  stopifnot(min(dims) >= 1)
  stopifnot(max(dims) <= ncol(object@reductions[[reduction]]@cell.embeddings))
  
  # Return reduction
  reduc_arr = object@reductions[[reduction]]@cell.embeddings[,dims]
  return(as.numeric(reduc_arr))
}

PCA1 = function(object) GetReduction(object, "pca", 1)
PCA2 = function(object) GetReduction(object, "pca", 2)

UMAP1 = function(object) GetReduction(object, "umap", 1)
UMAP2 = function(object) GetReduction(object, "umap", 2)

TSNE1 = function(object) GetReduction(object, "tsne", 1)
TSNE2 = function(object) GetReduction(object, "tsne", 2)


