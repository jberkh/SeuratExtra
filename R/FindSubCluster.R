#' FindSubCluster
#'
#' Convenience wrapper for Seurat::FindSubCluster
#'
#' @param object A Seurat object to find subclusters in
#' @param cluster The cluster to be sub-clustered
#' @param graph.name The name of the graph to use for the clustering algorithm. Defaults to RNA_snn for convenience 
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @return A Seurat object with subclusters stored in Idents(object) instead of a meta.data columns
#' @export
#' @examples
#' FindSubCluster(seuratobj, 'CD8+')
#' FindMarkers(FindSubCluster(seuratobj, "CD4+"), ident.1 = "CD4+_0", ident.2 = "CD4+_1")


FindSubCluster <- function(object, cluster, ...) {
  # Input validation
  stopifnot(class(object) == "Seurat")
  stopifnot(cluster %in% levels(object))
  
  # Preparing function cal arguments
  # Capture all arguments in a list
  args <- list(...)  
  # Put user inputs in args
  args$object = object
  args$cluster = cluster
  # If not explicitly defined, use RNA_snn as default
  if (!"graph.name" %in% names(args)) args$graph.name = "RNA_snn"
  # Override subcluster.name since subclusters are stored in Idents()
  args$subcluster.name = "subcluster"
  
  # Get subcluster output
  ids = do.call(Seurat::FindSubCluster, args)$subcluster
  
  # Fix factor levels
  # Create all factor levels w/ loop
  lvls = list(); i = 0
  for (lvl in levels(object)) lvls[[i<-i+1]] = c(lvl, paste0(lvl, "_", 0:99))
  # Retain existing levels from generated; Result == ordered
  lvls = unlist(lvls)[unlist(lvls) %in% unique(ids)]
  # Set object Idents() and return
  Idents(object) = factor(ids, levels = lvls)
  return(object)
}