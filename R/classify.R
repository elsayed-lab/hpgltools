## Some functions to help with classification

#' Given an n-dimensional matrix, try some KNN-esque clustering on it.
#'
#' I want some functions to help me understand clustering.  This is a
#' first pass at that goal.
#'
#' @param mtrx Matrix to cluster, usually 2d from a point plot.
#' @param resolution Used after cluster generation for making neighbor
#'  groups.
#' @param k Used during cluster generation.
#' @param type Define the type of clustering to perform, currently
#'  only KNN/SNN
#' @param full Get the full set of metrics from bluster.
#' @param merge_to Use the neighborhood collapse function to set a
#'  hard ceiling on the number of clusters in the final result.
generate_nn_groups <- function(mtrx, resolution = 1, k = 10, type = "snn",
                                full = TRUE, merge_to = NULL, ...) {
  params <- bluster::SNNGraphParam(k = k, ...)
  if (type == "knn") {
    params <- bluster::KNNGraphParam(k = k)
  }
  clusters <- bluster::clusterRows(mtrx, params, full = full)
  merged <- NULL
  groups <- as.factor(paste0("g", clusters$clusters))
  message("After clustering, there are: ", length(levels(groups)), " groups.")
  if (!is.null(merge_to)) {
    merged <- bluster::mergeCommunities(clusters[["objects"]][["graph"]],
                                        clusters[["clusters"]],
                                        steps = merge_to)
  }
  merged_groups <- as.factor(paste0("m", merged$clusters))
  message("After merging, there are: ", length(levels(merged_groups)), " groups.")

  retlist <- list(
    "clusters" = clusters,
    "merged" = merged,
    "start_groups" = groups,
    "merged_groups" = merged_groups)
  return(retlist)
}
