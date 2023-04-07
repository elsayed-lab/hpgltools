## Some functions to help with classification

generate_nn_groups <- function(mtrx, resolution = 1, k = 10, type = "snn",
                                full = TRUE, merge_to = NULL, ...) {
  params <- SNNGraphParam(k = k, ...)
  if (type == "knn") {
    params <- KNNGraphParam(k = k)
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
