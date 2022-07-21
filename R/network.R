## Attempt to standardize some network-formation functions.

network_from_matrix <- function(score_file, metadata = NULL, type = "distcor", simplify = TRUE,
                                mode = "undirected", weighted = TRUE, diag = FALSE) {
  scores <- as.matrix(read.table(score_file, header = TRUE, sep = "\t", row.names = 1))
  input <- scores
  if (type == "cor") {
    message("Calculating correlation matrix.")
    input <- as.matrix(cor(input))
  } else if (type == "dist") {
    message("Calculating distance matrix.")
    input <- as.matrix(as.dist(input))
  } else if (type == "distcor") {
    message("Calculating correlation matrix.")
    cor_mtrx <- cor(input)
    message("Calculating distance matrix of correlations.")
    input <- as.matrix(as.dist(cor_mtrx))
  }

  ## initial <- igraph::graph.adjacency(input, mode="undirected", weighted=TRUE, diag=FALSE)
  initial <- igraph::graph.adjacency(input, mode = mode,
                                     weighted = weighted, diag = diag)
  if (isTRUE(simplify)) {
    initial <- igraph::simplify(initial, remove.multiple = TRUE,
                                remove.loops = TRUE)
  }
  return(initial)
}

annotate_network <- function(network, names, color = NULL, default = NULL,
                             annot_name = "type", annot_value = "high") {
  net_names <- igraph::V(network)$name
  if (!is.null(default)) {
    network <- igraph::set_vertex_attr(graph = network, name = annot_name,
                                       value = default)
  }
  for (name in names) {
    wanted_names <- grepl(x = net_names, pattern = name)
    message("Network name: ", name, " was found ", sum(wanted_names), " times.")
    network <- igraph::set_vertex_attr(graph = network, name = annot_name,
                                       index = wanted_names, value = annot_value)
    if (!is.null(color)) {
      network <- igraph::set_vertex_attr(graph = network, name = "color",
                                         index = wanted_names, value = color)
    }
  }
  return(network)
}

annotate_network_df <- function(network, df, default = NULL) {
  new <- network
  net_names <- igraph::V(network)$name
  df <- as.data.frame(df)
  for (col in colnames(df)) {
    message("Starting annotation of ", col, ".")
    values <- df[[col]]
    names(values) <- rownames(df)
    if (!is.null(default)) {
      new <- igraph::set_vertex_attr(graph = new, name = col, value = default)
    }
    possibilities <- as.factor(values)
    for (pos in levels(possibilities)) {
      message("Setting vertices to: ", pos, ".")
      pos_idx <- values == pos
      pos_names <- names(values)[pos_idx]
      pos_positive <- net_names %in% pos_names
      positive_vertices <- igraph::V(network)[pos_positive]
      new <- igraph::set_vertex_attr(graph = new, name = pos,
                                     index = positive_vertices, value = pos)
    }
  }
  return(new)
}

prune_network <- function(network, min_weight = 0.4, min_connectivity = 1) {
  start_vertices <- length(igraph::V(network))
  start_edges <- length(igraph::E(network))
  low_nodes <- igraph::E(network)$weight <= min_weight
  low_edges <- igraph::E(network)[low_nodes]
  pruned_edges <- igraph::delete_edges(network, low_edges)
  low_degree <- igraph::degree(pruned_edges) <= min_connectivity
  if (sum(low_degree) == start_vertices) {
    stop("Every vertex is removed by this minimum connectivity.")
  }
  pruned_vertices <- igraph::delete.vertices(pruned_edges, low_degree)
  end_vertices <- length(igraph::V(pruned_vertices))
  end_edges <- length(igraph::E(pruned_vertices))
  delta_vertices <- start_vertices - end_vertices
  delta_edges <- start_edges - end_edges
  message("Network pruning vertices, start: ", start_vertices, ", end: ",
          end_vertices, ", delta: ", delta_vertices, ".")
  if (start_edges > 0) {
    message("Network pruning edges, start: ", start_edges, ", end: ",
            end_edges, ", delta: ", delta_edges, ".")
  }
  return(pruned_vertices)
}
