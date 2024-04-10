## Attempt to standardize some network-formation functions.

#' Given a matrix of scores (bit score, e-value, etc), create an adjacency graph.
#'
#' I am hoping to use this as the starting point for a generic network generator.
#' In its current form it takes a matrix of pairwise scores and
#' generates an adjacency graph of those scores.
#'
#' @param scores tsv or matrix of scores with column and row names containing IDs.
#' @param metadata Currently unused, but intended to provide a
#'  starting point for annotating the resulting adjacency network.
#'  When implemented, it should make use of the annotate_network()
#'  functions which follow.
#' @param type Currently I only know of networks which use
#'  correlation, distance, and distcor matrices of the original
#'  scores; but I suspect a cursory glance at the WGCNA documentation
#'  will teach me that there are many more possibilities.
#' @param simplify Return a simplified matrix without loops and
#'  redundancies?
#' @param mode Network type to create, I don't yet understand the implications
#'  of changing this.
#' @param weighted Add weights to the nodes?  I also don't yet
#'  understand what happens when you mess with this.
#' @param diag Include the matrix-diagonal nodes?  I do not know when
#'  one would want these.
#' @return igraph adjacency network.
network_from_matrix <- function(scores, metadata = NULL, type = "distcor", simplify = TRUE,
                                mode = "undirected", weighted = TRUE, diag = FALSE) {
  if ("character" %in% class(scores)) {
    scores <- as.matrix(read.table(scores, header = TRUE, sep = "\t",
                                   row.names = 1))
  } else {
    scores <- as.matrix(scores)
  }

  input <- scores
  if (type == "cor") {
    mesg("Calculating correlation matrix.")
    input <- as.matrix(cor(input))
  } else if (type == "dist") {
    mesg("Calculating distance matrix.")
    input <- as.matrix(stats::as.dist(input))
  } else if (type == "distcor") {
    mesg("Calculating correlation matrix.")
    cor_mtrx <- cor(input)
    mesg("Calculating distance matrix of correlations.")
    input <- as.matrix(stats::as.dist(cor_mtrx))
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

#' Use grep to add a vector of annotations/colors to a network.
#'
#' The igraph syntaxes are a little clunky, but the set_attr()
#' functions mostly make sense.
#'
#' @param network Input network
#' @param names set of node-names to which to add annotations.
#' @param color Color to attach to the added annotation.
#' @param default Set a default annotation for this name to all nodes.
#' @param annot_name Annotation name to attach to the nodes.
#' @param annot_value and the associated value.
#' @return a new network!
annotate_network <- function(network, names, color = NULL, default = NULL,
                             annot_name = "type", annot_value = "high") {
  net_names <- igraph::vertex_attr(network, name = "name")
  if (!is.null(default)) {
    network <- igraph::set_vertex_attr(graph = network, name = annot_name,
                                       value = default)
  }
  for (name in names) {
    wanted_names <- grepl(x = net_names, pattern = name)
    mesg("Network name: ", name, " was found ", sum(wanted_names), " times.")
    network <- igraph::set_vertex_attr(graph = network, name = annot_name,
                                       index = wanted_names, value = annot_value)
    if (!is.null(color)) {
      network <- igraph::set_vertex_attr(graph = network, name = "color",
                                         index = wanted_names, value = color)
    }
  }
  return(network)
}

#' A version of annotate_network, but which uses a dataframe as input.
#'
#' The annotate_network() function uses a vector of values, this
#' extends that logic to add every column of a dataframe.  I would
#' like to make this function a little more fun vis a vis abilities to
#' add colors and such.
#'
#' @param network input network.
#' @param df input dataframe, columns are the new metadata, rows are
#'  the node-strings to search on.
#' @param default Set a default?
annotate_network_df <- function(network, df, default = NULL) {
  new <- network
  net_names <- igraph::vertex_attr(network, name = "name")
  df <- as.data.frame(df)
  for (col in colnames(df)) {
    mesg("Starting annotation of ", col, ".")
    values <- df[[col]]
    names(values) <- rownames(df)
    if (!is.null(default)) {
      new <- igraph::set_vertex_attr(graph = new, name = col, value = default)
    }
    possibilities <- as.factor(values)
    for (pos in levels(possibilities)) {
      mesg("Setting vertices to: ", pos, ".")
      pos_idx <- values == pos
      pos_names <- names(values)[pos_idx]
      pos_positive <- net_names %in% pos_names
      positive_vertices <- igraph::V(network)[pos_positive]
      new <- igraph::set_vertex_attr(graph = new, name = col,
                                     index = positive_vertices, value = pos)
    }
  }
  return(new)
}

annotate_network_nodes <- function(network, df, column = "assemblyxls", col_number = 5,
                                   column_name = "interpfamvalue") {
  entries <- rownames(df)
  net_names <- igraph::vertex_attr(network, name = "name")
  nodes_annotated <- 0
  for (e in seq_len(length(entries))) {
    entry <- entries[e]
    gene_annotation_file <- df[e, column]
    gene_annotations <- try(sm(extract_metadata(gene_annotation_file,
                                                fill = "")), silent=TRUE)
    if ("try-error" %in% class(gene_annotations)) {
      next
    }
    pfam_names <- as.data.frame(stringr::str_split_fixed(gene_annotations[[column_name]],
                                                         ",", col_number))[[col_number]]
    names(pfam_names) <- rownames(gene_annotations)

    defined_names <- pfam_names != ""
    pfam_defined <- pfam_names[defined_names]
    entries_found <- 0
    nodes_found <- 0
    for (d in seq_len(length(pfam_defined))) {
      defined_name <- names(pfam_defined)[d]
      defined_value <- pfam_defined[d]
      pfam_nodes <- net_names == defined_name
      times_found <- sum(pfam_nodes)
      if (times_found > 0) {
        entries_found <- entries_found + 1
        nodes_found <- nodes_found + times_found
        nodes_annotated <- nodes_annotated + times_found
        pfam_vertices <- igraph::V(network)[pfam_nodes]
        network <- igraph::set_vertex_attr(
          graph = network, name = column_name,
          index = pfam_vertices, value = defined_value)
      }
      mesg("Finished searching sample: ", entry, ", found ", entries_found, " genes across ",
           nodes_found, " network vertices.")
    }
  }
  message("Finished iterating over annotations, annotated: ", nodes_annotated,
          " vertices out of ", length(net_names), ".")
  return(network)
}

#' Exclude nodes from a network which are not well connected.
#'
#' @param network input network to prune.
#' @param min_weight Minimum acceptable weight.
#' @param min_connectivity Minimum number of nodes to which to be connected.
#' @return A hopefully smaller, but not too small network.
prune_network <- function(network, min_weight = 0.4, min_connectivity = 1) {
  start_vertices <- length(igraph::V(network))
  start_edges <- length(igraph::E(network))
  low_nodes <- igraph::edge_attr(network, name = "weight") <= min_weight
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
  mesg("Network pruning vertices, start: ", start_vertices, ", end: ",
       end_vertices, ", delta: ", delta_vertices, ".")
  if (start_edges > 0) {
    mesg("Network pruning edges, start: ", start_edges, ", end: ",
         end_edges, ", delta: ", delta_edges, ".")
  }
  return(pruned_vertices)
}

#' Just putting down an example from Alejandro's work and
#' https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
#' so that I have someplace to remember the general path wgcna takes and
#' as a starting point to explore further.
#'
#' @param expt Input expressionset.
wgcna_network <- function(expt) {
  chosen_power <- 8
  ## WGCNA calls cor() without specifying its own namespace, so overwrite cor for the moment.
  cor <- WGCNA::cor
  start_exprs <- t(exprs(expt))
  l2input <- t(exprs(normalize_expt(expt, transform = "log2")))

  initial_modules <- WGCNA::blockwiseModules(
    start_exprs, maxBlockSize = 11000, TOMType = "signed",
    power = chosen_power, mergeCutHeight = 0.25, numericLabels = FALSE,
    verbose = 3)
  cor <- stats::cor

  initial_eigen <- initial_modules[["MEs"]]
  network_colors <- WGCNA::labels2colors(initial_modules[["colors"]])
  WGCNA::plotDendroAndColors(
    initial_modules[["dendrograms"]][[1]],
    network_colors[initial_modules[["blockGenes"]][[1]]],
    "Modules",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05)

  meta_numeric <- data.frame(
    "cf_numeric" = as.numeric(as.factor(pData(l2input)[["finaloutcome"]])),
    "visit_numeric" = as.numeric(as.factor(pData(l2input)[["visitnumber"]])))
  rownames(meta_numeric) <- rownames(pData(l2input))

  meta_factors <- pData(l2input)[, c("finaloutcome", "visitnumber")]
  meta_eigen <- merge(initial_eigen, meta_factors, by = "row.names")
  rownames(meta_eigen) <- meta_eigen[["Row.names"]]
  meta_eigen[["Row.names"]] <- NULL

  module_trait_corr <- stats::cor(initial_eigen, meta_numeric, use = "p")
  module_trait_corr
  module_trait_pvalues <- WGCNA::corPvalueStudent(module_trait_corr, nrow(meta_numeric))
  module_trait_pvalues

  wanted <- initial_modules[["colors"]] == "turqoise" | initial_modules[["colors"]] == "pink"
  sum(wanted)
  interesting_genes <- names(initial_modules[["colors"]])[wanted]

  fData(l2input)[interesting_genes, "hgnc_symbol"]


  ## Note that we can do similarity matrices on the samples too in order to get
  ## dendrograms which may get interesting groups of samples?
  not_grey <- initial_modules[["colors"]] != "grey"
  not_grey_exprs <- t(exprs(l2input))[, not_grey]
  dim(not_grey_exprs)
  not_grey_genes <- colnames(not_grey_exprs)
  dist_tom <- 1 - WGCNA::TOMsimilarityFromExpr(
    not_grey_exprs,
    power = chosen_power)
  colnames(dist_tom) <- not_grey_genes
  rownames(dist_tom) <- colnames(dist_tom)

  similarity_cluster <- flashClust::flashClust(
    stats::as.dist(dist_tom), method = "average")

  WGCNA::plotDendroAndColors(
    similarity_cluster,
    colors = initial_modules[["colors"]][not_grey_genes],
    "ME",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    cex.colorLabels = 1.5,
    main = ("Cluster Dendrogram (WGCNA)"))

}

## EOF
