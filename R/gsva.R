#' Change gene IDs to the format expected by gsva using an orgdb.
#'
#' Though it is possible to use gsva without ENTREZ IDs, it is not trivial.
#' This function attempts to ensure that the IDs in one's expressionset are
#' therefore entrez IDs. It is possible that this function is at least partially
#' redundant with other functions in this package and should be replaced.
#'
#' @param ids Vector of IDS to modify.
#' @param from Change from this format.
#' @param to Change to this format.
#' @param orgdb Using this orgdb instance.
#' @return New vector of ENTREZ IDs.
#' @seealso [AnnotationDbi]
#' @export
convert_ids <- function(ids, from = "ENSEMBL", to = "ENTREZID", orgdb = "org.Hs.eg.db") {
  lib_result <- sm(requireNamespace(orgdb))
  att_result <- sm(try(attachNamespace(orgdb), silent = TRUE))
  new_ids <- sm(AnnotationDbi::select(x = get0(orgdb),
                                      keys = ids,
                                      keytype = from,
                                      columns = to))
  new_idx <- complete.cases(new_ids)
  new_ids <- new_ids[new_idx, ]
  message("Before conversion: ", length(ids),
          ", after conversion: ", length(rownames(new_ids)), ".")
  return(new_ids)
}

#' Extract the GeneSets corresponding to the provided name(s).
#'
#' Many of the likely GSCs contain far more gene sets than one actually wants to
#' deal with.  This will subset them according to a the desired 'requests'.
#'
#' @param sig_data The pile of GeneSets, probably from GSVAdata.
#' @param requests Character list of sources to keep.
#' @return Whatever GeneSets remain.
#' @export
get_gsvadb_names <- function(sig_data, requests = NULL) {
  requests <- toupper(requests)
  categories <- names(sig_data)
  prefixes <- gsub(pattern = "^(.*?)_.*", replacement = "\\1",
                   x = categories)
  tab <- table(prefixes)
  if (is.null(requests)) {
    new_order <- order(as.numeric(tab), decreasing = TRUE)
    tab <- tab[new_order]
    return(tab)
  }

  kept_requests <- requests %in% names(tab)
  requests <- requests[kept_requests]
  message("Subsetting the datasets to extract the: ", toString(requests), " data.")

  kept_idx <- rep(FALSE, length(names(sig_data)))
  for (kept in requests) {
    keepers <- grepl(x = names(sig_data), pattern = glue("^{kept}"))
    kept_idx <- kept_idx | keepers
  }

  message("After subsetting, ", length(remaining), " entries remain.")
  remaining <- sig_data[kept_idx]
  return(remaining)
}

#' Create dataframe which gets the maximum within group mean gsva score for each gene set
#'
#' @param gsva_scores Result from simple_gsva()
#' @param groups list of groups for which to calculate the means
#' @param keep_single Keep categories with only 1 element.
#' @param method mean or median?
#' @return dataframe containing max_gsva_score, and within group means for gsva scores
#' @seealso [simple_gsva()]
#' @export
get_group_gsva_means <- function(gsva_scores, groups, keep_single = TRUE, method = "mean") {
  start_gsva_result <- exprs(gsva_scores)

  groupMeans <- data.frame(row.names = rownames(start_gsva_result))
  groupInd <- list()
  ## Added a little logic in case there are categories with <= 1 column.
  for (group in groups) {
    ## get columns for that group
    ind <- pData(gsva_scores)[["condition"]] == group
    if (sum(ind) < 1) {
      next
    } else if (sum(ind) == 1) {
      if (isTRUE(keep_single)) {
        groupMeans[[group]] <- as.data.frame(start_gsva_result[, ind])
        colnames(groupMeans[[group]]) <- colnames(pData(gsva_scores))[ind]
      } else {
        next
      }
    } else {
      subset <- start_gsva_result[, ind]
      if (method == "mean") {
        groupMeans[[group]] <- abs(rowMeans(subset))
      } else if (method == "median") {
        groupMeans[[group]] <- abs(Biobase::rowMedians(subset))
      } else {
        stop("I do not know this method: ", method, ".")
      }
    }
    groupInd[[group]] <- ind
  }
  return(list("Means" = groupMeans, "Index" = groupInd))
}

#' Attempt to score the results from simple_gsva()
#'
#' This function uses a couple of methods to try to get an idea of whether the
#' results from gsva are actually interesting.  It does so via the following
#' methods:
#'   1.  Use limma on the expressionset returned by simple_gsva(), this might
#' provide an idea of if there are changing signatures among the sample types.
#'   2.  Perform a simplified likelihood estimate to get a sense of the
#' significant categories.
#'
#' @param gsva_result Result from simple_gsva()
#' @param cutoff Significance cutoff
#' @param excel Excel file to write the results.
#' @param model_batch Add batch to limma's model.
#' @param factor_column When extracting significance information, use this
#'  metadata factor.
#' @param factor Use this metadata factor as the reference.
#' @param label_size Used to make the category names easier to read at the expense
#'  of dropping some.
#' @param col_margin Attempt to make heatmaps fit better on the screen with this and...
#' @param row_margin this parameter
#' @param type Either mean or median of the scores to return.
#' @return List containing the gsva results, limma results, scores, some plots, etc.
#' @seealso [score_gsva_likelihoods()] [get_group_gsva_means()] [limma_pairwise()]
#'  [simple_gsva()]
#' @export
get_sig_gsva_categories <- function(gsva_result, cutoff = 0.95, excel = "excel/gsva_subset.xlsx",
                                    model_batch = FALSE, factor_column = "condition", factor = NULL,
                                    label_size = NULL, col_margin = 6, row_margin = 12,
                                    type = "mean") {
  gsva_scores <- gsva_result[["expt"]]

  ## FIXME: If one uses factor_column in this current function, that will likely lead to
  ## incorrect results because limma is using the 'condition' metadata factor; but
  ## median_by_factor() is using this new column 'factor_column'.  As
  ## a result our mean gsva scores will no longer have any connection
  ## to the results from limma.
  ## I think this may be trivially fixed though? ...
  if (factor_column != "condition") {
    gsva_scores <- set_expt_conditions(gsva_scores, fact = factor_column)
  }

  ## Use limma on the gsva result
  gsva_limma <- limma_pairwise(gsva_scores, model_batch = model_batch,
                               which_voom = "none")

  ## Combine gsva max(scores) with limma results
  ### get gsva within group means
  groups <- levels(gsva_scores[["conditions"]])
  gsva_score_means <- median_by_factor(data = gsva_scores, fact = factor_column, fun = type)
  ## gsva_score_means <- get_group_gsva_means(gsva_scores, groups)
  num_den_string <- strsplit(x = names(gsva_limma[["all_tables"]]), split = "_vs_")

  for (t in seq_along(gsva_limma[["all_tables"]])) {
    table_name <- names(gsva_limma[["all_table"]])[t]
    table <- gsva_limma[["all_tables"]][[t]]
    contrast <- num_den_string[[t]]
    ## get means from gsva_score_means for each contrast
    numerator <- contrast[1]
    denominator <- contrast[2]
    ## first <- gsva_score_means[["Index"]][numerator][[1]]
    ## second <- gsva_score_means[["Index"]][denominator][[1]]
    numerator_samples <- gsva_score_means[["indexes"]][[numerator]]
    if (is.null(numerator_samples)) {
      next
    }
    denominator_samples <- gsva_score_means[["indexes"]][[denominator]]
    if (is.null(denominator_samples)) {
      next
    }
    ## get maximum value of group means in each contrast
    ## maxs <- apply(exprs(gsva_scores)[, first | second], 1, max)
    max_values <- apply(exprs(gsva_scores)[, numerator_samples | denominator_samples], 1, max)
    table[["gsva_score_max"]] <- max_values
    varname1 <- paste0("Mean_", numerator)
    varname2 <- paste0("Mean_", denominator)
    ## table[[varname1]] <- gsva_score_means[["Means"]][[contrasts[1]]]
    table[[varname1]] <- gsva_score_means[["medians"]][[numerator]]
    ## table[[varname2]] <- gsva_score_means[["Means"]][[contrasts[2]]]
    table[[varname2]] <- gsva_score_means[["medians"]][[denominator]]
    gsva_limma[["all_tables"]][[t]] <- table
  }

  ## FIXME: This is not likely needed anymore.
  gsva_eset <- gsva_scores[["expressionset"]]
  ## Go from highest to lowest score, using the first sample as a guide.
  values <- as.data.frame(exprs(gsva_eset))
  annot <- fData(gsva_eset)
  meta <- pData(gsva_eset)

  ## Choose the reference factor
  clevels <- levels(as.factor(meta[[factor_column]]))
  fact <- factor
  if (is.null(factor)) {
    fact <- clevels[1]
  }

  ## Copy the gsva expressionset and use that to pull the 'significant' entries.
  subset_eset <- gsva_eset
  gl <- score_gsva_likelihoods(gsva_result, factor = fact, label_size = label_size)
  likelihoods <- gl[["likelihoods"]]
  keep_idx <- likelihoods[[fact]] >= cutoff
  scored_ht <- subset_table <- scored_ht_plot <- NULL
  if (sum(keep_idx) > 1) {
    subset_eset <- subset_eset[keep_idx, ]
    jet_colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                                "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    tmp_file <- tmpmd5file(pattern = "heat", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    if (is.null(label_size)) {
      scored_ht <- heatmap.3(exprs(subset_eset), trace = "none", col = jet_colors,
                             margins = c(col_margin, row_margin))
    } else {
      scored_ht <- heatmap.3(exprs(subset_eset), trace = "none", col = jet_colors,
                             margins = c(col_margin, row_margin),
                             cexRow = label_size)
    }
    scored_ht_plot <- grDevices::recordPlot()
    dev.off()
    removed <- suppressWarnings(file.remove(tmp_file))
    removed <- unlink(dirname(tmp_file))
    subset_order <- rev(scored_ht[["rowInd"]])
    subset_rownames <- rownames(exprs(subset_eset))[subset_order]

    ## Here is a bizarre little fact: the rownames of fData(subset_mtrx)
    ## are not the same as the rownames of exprs(subset_mtrx)
    ## which AFAIK should not be possible, but clearly I was wrong.
    ## At least when I create an expressionset, the fData and exprs have the
    ## same rownames from beginning to end.
    ## I guess this does not really matter, since we can use the full annotation table.
    subset_tbl <- as.data.frame(exprs(subset_eset))
    order_column <- colnames(subset_tbl)[1]
    subset_table <- merge(annot, subset_tbl, by = "row.names")
    rownames(subset_table) <- subset_table[["Row.names"]]
    subset_table[["Row.names"]] <- NULL
    ## Set the table row order to the same as the hclust from the heatmap.
    subset_table <- subset_table[subset_rownames, ]
  } else {
    mesg("There are not many entries which pass the cutoff for significance.")
  }

  order_column <- colnames(values)[1]
  gsva_table <- merge(annot, values, by = "row.names")
  rownames(gsva_table) <- gsva_table[["Row.names"]]
  gsva_table[["Row.names"]] <- NULL
  reordered_gsva_idx <- order(gsva_table[[order_column]], decreasing = TRUE)
  gsva_table <- gsva_table[reordered_gsva_idx, ]

  gl_tbl <- as.data.frame(gl[["likelihoods"]])
  order_column <- colnames(gl_tbl)[1]
  likelihood_table <- merge(annot, gl_tbl, by = "row.names")
  rownames(likelihood_table) <- likelihood_table[["Row.names"]]
  likelihood_table[["Row.names"]] <- NULL
  likelihood_table_idx <- order(likelihood_table[[order_column]], decreasing = TRUE)
  likelihood_table <- likelihood_table[likelihood_table_idx, ]

  retlist <- list(
      ## Everything provided by simple_gsva()
      "input" = gsva_result,
      ## The table from simple_gsva merged with the annotations.
      "gsva_table" = gsva_table,
      ## Heatmap of gsva result.
      "raw_plot" = gl[["raw_plot"]],
      ## The result from score_gsva_likelihoods, which compares condition vs. others.
      "likelihood_table" = likelihood_table,
      ## Corresponding plot from score_gsva_likelihoods
      "score_plot" = gl[["likelihood_plot"]],
      ## The subset of gsva scores deemed 'significant' by score_gsva_likelihoods.
      "subset_table" = subset_table,
      ## The corresponding plot for the subset.
      "subset_plot" = scored_ht_plot,
      "scores" = gl,
      "score_pca" = gl[["pca"]][["plot"]],
      "subset_expt" = subset_eset,
      "gsva_limma" = gsva_limma)

  if (!is.null(excel)) {
    retlist[["excel"]] <- write_gsva(retlist, excel)
  }
  class(retlist) <- "gsva_sig"
  return(retlist)
}

#' Take a result from simple_gsva(), a list of gene IDs, and intersect them.
#'
#' Najib is curious about the relationship of genes in sets, the sets, and the
#' genes that comprise those sets.  This is pushing gsva towards a oroborous-ish
#' state.
#'
#' @param gsva_result Result from simple_gsva().
#' @param lst List of genes of interest.
#' @param freq_cutoff Minimum number of observations to be counted.
#' @param sig_weights When making venn diagrams, weight them?
#' @param gene_weights When venning genes, weight them?
#' @return List containing some venns, lists, and such.
#' @seealso [Vennerable] [simple_gsva()]
#' @export
intersect_signatures <- function(gsva_result, lst, freq_cutoff = 2,
                                 sig_weights = TRUE, gene_weights = TRUE) {
  sig_venn <- Vennerable::Venn(Sets = lst)
  Vennerable::plot(sig_venn, doWeights = sig_weights)
  sig_plot <- grDevices::recordPlot()
  sig_int <- sig_venn@IntersectionSets
  annot <- fData(gsva_result[["expt"]])
  sig_genes <- list()
  gene_venn_lst <- list()
  venn_names <- list()
  ## Skip the non-existant set of 00, thus 2:length()
  ## Top level loop iterates through the observed intersections/unions from Vennerable.
  for (i in seq(from = 2, to = length(sig_int))) {
    name <- names(sig_int)[i]
    ## Make a human readable version of the venn names.
    name_chars <- strsplit(x = name, split = "")[[1]]
    venn_name <- ""
    for (c in seq_along(name_chars)) {
      char <- name_chars[[c]]
      if (char == "1") {
        venn_name <- glue("{venn_name}_{names(lst)[c]}")
      }
    }
    venn_name <- gsub(pattern = "^_", replacement = "", x = venn_name)

    sigs <- sig_int[[i]]
    sig_annot <- annot[sigs, ]
    gene_ids <- sig_annot[["ids"]]
    internal_ret <- list()
    ## This loop iterates through the set of observed gene IDs in each intersection
    for (j in seq_along(gene_ids)) {
      id_lst <- gene_ids[j]
      ids <- strsplit(id_lst, ", ")[[1]]
      ## Finally, we count how many times each id is observed in each signature
      for (k in seq_along(ids)) {
        element <- ids[k]
        if (is.null(internal_ret[[element]])) {
          internal_ret[[element]] <- 1
        } else {
          internal_ret[[element]] <- internal_ret[[element]] + 1
        }
      }
    }
    internal_ret <- internal_ret[order(as.numeric(internal_ret), decreasing = TRUE)]
    ret <- as.numeric(internal_ret)
    names(ret) <- names(internal_ret)
    sig_genes[[venn_name]] <- ret
    gene_venn_lst[[venn_name]] <- names(ret[ret >= freq_cutoff])
  }

  gene_venn <- Vennerable::Venn(Sets = gene_venn_lst)
  Vennerable::plot(gene_venn, doWeights = gene_weights)
  gene_venn_plot <- grDevices::recordPlot()
  gene_int <- gene_venn@IntersectionSets

  retlst <- list(
      "signature_venn" = sig_venn,
      "signature_intersection" = sig_int,
      "signature_venn_plot" = sig_plot,
      "signature_genes" = sig_genes,
      "gene_venn" = gene_venn,
      "gene_intersection" = gene_int,
      "gene_venn_plot" = gene_venn_plot)
  return(retlst)
}

#' Score the results from simple_gsva().
#'
#' Yeah, this is a bit meta, but the scores from gsva seem a bit meaningless to
#' me, so I decided to look at the distribution of observed scores in some of my
#' data; I quickly realized that they follow a nicely normal distribution.
#' Therefore, I thought to calculate some scores of gsva() using that
#' information.
#'
#' The nicest thing in this, I think, is that it provides its scoring metric(s)
#' according to a few different possibilities, including:
#'   * the mean of samples found in an experimental factor
#'   * All provided scores against the distribution of observed scores as
#'     z-scores.
#'   * A single score against all scores.
#'   * Rows (gene sets) against the set of all gene sets.
#'
#' @param gsva_result Input result from simple_gsva()
#' @param score What type of scoring to perform, against a value, column, row?
#' @param category What category to use as baseline?
#' @param factor Which experimental factor to compare against?
#' @param sample Which sample to compare against?
#' @param factor_column When comparing against an experimental factor, which design
#'  column to use to find it?
#' @param method mean or median when when bringing together values?
#' @param label_size By default, enlarge the labels to readable at the cost of losing some.
#' @param col_margin Attempt to make heatmaps fit better on the screen with this and...
#' @param row_margin this parameter
#' @param cutoff Highlight only the categories deemed more significant than this.
#' @return The scores according to the provided category, factor, sample, or
#'  score(s).
#' @seealso [simple_gsva()]
#' @export
score_gsva_likelihoods <- function(gsva_result, score = NULL, category = NULL,
                                   factor = NULL, sample = NULL, factor_column = "condition",
                                   method = "mean", label_size = NULL,
                                   col_margin = 6, row_margin = 12, cutoff = 0.95) {
  values <- exprs(gsva_result[["expt"]])
  design <- pData(gsva_result[["expt"]])
  gsva_pca <- plot_pca(gsva_result[["expt"]])

  ## Start off with a plot of the gsva return values.
  color_range <- c("#00007F", "blue", "#007FFF", "cyan",
                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
  jet_colors <- grDevices::colorRampPalette(color_range)
  starting_ht <- NULL
  tmp_file <- tmpmd5file(pattern = "heat", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  if (is.null(label_size)) {
    starting_ht <- heatmap.3(values, trace = "none", col = jet_colors,
                             margins = c(col_margin, row_margin))
  } else {
    starting_ht <- heatmap.3(values, trace = "none", col = jet_colors,
                             margins = c(col_margin, row_margin),
                             cexCol = label_size, cexRow = label_size)
  }
  starting_ht_plot <- grDevices::recordPlot()
  dev.off()
  removed <- suppressWarnings(file.remove(tmp_file))
  removed <- unlink(dirname(tmp_file))
  tests <- test_values <- against_values <- NULL
  choice <- NULL
  if (is.null(score) & is.null(category) & is.null(sample) & is.null(factor)) {
    message("Nothing was requested, examining the first column scores: ",
            colnames(values)[1], ".")
    sample <- 1
    tests <- as.numeric(values[, sample])
    choice <- "column"
  } else if (!is.null(factor_column)) {
    choice <- "against"
  } else if (!is.null(score)) {
    message("Examining the score: ", score, " against the data.")
    tests <- score
    choice <- "value"
  } else if (!is.null(category)) {
    message("Examining row: ", category, " against the data.")
    tests <- as.numeric(values[category, ])
    choice <- "row"
  } else {
    tests <- as.numeric(values[, sample])
    choice <- "column"
  }

  population_values <- values
  ## ok, this function will give 'better' scores for higher
  ## gsva scores, so that is good, I will therefore just say 'higher is better'.
  cheesy_likelihood <- function(test) {
    gsva_mean <- mean(population_values)
    gsva_sd <- sd(population_values)
    num_values <- length(population_values)
    pop_sd <- gsva_sd * sqrt((num_values - 1) / num_values)
    notz <- (test - gsva_mean) / pop_sd
    likelihood <- pnorm(notz)
    return(likelihood)
  }

  results <- test_results <- against_results <- t_vs_a_results <- score_plot <- NULL
  if (choice == "against") {
    message("Testing each factor against the others.")
    fact_lvls <- levels(as.factor(design[[factor_column]]))
    result_df <- data.frame()
    for (f in seq_along(fact_lvls)) {
      fact <- fact_lvls[f]
      message("Scoring ", fact, " against everything else.")
      sample_idx <- design[[factor_column]] == fact
      if (sum(sample_idx) == 0) {
        ## Just in case the factor has empty levels.
        next
      }

      test_values <- values[, sample_idx]
      against_values <- values[, !sample_idx]
      population_values <- against_values
      if (method == "mean") {
        if (sum(sample_idx) == 1) {
          tests <- test_values
        } else {
          tests <- rowMeans(test_values)
        }
        againsts <- rowMeans(against_values)
      } else {
        if (sum(sample_idx) == 1) {
          tests <- test_values
        } else {
          tests <- Biobase::rowMedians(test_values)
        }
        againsts <- Biobase::rowMedians(against_values)
      }
      a_column <- sapply(X = tests, FUN = cheesy_likelihood)
      if (f == 1) {
        result_df <- as.data.frame(a_column)
      } else {
        result_df <- cbind(result_df, a_column)
      }
    } ## End iterating over every level in the chosen factor.
    colnames(result_df) <- fact_lvls
    heat_colors <- grDevices::colorRampPalette(c("white", "black"))
    tmp_file <- tmpmd5file(pattern = "heat", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    ht_result <- heatmap.3(as.matrix(result_df), trace = "none", col = heat_colors,
                           margins = c(col_margin, row_margin), Colv = FALSE,
                           dendrogram = "row", cexCol = label_size, cexRow = label_size)
    score_plot <- grDevices::recordPlot()
    dev.off()
    removed <- suppressWarnings(file.remove(tmp_file))
    removed <- unlink(dirname(tmp_file))
    test_results <- result_df
  } else if (choice == "column") {
    test_results <- sapply(X = tests, FUN = cheesy_likelihood)
    names(test_results) <- rownames(values)
    score_plot <- plot_histogram(test_results)
  } else if (choice == "row") {
    test_results <- sapply(X = tests, FUN = cheesy_likelihood)
    names(test_results) <- colnames(values)
    score_plot <- plot_histogram(test_results)
  } else {
    names(test_results) <- "value"
    score_plot <- plot_histogram(test_results)
  }

  retlist <- list(
      "pca" = gsva_pca,
      "raw_plot" = starting_ht_plot,
      "likelihoods" = test_results,
      "likelihood_plot" = score_plot)
  return(retlist)
}

#' Provide some defaults and guidance when attempting to use gsva.
#'
#' gsva seems to hold a tremendous amount of potential.  Unfortunately, it is
#' somewhat opaque and its requirements are difficult to pin down.  This
#' function will hopefully provide some of the requisite defaults and do some
#' sanity checking to make it more likely that a gsva analysis will succeed.
#'
#' @param expt Expt object to be analyzed.
#' @param signatures Provide an alternate set of signatures (GeneSetCollections)
#' @param data_pkg What package contains the requisite dataset?
#' @param signature_category Specify a subset category to extract from the signatures database.
#' @param cores How many CPUs to use?
#' @param current_id Where did the IDs of the genes come from?
#' @param required_id gsva (I assume) always requires ENTREZ IDs, but just in
#'  case this is a parameter.
#' @param min_catsize Minimum category size to consider interesting (passed to gsva()).
#' @param orgdb What is the data source for the rownames()?
#' @param method Which gsva method to use? Changed this from gsva to ssgsea
#'  because it was throwing segmentation faults.
#' @param kcdf Options for the gsva methods.
#' @param ranking another gsva option.
#' @param msig_db File contining msigdb annotations.
#' @param wanted_meta Desired metadata elements from the mxig_xml file.
#' @param mx_diff Passed to gsva(), I do not remember what it does.
#' @param verbose Print some information while running?
#' @param id_type Specify the ID type when loading the signature database.
#' @return List containing three elements: first a modified expressionset using
#'  the result of gsva in place of the original expression data; second the
#'  result from gsva, and third a data frame of the annotation data for the
#'  gene sets in the expressionset.  This seems a bit redundant, perhaps I
#'  should revisit it?
#' @seealso [GSEABase] [load_gmt_signatures()] [create_expt()] [GSVA]
#' @export
simple_gsva <- function(expt, signatures = "c2BroadSets", data_pkg = "GSVAdata",
                        signature_category = "c2", cores = NULL, current_id = "ENSEMBL",
                        required_id = "ENTREZID", min_catsize = 5, orgdb = "org.Hs.eg.db",
                        method = "ssgsea", kcdf = NULL, ranking = FALSE, msig_db = NULL,
                        ## wanted_meta = c("ORGANISM", "DESCRIPTION_BRIEF", "AUTHORS", "PMID"),
                        wanted_meta = "all", mx_diff = TRUE, verbose = FALSE, id_type = "entrez") {
  if (is.null(kcdf)) {
    if (expt[["state"]][["transform"]] == "raw") {
      kcdf <- "Poisson"
    } else {
      kcdf <- "Gaussian"
    }
  }

  if (is.null(cores)) {
    cores <- min(parallel::detectCores() - 1, 16)
  }

  if (!is.null(msig_db)) {
    if (!file.exists(msig_db)) {
      stop("The msig_db parameter was defined, but the file does not exist.")
    }
  }

  ## Make sure some data is loaded.  I will no longer assume anything here.
  ## Here is how I will decide:
  ## 1.  If signatures is a (string)filename ending in '.gmt', then extract the genesetlists and use it.
  ## 2.  If signatures is a string, then load the data_pkg, presumably GSVAdata.
  ## 3.  If signatures is not a string, assume it is a genesetlist/geneset and use that.

  ## Assume the desired category is c2 unless specified.
  if (is.null(signature_category)) {
    signature_category <- "c2"
  }
  signature_data <- load_gmt_signatures(signatures = signatures, data_pkg = data_pkg,
                                            signature_category = signature_category,
                                            id_type = id_type)

  ## The expressionset must have the annotation field filled in for gsva to
  ## work.
  eset_annotation <- annotation(expt)
  eset_pattern <- grepl(pattern = "Fill me in", x = annotation(expt))
  if (length(eset_annotation) == 0 | isTRUE(eset_pattern)) {
    message("gsva requires the annotation field to be filled in. Setting it to orgdb given.")
    annotation(expt) <- orgdb
  }

  eset <- expt
  ## The rownames() of the expressionset must be in ENTREZIDs for gsva to work.
  if (current_id != required_id | !is.integer(grep("ENSG", rownames(exprs(expt))))) {
    message("Converting the rownames() of the expressionset to ENTREZID.")
    ##tt <- sm(library(orgdb, character.only = TRUE))
    lib_result <- sm(requireNamespace(orgdb))
    att_result <- sm(try(attachNamespace(orgdb), silent = TRUE))
    old_ids <- rownames(exprs(expt))
    new_ids <- sm(AnnotationDbi::select(x = get0(orgdb),
                                        keys = old_ids,
                                        keytype = current_id,
                                        columns = c(required_id)))
    new_idx <- complete.cases(new_ids)
    if (!all(new_idx)) {
      message(sum(new_idx == FALSE),
              " ENSEMBL ID's didn't have a matching ENTEREZ ID. Dropping them now.")
    }
    duplicate_current <- duplicated(new_ids[[current_id]])
    new_ids <- new_ids[!duplicate_current, ]
    duplicate_required <- duplicated(new_ids[[required_id]])
    new_ids <- new_ids[!duplicate_required, ]
    rownames(new_ids) <- new_ids[[current_id]]
    dropme <- is.na(new_ids[[current_id]])
    new_ids <- new_ids[!dropme, ]
    dropme <- is.na(new_ids[[required_id]])
    new_ids <- new_ids[!dropme, ]
    message("Before conversion, the expressionset has ", length(rownames(exprs(expt))),
            " entries.")
    new_rownames <- new_ids[[required_id]]
    old_rownames <- new_ids[[current_id]]
    sub_expt <- expt[old_rownames, ]
    converted_expt <- set_expt_genenames(sub_expt, new_rownames)
    message("After conversion, the expressionset has ",
            length(rownames(exprs(converted_expt))),
            " entries.")
    rownames(exprs(converted_expt)) <- gsub(pattern = "^X", replacement = "",
                                     x = rownames(exprs(converted_expt)))
    eset <- converted_expt
    fData(eset)[[required_id]] <- rownames(fData(eset))
  }

  gsva_result <- NULL
  if (isTRUE(verbose)) {
    gsva_result <- GSVA::gsva(eset[["expressionset"]], signature_data, verbose = TRUE,
                              method = method, min.sz = min_catsize, kcdf = kcdf,
                              abs.ranking = ranking, parallel.sz = cores, mx.diff = mx_diff)
  } else {
    gsva_result <- sm(GSVA::gsva(eset[["expressionset"]], signature_data, verbose = TRUE,
                                 method = method, min.sz = min_catsize, kcdf = kcdf,
                                 abs.ranking = ranking, parallel.sz = cores, mx.diff = mx_diff))
  }
  fdata_df <- data.frame(row.names = rownames(exprs(gsva_result)))

  fdata_df[["description"]] <- ""
  fdata_df[["ids"]] <- ""
  message("Adding descriptions and IDs to the gene set annotations.")
  shared_sig_data <- signature_data[rownames(fdata_df)]
  for (i in seq_along(names(shared_sig_data))) {
    fdata_df[i, "description"] <- description(shared_sig_data[[i]])
    fdata_df[i, "ids"] <- toString(GSEABase::geneIds(shared_sig_data[[i]]))
  }

  if (!is.null(msig_db)) {
    message("Adding annotations from ", msig_db, ".")
    improved <- get_msigdb_metadata(msig_db = msig_db, wanted_meta = wanted_meta)
    ## Add improved to fData
    fdata_df <- merge(fdata_df, improved, by = "row.names", all.x = TRUE)
    rownames(fdata_df) <- fdata_df[["Row.names"]]
    fdata_df[["Row.names"]] <- NULL
  }

  fData(gsva_result) <- fdata_df
  new_expt <- expt
  new_expt[["expressionset"]] <- gsva_result

  retlist <- list(
      "method" = method,
      "signatures" = signatures,
      "signature_category" = signature_category,
      "required_id" = required_id,
      "min_catsize" = min_catsize,
      "expt" = new_expt,
      "gsva" = gsva_result,
      "fdata" = fdata_df)
  class(retlist) <- "gsva_result"
  return(retlist)
}

#' Invoke xCell and pretty-ify the result.
#'
#' I initially thought xCell might prove the best tool/method for exploring cell
#' deconvolution.  I slowly figured out its limitations, but still think it
#' seems pretty nifty for its use case.  Thus this function is intended to make
#' invoking it easier/faster.
#'
#' @param expt Expressionset to query.
#' @param signatures Alternate set of signatures to use.
#' @param genes Subset of genes to query.
#' @param spill The xCell spill parameter.
#' @param expected_types Set of assumed types in the data.
#' @param label_size How large to make labels when printing the final heatmap.
#' @param col_margin Used by par() when printing the final heatmap.
#' @param row_margin Ibid.
#' @param sig_cutoff Only keep celltypes with a significance better than this.
#' @param verbose Print some extra information during runtime.
#' @param cores How many CPUs to use?
#' @param ... Extra arguments when normalizing the data for use with xCell.
#' @return Small list providing the output from xCell, the set of signatures,
#'  and heatmap.
#' @seealso [xCell]
#' @export
simple_xcell <- function(expt, signatures = NULL, genes = NULL, spill = NULL,
                         expected_types = NULL, label_size = NULL, col_margin = 6,
                         row_margin = 12, sig_cutoff = 0.2, verbose = TRUE, cores = 4, ...) {
  arglist <- list(...)
  xcell_annot <- load_biomart_annotations()
  xref <- xcell_annot[["annotation"]][, c("ensembl_gene_id", "hgnc_symbol")]
  expt_state <- expt[["state"]][["conversion"]]
  xcell_eset <- NULL
  if (expt_state != "rpkm") {
    message("xCell strongly perfers rpkm values, re-normalizing now.")
    xcell_eset <- normalize_expt(expt, convert = "rpkm",
                                 ...)
  } else {
    xcell_eset <- normalize_expt(expt, norm = arglist[["norm"]], convert = arglist[["convert"]],
                                 filter = arglist[["filter"]], batch = arglist[["batch"]])
  }
  xcell_mtrx <- exprs(xcell_eset)
  xcell_na <- is.na(xcell_mtrx)
  xcell_mtrx[xcell_na] <- 0
  xcell_input <- merge(xcell_mtrx, xref, by.x = "row.names", by.y = "ensembl_gene_id")
  rownames(xcell_input) <- make.names(xcell_input[["hgnc_symbol"]], unique = TRUE)
  xcell_input[["Row.names"]] <- NULL
  xcell_input[["hgnc_symbol"]] <- NULL

  xCell.data <- NULL
  tt <- requireNamespace("xCell")
  data("xCell.data", package = "xCell")
  if (is.null(signatures)) {
    signatures <- xCell.data[["signatures"]]
  }
  if (is.null(genes)) {
    genes <- xCell.data[["genes"]]
  }
  if (is.null(spill)) {
    spill <- xCell.data[["spill"]]
  }

  xcell_result <- NULL
  if (isTRUE(verbose)) {
    xcell_result <- xCell::xCellAnalysis(expr = xcell_input, signatures = signatures,
                                         genes = genes, spill = spill,
                                         cell.types = expected_types, parallel.sz = cores)
  } else {
    xcell_result <- sm(xCell::xCellAnalysis(expr = xcell_input, signatures = signatures,
                                            genes = genes, spill = spill, parallel.sz = cores,
                                            cell.types = expected_types))
  }

  jet_colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  tmp_file <- tmpmd5file(pattern = "heat", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  if (is.null(label_size)) {
    ht <- heatmap.3(xcell_result, trace = "none", col = jet_colors,
                    margins = c(col_margin, row_margin))
  } else {
    ht <- heatmap.3(xcell_result, trace = "none", col = jet_colors,
                    margins = c(col_margin, row_margin),
                    cexCol = label_size, cexRow = label_size)
  }
  ht_plot <- grDevices::recordPlot()
  dev.off()
  ## sometimes the file does not get created.
  removed <- suppressWarnings(file.remove(tmp_file))
  removed <- unlink(dirname(tmp_file))

  sig_idx <- Biobase::rowMax(xcell_result) >= sig_cutoff
  sig_plot <- NULL
  sig_result <- NULL
  tmp_file <- tmpmd5file(pattern = "heat", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  if (sum(sig_idx) > 1) {
    sig_result <- as.matrix(xcell_result[sig_idx, ])
    if (is.null(label_size)) {
      ht <- heatmap.3(sig_result, trace = "none", col = jet_colors,
                      margins = c(col_margin, row_margin))
    } else {
      ht <- heatmap.3(sig_result, trace = "none", col = jet_colors,
                      margins = c(col_margin, row_margin),
                      cexCol = label_size, cexRow = label_size)
    }
    sig_plot <- grDevices::recordPlot()
  }
  dev.off()
  removed <- suppressWarnings(file.remove(tmp_file))
  removed <- unlink(dirname(tmp_file))

  retlist <- list(
      "xcell_input" = xcell_input,
      "xcell_result" = xcell_result,
      "signatures" = xCell.data[["signatures"]],
      "heatmap" = ht_plot,
      "sig_result" = sig_result,
      "sig_plot" = sig_plot)
  return(retlist)
}

#' Write out my various attempts at making sense of gsva.
#'
#' While I am trying to make sense of gsva, I will use this function to write
#' out the results I get so I can pass them to Najib/Maria Adelaida/Theresa to
#' see if I am making sense.
#'
#' @param retlist Result from running get_sig_gsva
#' @param excel Excel file to write
#' @param plot_dim Plot dimensions, likely needs adjustment.
#' @seealso [simple_gsva()] [score_gsva_likelihoods()] [get_sig_gsva_categories()]
#' @export
write_gsva <- function(retlist, excel, plot_dim = 6) {
  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]

  methods <- list(
      ## Using the correct character results in a warning from R CMD check... what to do?
      ## "gsva" = "Hänzelmann et al, 2013",
      "gsva" = "Hanzelmann et al, 2013",
      "ssgsea" = "Barbie et al, 2009",
      "zscore" = "Lee et al, 2008",
      "plage" = "Tomfohr et al, 2005")

  db_used <- retlist[["input"]][["signatures"]]
  if (class(db_used)[1] != "character") {
    db_used <- "user provided"
  }

  method <- retlist[["input"]][["method"]]
  if (method %in% names(methods)) {
    method <- paste0(method, ": ", methods[method])
  }

  ## Write the legend.
  legend <- data.frame(rbind(
      c("Signature database used:", db_used),
      c("Database subset used:", retlist[["input"]][["signature_category"]]),
      c("Required ID type:", retlist[["input"]][["required_id"]]),
      c("Minimum category size:", retlist[["input"]][["min_catsize"]]),
      c("GSVA method used:", method),
      c("", ""),
      c("Sheet 1: gsva_scores", "All scores as provided by gsva()."),
      c("Sheet 2: score_gsva_likelihoods", "All likelihood scores calculated using pnorm() of the values."),
      c("Sheet 3: factor_likelihoods", "Likelihood values for each experimental factor."),
      c("Sheet 4: subset", "GSVA scores for the categories deemed 'significant' using sheet 2/3."),
      c("Sheet 5 on:", "Limma scoring of differential signatures.")),
      stringsAsFactors = FALSE)
  colnames(legend) <- c("Term", "Definition")
  xls_result <- write_xlsx(
      wb, data = legend, sheet = "legend", rownames = FALSE,
      title = "Summary and sheets in this workbook.")
  xl_result <- openxlsx::writeData(wb = wb, sheet = "legend", x = "PCA of categories vs sample type.",
                                   startRow = 1, startCol = 8)
  try_result <- xlsx_insert_png(retlist[["score_pca"]], wb = wb, sheet = "legend",
                              start_row = 2, start_col = 8,
                              width=(plot_dim * 3/2), height = plot_dim,
                              plotname = "gsva_pca", savedir = excel_basename)

  ## Write the result from gsva()
  xls_result <- write_xlsx(data = retlist[["gsva_table"]], wb = wb, sheet = "gsva_scores")
  current_column <- xls_result[["end_col"]] + 2
  current_row <- 1
  plot_width <- 8
  plot_height <- 16
  try_result <- xlsx_insert_png(a_plot = retlist[["raw_plot"]], wb = wb, sheet = "gsva_scores",
                              start_row = current_row, start_col = current_column,
                              width = plot_width, height = plot_height)

  ## Write the likelihoods
  xls_result <- write_xlsx(data = retlist[["likelihood_table"]], wb = wb, sheet = "likelihood_scores")
  current_column <- xls_result[["end_col"]] + 2
  current_row <- 1
  plot_width <- 6
  plot_height <- 15
  try_result <- xlsx_insert_png(a_plot = retlist[["score_plot"]], wb = wb, sheet = "likelihood_scores",
                              start_row = current_row, start_col = current_column,
                              width = plot_width, height = plot_height)

  ## Write the subset
  xls_result <- write_xlsx(data = retlist[["subset_table"]], wb = wb, sheet = "subset_table")
  current_column <- xls_result[["end_col"]] + 2
  current_row <- 1
  plot_width <- 6
  plot_height <- 6
  try_result <- xlsx_insert_png(a_plot = retlist[["subset_plot"]], wb = wb, sheet = "subset_table",
                              start_row = current_row, start_col = current_column,
                              width = plot_width, height = plot_height)

  limma_tables <- retlist[["gsva_limma"]][["all_tables"]]
  table_names <- names(limma_tables)
  for (i in seq_along(table_names)) {
    table_name <- table_names[i]
    title <- glue::glue("Result from using limma to compare {table_name}.")
    table <- limma_tables[[table_name]]
    table_idx <- order(table[["adj.P.Val"]], decreasing = FALSE)
    table <- table[table_idx, ]
    xls_result <- write_xlsx(data = table, wb = wb, sheet = table_name, title = title)
  }

  save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  if (class(save_result)[1] == "try-error") {
    message("Saving xlsx failed.")
  }
  return(save_result)
}

## EOF
