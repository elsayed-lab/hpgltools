#' Add binary state information to the scd.
#'
#' I am adding these only so that it is easier to visualize the cells
#' compared when performing FindAllMarkers(); e.g. it compares each
#' identity to all others; so I imagine it would be nice to see a
#' dimplot or something of each state vs. all others as a binary pair
#' rather than as n separate groups.
add_binary_states <- function(scd, column = NULL) {
  identity_levels <- levels(as.factor(Seurat::Idents(object = scd)))
  if (!is.null(column)) {
    identity_levels <- levels(as.factor(scd[[column]]))
  }
  cell_ids <- colnames(scd)

  for (group in identity_levels) {
    not_value <- paste0("Not ", group)
    if (is.null(column)) {
      ## A little lisp predicate happiness
      group_name <- paste0(group, "_identp")
      scd@meta.data[[group_name]] <- not_value
      true_idx <- Seurat::Idents(object = scd) == group
    } else {
      group_name <- paste0(group, "_", column, "p")
      scd@meta.data[[group_name]] <- not_value
      true_idx <- scd[[column]] == group
    }
    mesg("There are ", sum(true_idx), " cells in group ", group, ".")
    scd@meta.data[true_idx, group_name] <- group
  }
  return(scd)
}

#' Add VDJ information using some code I found.
#'
#' The original implementation of this idea resides at:
#' https://ucdavis-bioinformatics-training.github.io/2020-Advanced_Single_Cell_RNA_Seq/data_analysis/VDJ_Analysis_fixed
#'
#' @param obj Seurat object to which we will add some information.
#' @param start_path root of the 10x data in which the vdj information should reside.
#' @param type The type of VDJ we expect, heavy(B) or light(T).
#' @return The Seurat object with some new information.
#' @export
add_clonotype_annotations <- function(scd, start_path, type = "t") {
  vdj_dir <- file.path(start_path, paste0("vdj_", type))
  vdj_csv <- file.path(vdj_dir, "filtered_contig_annotations.csv")
  reference_csv <- file.path(vdj_dir, "clonotypes.csv")
  tcr <- readr::read_csv(vdj_csv, show_col_types = FALSE)
  ref <- readr::read_csv(reference_csv, show_col_types = FALSE)

  tcr_duplicate_barcode_idx <- duplicated(tcr[["barcode"]])
  tcr_nodup <- tcr[!tcr_duplicate_barcode_idx, ]
  mesg("  After deduplication, the vdj data has ", nrow(tcr_nodup), " rows.")

  both <- as.data.frame(merge(tcr_nodup, ref, by.x="raw_clonotype_id", by.y="clonotype_id"))
  rownames(both) <- both[["barcode"]]
  both[["barcode"]] <- NULL

  scd <- Seurat::AddMetaData(object = scd, metadata = both)
  return(scd)
}

#' Add a df of clonotype observations by cell cluster to the @misc slot of the Seurat data.
count_clonotype_by_cluster <- function(scd, column = "res0p2_clusters",
                                       clono_column = "raw_clonotype_id", add_sum = TRUE) {
  meta <- scd@meta.data
  retstring <- glue::glue("This result may be found in the {column} element of scd@misc.")
  message(retstring)
  cluster_levels <- levels(as.factor(meta[[column]]))
  ret <- data.frame(rownames = cluster_levels)
  ret[["cluster_name"]] <- ""
  ret[["cells"]] <- 0
  ret[["has_clono"]] <- 0
  for (cluster_num in seq_along(cluster_levels)) {
    cluster <- cluster_levels[cluster_num]
    ret[cluster_num, "cluster_name"] <- cluster
    ret[cluster_num, "cells"] <- sum(meta[[column]] == cluster)
    ret[cluster_num, "has_clono"] <- sum(meta[[column]] == cluster &
                                         !is.na(meta[["raw_clonotype_id"]]))
  }
  ret[["rownames"]] <- NULL
  ret[["proportion"]] <- ret[["has_clono"]] / ret[["cells"]]

  if (isTRUE(add_sum)) {
    sum_cells <- sum(ret[["cells"]])
    sum_clono <- sum(ret[["has_clono"]])
    sum_vector = c("cluster_sum", sum_cells, sum_clono, sum_clono / sum_cells)
    ret <- rbind(ret, sum_vector)
  }

  scd@misc[[column]] <- ret
  return(scd)
}

#' Create a combined seurat object from a sample sheet.
#'
#' I would like to have a simpler function for creating seurat data
#' structures similar to my create_expt().  This will try to do so.
#'
#' @param metadata Sample sheet.
#' @param expression_column Metadata column containing the base
#'  location of the cellranger outputs.
#' @param vdj_t_column Column, which if filled in, tells this to look
#'  for VDJ information specific to light chains.
#' @param vdj_b_column Column, which if filled in, tells this to look
#'  for VDJ information specific to heavy chains.
#' @param prefix Arbitrary prefix for the location information,
#'  included because I am messing with cellranger and have multiple
#'  output directories and want to be able to switch between them.
#' @param separate When true, this function should return a list
#'  comprised of the individual sample objects.
#' @return Either a list or merged seurat object(s).
#' @export
create_scd <- function(metadata, expression_column = "gexfile",
                       vdj_t_column = "vdjtcells", vdj_b_column = "vdjbcells",
                       prefix = NULL, separate = FALSE, types = "gex",
                       mito_pattern = "^mt-", ribo_pattern = "^Rp[sl]") {
  if ("character" %in% class(metadata)) {
    metadata <- extract_metadata(metadata)
  }

  count <- 1
  obj <- NULL
  merged <- NULL
  if (isTRUE(separate)) {
    merged <- list()
  }

  for (dirname in metadata[[expression_column]]) {
    path_name <- basename(dirname)
    start_path <- file.path(dirname, "outs", "per_sample_outs", path_name)
    mesg("Loading sample ", path_name, " from: ", start_path, ".")
    if (!is.null(prefix)) {
      start_path <- file.path(prefix, dirname, "outs", "per_sample_outs", path_name)
    }
    full_path <- file.path(start_path, "count", "sample_filtered_feature_bc_matrix")
    if (length(types) == 1) {
        obj <- Seurat::Read10X(full_path) %>%
            Seurat::CreateSeuratObject(project = path_name)
    } else {
        ## If there are more than 1 type, then obj should be a list.
        obj <- Seurat::Read10X(full_path)
        gex <- Seurat::CreateSeuratObject(obj[[1]])
        Seurat::Idents(object = gex) <- path_name
        ## gex@meta.data[["orig.ident"]] <- path_name
        ab <- Seurat::CreateAssayObject(obj[[2]])
        gex[["ab"]] <- ab
        obj <- gex
    }

    if (!is.null(metadata[[vdj_t_column]])) {
      mesg("  Adding light-chain clonotype information for sample: ",
           path_name, ".")
      obj <- add_clonotype_annotations(obj, type = "t",
                                       start_path = start_path)
    }
    if (!is.null(metadata[[vdj_b_column]])) {
      mesg("  Adding heavy-chain clonotype information for sample: ",
           path_name, ".")
      obj <- add_clonotype_annotations(obj, type = "b",
                                       start_path = start_path)
    }

    if (isTRUE(separate)) {
      merged[[path_name]] <- obj
    } else {
      if (count == 1) {
        merged <- obj
      } else {
        merged <- merge(merged, obj)
      }
    }
    count <- count + 1
  } ## End of loop

  merged@misc[["sample_metadata"]] <- metadata

  merged[["pct_mito"]] <- Seurat::PercentageFeatureSet(merged, pattern = mito_pattern)
  merged[["pct_ribo"]] <- Seurat::PercentageFeatureSet(merged, pattern = ribo_pattern)

  return(merged)
}

#' Perform a series of filters on a single-cell dataset.
#'
#' This function should perform a series of relatively consistent
#' filters on a single-cell dataset, with options to play with the
#' various filters and their parameters.
#'
#' @param scd Single Cell Dataset to filter.
#' @param min_num_rna Drop cells with fewer than this number of
#'  observed RNA species.
#' @param min_pct_ribo Drop cells with less than this percentage of
#'  ribosomal protein RNAs observed.
#' @param max_pct_ribo Drop cells with more than this percentage of
#'  ribosomal protein RNAs observed.
#' @param max_pct_mito Drop cells with more than this percentage of
#'  mitochondrial RNA observed.
#' @param min_pct_mito Drop cells with less than this percentage of
#'  mitochondrial RNA observed.
#' @param mito_pattern Regex pattern to search RNA symbols for
#'  mitochondrial species.
#' @param ribo_pattern Regex pattern to search RNA symbols for
#'  ribosomal protein species.
#' @param min_gene_counts Drop genes across cells which are observed
#'  less than this number of times, I don't expect many of these.
#' @return Filtered scd
#' @export
filter_scd <- function(scd, min_num_rna = 200, max_num_rna = NULL,
                       min_pct_ribo = 5, max_pct_ribo = NULL, remerge = NULL,
                       max_pct_mito = 15, min_pct_mito = NULL, mito_pattern = "^mt-",
                       ribo_pattern = "^Rp[sl]", min_gene_counts = 3, verbose = FALSE) {
  current_cells <- ncol(scd)
  current_genes <- nrow(scd)
  current_counts <- sum(BiocGenerics::colSums(scd))
  start_cells <- current_cells
  start_genes <- current_genes
  start_counts <- current_counts

  remerge_assay <- NULL
  if (!is.null(remerge)) {
    remerge_assay <- scd[[remerge]]
  }

  mesg("We are starting with ", current_cells, " cells, ", current_genes,
       " genes, and ", current_counts, " counts.")

  if (is.null(scd@meta.data[["pct_mito"]])) {
    scd[["pct_mito"]] <- Seurat::PercentageFeatureSet(scd, pattern = mito_pattern)
  }
  if (is.null(scd@meta.data[["pct_ribo"]])) {
    scd[["pct_ribo"]] <- Seurat::PercentageFeatureSet(scd, pattern = ribo_pattern)
  }
  filt_scd <- scd

  pre_plots <- plot_seurat_scatter(scd)
  sufficient_rna <- nFeature_RNA <- NULL
  if (!is.null(min_num_rna) & !is.null(max_num_rna)) {
    mesg("Filtering both less than minimum and more than maximum.")
    sufficient_rna <- Seurat::WhichCells(filt_scd, expression = nFeature_RNA >= min_num_rna) &
      Seurat::WhichCells(filt_scd, expression = nFeature_RNA <= max_num_rna)
  } else if (!is.null(min_num_rna)) {
    mesg("Filtering less than minimum number of RNAs observed.")
    sufficient_rna <- Seurat::WhichCells(filt_scd, expression = nFeature_RNA >= min_num_rna)
  } else if (!is.null(max_num_rna)) {
    mesg("Filtering more than maximum number of RNAs observed.")
    sufficient_rna <- Seurat::WhichCells(filt_scd, expression = nFeature_RNA >= max_num_rna)
  }

  if (is.null(sufficient_rna)) {
    mesg("Not filtering based on sufficient number of genes observed.")
  } else if (length(sufficient_rna) == current_cells) {
    mesg("The provided filters on sufficient genes observed did not remove any cells.")
  } else {
    filt_scd <- subset(filt_scd, cells = sufficient_rna)
    new_cells <- ncol(filt_scd)
    new_genes <- nrow(filt_scd)
    new_counts <- sum(BiocGenerics::colSums(filt_scd))
    delta_cells <- current_cells - new_cells
    delta_genes <- current_genes - new_genes
    delta_counts <- current_counts - new_counts
    mesg("The sufficient genes observed filter removed: ", delta_cells, " cells and ",
         delta_counts, " reads.")
    current_cells <- new_cells
    current_genes <- new_genes
    current_counts <- new_counts
  }

  ##mesg("We now have ", current_cells, " cells, ", current_genes,
  ##     " genes, and ", current_counts, " counts.")
  if (is.null(min_gene_counts)) {
    mesg("Not filtering based on the minimum number of counts across cells observed.")
  } else {
    sufficiently_observed <- BiocGenerics::rowSums(filt_scd) >= min_gene_counts
    sufficient_genes <- rownames(filt_scd)[sufficiently_observed]
    if (current_genes == length(sufficient_genes)) {
      mesg("The minimum number of counts/gene across cells filter did not remove any genes.")
    } else {
      filt_scd <- subset(filt_scd, features = sufficient_genes)
      new_cells <- ncol(filt_scd)
      new_genes <- nrow(filt_scd)
      new_counts <- sum(BiocGenerics::colSums(filt_scd))
      delta_cells <- current_cells - new_cells
      delta_genes <- current_genes - new_genes
      delta_counts <- current_counts - new_counts
      mesg("The sufficient counts/gene across cells filter removed: ", delta_cells, " cells, ",
           delta_genes, " genes, and ", delta_counts, " reads.")
      current_cells <- new_cells
      current_genes <- new_genes
      current_counts <- new_counts
    }
  }

  ##mesg("We now have ", current_cells, " cells, ", current_genes,
  ##     " genes, and ", current_counts, " counts.")
  sufficient_ribo <- pct_ribo <- pct_mito <- NULL
  if (!is.null(min_pct_ribo) & !is.null(max_pct_ribo)) {
    sufficient_ribo <- Seurat::WhichCells(filt_scd, expression = pct_ribo >= min_pct_ribo) &
      Seurat::WhichCells(filt_scd, expression = pct_ribo <= max_pct_ribo)
  } else if (!is.null(min_pct_ribo)) {
    sufficient_ribo <- Seurat::WhichCells(filt_scd, expression = pct_ribo >= min_pct_ribo)
  } else if (!is.null(max_pct_ribo)) {
    sufficient_ribo <- Seurat::WhichCells(filt_scd, expression = pct_ribo <= max_pct_ribo)
  }
  if (is.null(sufficient_ribo)) {
    mesg("Not filtering based on sufficient ribosomal protein percentage observed.")
  } else if (length(sufficient_ribo) == current_cells) {
    mesg("The sufficient ribosomal percentage filter did not remove any cells.")
  } else {
    filt_scd <- subset(filt_scd, cells = sufficient_ribo)
    new_cells <- ncol(filt_scd)
    new_genes <- nrow(filt_scd)
    new_counts <- sum(BiocGenerics::colSums(filt_scd))
    delta_cells <- current_cells - new_cells
    delta_genes <- current_genes - new_genes
    delta_counts <- current_counts - new_counts
    mesg("The ribosomal protein filter removed: ", delta_cells, " cells and ",
         delta_counts, " reads.")
    current_cells <- new_cells
    current_genes <- new_genes
    current_counts <- new_counts
  }

  ##mesg("We now have ", current_cells, " cells, ", current_genes,
  ##     " genes, and ", current_counts, " counts.")
  sufficient_mito <- NULL
  if (!is.null(min_pct_mito) & !is.null(max_pct_mito)) {
    sufficient_mito <- Seurat::WhichCells(filt_scd, expression = pct_mito >= min_pct_mito) &
      Seurat::WhichCells(filt_scd, expression = pct_mito <= max_pct_mito)
  } else if (!is.null(min_pct_mito)) {
    sufficient_mito <- Seurat::WhichCells(filt_scd, expression = pct_mito >= min_pct_mito)
  } else if (!is.null(max_pct_mito)) {
    sufficient_mito <- Seurat::WhichCells(filt_scd, expression = pct_mito <= max_pct_mito)
  }
  if (is.null(sufficient_mito)) {
    mesg("Not filtering based on excessive mitochondrial RNA percentage observed.")
  } else if (length(sufficient_mito) == current_cells) {
    mesg("The mitochondrial filter did not remove any cells.")
  } else {
    filt_scd <- subset(filt_scd, cells = sufficient_mito)
    new_cells <- ncol(filt_scd)
    new_genes <- nrow(filt_scd)
    new_counts <- sum(BiocGenerics::colSums(filt_scd))
    delta_cells <- current_cells - new_cells
    delta_genes <- current_genes - new_genes
    delta_counts <- current_counts - new_counts
    mesg("The mitochondrial filter removed: ", delta_cells, " cells and ",
         delta_counts, " reads.")
    current_cells <- new_cells
    current_genes <- new_genes
    current_counts <- new_counts
  }

  filt_scd <- record_seurat_samples(filt_scd, type = "num_cells",
                                    column_name = "filt_num_cells", verbose = verbose) %>%
    record_seurat_samples(type = "nFeature_RNA", column_name = "filt_nfeature",
                          verbose = verbose) %>%
    record_seurat_samples(type = "nCount_RNA", column_name = "filt_ncount",
                          verbose = verbose) %>%
    record_seurat_samples(type = "pct_mito", column_name = "filt_pct_mito",
                          pattern = "^mt-", verbose = verbose) %>%
    record_seurat_samples(type = "pct_ribo", column_name = "filt_pct_ribo",
                          pattern = "^Rp[sl]", verbose = verbose)

  post_plots <- plot_seurat_scatter(filt_scd)

  message(sprintf("All filters removed %d (%f%%) cells, %d (%f%%) genes, %d (%f%%) counts.",
          start_cells - current_cells, (1.0 -(current_cells / start_cells)) * 100.0,
          start_genes - current_genes, (1.0 - (current_genes / start_genes)) * 100.0,
          start_counts - current_counts, (1.0 - (current_counts / start_counts)) * 100.0))

  if (!is.null(remerge)) {
    remerge_idx <- colnames(remerge_assay) %in% colnames(filt_scd)
    dropped <- sum(!remerge_idx)
    message("Filtering removed ", dropped, " cells from the ", remerge, " assay.")
    kept <- remerge_assay[, remerge_idx]
    wanted_cells <- colnames(remerge_assay)[remerge_idx]
    kept_assay <- subset(remerge_assay, cells = wanted_cells)
    filt_scd[[remerge]] <- kept_assay
  }

  filt_scd@misc[["pre_plots"]] <- pre_plots
  filt_scd@misc[["post_plots"]] <- post_plots
  return(filt_scd)
}

#' Add into the miscellaneous SCD slot a dataframe with some summary stats.
#'
#' There are some simple summaries which are nice to have on hand
#' regarding the number of RNAs, cells, rProteins, rmito observed.
#' This function collects them and drops them into a dataframe within
#' the slot 'misc' of the SCD.  I may also print to screen some
#' pretty skims of the results.
#'
#' @param scd Single Cell Dataset to query.
#' @param type Type of column to add to the metadata df, named for the column in the
#'  Cell-annotation table to query.
#' @param pattern Pattern used for regex-based queries.
#' @param column_name Name for the new column.
#' @param column_prefix Prefix added to the new column.
#' @param verbose Print the summaries to screen?
#' @return Give back the SCD with some new information.
#' @export
record_seurat_samples <- function(scd, type = "num_cells", pattern = NULL,
                                  column_name = NULL, column_prefix = NULL,
                                  verbose = TRUE, group = "Idents",
                                  assay = "RNA") {
  scd_meta <- scd@meta.data
  if (is.null(scd_meta[[group]])) {
    scd[["Idents"]] <- Seurat::Idents(scd)
    group <- "Idents"
  }

  sample_names <- levels(as.factor(scd[[group]]))
  test_slot <- "sample_metadata"
  if (group != "Idents") {
    test_slot <- paste0(group, "_metadata")
  }

  sample_meta <- data.frame()
  ## Check that we have existing information to append
  if (is.null(scd@misc[[test_slot]])) {
    sample_meta <- data.frame(rownames = sample_name)
  } else {
    sample_meta <- scd@misc[[test_slot]]
  }

  ## If a regex pattern is provided, use PercentageFeatureSet to add
  ## the percentage counts with genes matching that pattern (like
  ## ribosomal proteins or mitochondrial RNA)
  if (!is.null(pattern)) {
    ## scd[[type]] <- Seurat::PercentageFeatureSet(scd, pattern = pattern, assay = assay)
    test <- Seurat::PercentageFeatureSet(scd, pattern = pattern, assay = assay)
    if (class(test)[1] == "data.frame") {
      scd@meta.data[[type]] <- test[[1]]
    } else {
      scd@meta.data[[type]] <- test
    }
    scd_meta <- scd@meta.data
  }

  ## Everything from here down are methods to extract the information
  ## of interest.
  if (type == "num_cells") {
    sample_meta[[type]] <- 0
    for (s in 1:nrow(sample_meta)) {
      sample_name <- rownames(sample_meta)[s]
      sample_meta[s, type] <- sum(scd_meta[[group]] == sample_name)
    }
  } else {
    sample_meta <- skim_seurat_metadata(
        sample_meta, scd_meta, meta_query = type,
        group_column = group, summary_query = "numeric.mean",
        column_name = column_name, column_prefix = column_prefix,
        verbose = verbose)
  }

  ## Now put the new information into the misc list of the scd.
  if (group == "orig.ident") {
    scd@misc[["sample_metadata"]] <- sample_meta
  } else {
    misc_name <- paste0(group, "_metadata")
    scd@misc[[misc_name]] <- sample_meta
  }
  return(scd)
}

#' Make a few of the likely scatterplots provided by FeatureScatter.
#'
#' It seems I have used the same couple of scatter plots more often than others.
#'
#' @param scd SCD to plot.
#' @param set List of plots, use my favorites when NULL.
#' @return List of plots.
#' @export
plot_seurat_scatter <- function(scd, set = NULL) {
  plots <- list()
  if (is.null(set)) {
    set <- list(
        "ribo_vs_mito" = c("pct_ribo", "pct_mito"),
        "count_vs_genes" = c("nCount_RNA", "nFeature_RNA"),
        "count_vs_ribo" = c("nCount_RNA", "pct_ribo"),
        "count_vs_mito" = c("nCount_RNA", "pct_mito"))
  }

  for (p in 1:length(set)) {
    name <- names(set)[p]
    xy <- set[[p]]
    x_column <- xy[1]
    y_column <- xy[2]
    if (is.null(scd@meta.data[[x_column]])) {
      message("The x-axis data does not exist: ", x_column, ".")
      next
    }
    if (is.null(scd@meta.data[[y_column]])) {
      message("The y-axis data does not exist: ", y_column, ".")
      next
    }
    plots[[name]] <- Seurat::FeatureScatter(scd, x_column, y_column)
  }
  return(plots)
}

#' Summarize scores across observed clusters in a scd.
#'
#' Currently this assumes the set of outputs produced by Seurat's
#' AddModuleScore() for a gsc.  It summarizes those scores for each
#' cluster and gives back the mean, sd, and z.
#'
#' @param scd Input dataset.
#' @param fx Function to summarize, this may change.
#' @param column_prefix Prefix for the scores of interest.
#' @param column_range Explicitly set the range of interested columns.
#' @param cluster_column The column containing the information about
#'  cluster occupancy.
#' @param real_column_names The original columns get names like bob1
#'  to bobn, this can be used to make them more informative.
#' @param abbreviate When using mSigDB information, the category names
#'  are exceedingly long with often a consistent prefix.
#' @param min_mean Currently unused, but intended to filter out gsc
#'  which are not observed to any significant degree.
summarize_scd_clusters <- function(scd, fx = "mean", column_prefix = "descartes",
                                   column_range = NULL, cluster_column = "cluster_sample",
                                   real_column_names = NULL, abbreviate = TRUE,
                                   min_mean = NULL) {

  fx_df <- scd@meta.data
  if (is.null(fx_df[[cluster_column]])) {
    message("The column ", cluster_column, " is missing.")
    message("Here are the actual columns: ")
    print(colnames(fx_df))
    stop("Cannot find the cluster column.")
  }
  test_measure_column <- paste0(column_prefix, "1")
  if (is.null(fx_df[[test_measure_column]])) {
    stop("Unable to find the first measured column with prefix: ", column_prefix, ".")
  }
  column_pattern <- glue::glue("^{column_prefix}\\d+")
  wanted_columns <- grepl(x = colnames(fx_df), pattern = column_pattern)
  cluster_idx <- which(colnames(fx_df) == cluster_column)
  cluster_vector <- fx_df[[cluster_idx]]

  summary_df <- fx_df[, wanted_columns]
  if (is.null(real_column_names)) {
    real_column_names <- colnames(summary_df)
  } else {
    colnames(summary_df) <- real_column_names
  }
  summary_df <- cbind(cluster_vector, summary_df)
  summary_df[[1]] <- as.factor(summary_df[[1]])
  colnames(summary_df)[1] <- cluster_column

  mean_df <- summary_df %>%
    dplyr::group_by(!!sym(cluster_column)) %>%
    dplyr::summarise_at(vars(all_of(real_column_names)), list(name = mean))
  sd_df <- summary_df %>%
    dplyr::group_by(!!sym(cluster_column)) %>%
    dplyr::summarise_at(vars(all_of(real_column_names)), list(name = sd))

  mean_df <- as.data.frame(mean_df)
  rownames(mean_df) <- mean_df[[1]]
  mean_df[[1]] <- NULL
  colnames(mean_df) <- gsub(x = colnames(mean_df), pattern = "_name$",
                            replacement = "")
  sd_df <- as.data.frame(sd_df)
  rownames(sd_df) <- sd_df[[1]]
  sd_df[[1]] <- NULL
  colnames(sd_df) <- gsub(x = colnames(sd_df), pattern = "_name$",
                            replacement = "")
  z_df <- mean_df / sd_df
  z_df <- as.matrix(z_df)

  if (isTRUE(abbreviate)) {
    colnames(mean_df) <- abbreviate(colnames(mean_df), minlength = 10)
    rownames(mean_df) <- abbreviate(rownames(mean_df), minlength = 10)
    colnames(sd_df) <- abbreviate(colnames(sd_df), minlength = 10)
    rownames(sd_df) <- abbreviate(rownames(sd_df), minlength = 10)
    colnames(z_df) <- abbreviate(colnames(z_df), minlength = 10)
    rownames(z_df) <- abbreviate(rownames(z_df), minlength = 10)
  }

  mean_df <- as.matrix(mean_df)
  sd_df <- as.matrix(sd_df)
  z_df <- as.matrix(z_df)

  retlist <- list(
      "mean_df" = mean_df,
      "sd_df" = sd_df,
      "z_df" = z_df)
  return(retlist)
}

#' Use skimr to make pretty summaries of cell-data.
#'
#' I think I want to expand this to handle RNA summaries as well.
#'
#' @param sample_meta df of the known samples by name.
#' @param obj_meta The 'meta.data' slot of a SCD
#' @param meta_query Column to query.
#' @param group_column Column used to group the cells.
#' @param summary_query Which of the various data produced by skimr should be extracted?
#' @param column_name Add the new column with this name.
#' @param column_prefix And this prefix.
#' @param verbose Print the pretty skimr table?
#' @return df with some new meta(meta?)data.
#' @export
skim_seurat_metadata <- function(sample_meta, obj_meta, meta_query = "nCount_RNA",
                                 group_column = NULL, summary_query = "numeric.mean",
                                 column_name = NULL, column_prefix = NULL, verbose = TRUE) {
  if (is.null(column_name)) {
    column_name <- meta_query
  }
  if (is.null(column_prefix)) {
    if (grep(x = summary_query, pattern = "\\.")) {
      column_prefix <- gsub(x = summary_query, pattern = "^.*\\.(.*)?$", replacement = "\\1")
    }
  }

  if (!is.null(column_prefix)) {
    column_name <- paste0(column_prefix, "_", column_name)
  }

  if (isTRUE(verbose)) {
    sample_meta[[column_name]] <- obj_meta %>%
      group_by(!!rlang::sym(group_column)) %>%
      skimr::skim_tee(meta_query) %>%
      skimr::skim(meta_query) %>%
      dplyr::select(summary_query) %>%
      unlist()
  } else {
    sample_meta[[column_name]] <- obj_meta %>%
      group_by(!!rlang::sym(group_column)) %>%
      skim(meta_query) %>%
      dplyr::select(summary_query) %>%
      unlist()
  }
  return(sample_meta)
}
