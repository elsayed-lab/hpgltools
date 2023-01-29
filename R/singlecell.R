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
add_clonotype_annotations <- function(obj, start_path, type = "t") {
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

  obj <- Seurat::AddMetaData(object = obj, metadata = both)
  return(obj)
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
                               prefix = NULL, separate = FALSE) {
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
    obj <- Seurat::Read10X(full_path) %>%
      Seurat::CreateSeuratObject(project = path_name)

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
filter_seurat_data <- function(scd, min_num_rna = 200, min_pct_ribo = 5, max_pct_ribo = NULL,
                               max_pct_mito = 15, min_pct_mito = NULL, mito_pattern = "^mt-",
                               ribo_pattern = "^Rp[sl]", min_gene_counts = 3) {
  scd[["percent_mito"]] <- Seurat::PercentageFeatureSet(scd, pattern = mito_pattern)
  scd[["percent_ribo"]] <- Seurat::PercentageFeatureSet(scd, pattern = ribo_pattern)
  filt_scd <- scd

  pre_plots <- scatter_seurat(scd)
  sufficient_rna <- NULL
  if (!is.null(min_num_rna)) {
    sufficient_rna <- Seurat::WhichCells(scd, expression = nFeature_RNA >= min_num_rna)
  } else if (!is.null(max_num_rna)) {
    sufficient_rna <- Seurat::WhichCells(scd, expression = nFeature_RNA <= max_num_rna)
  }
  if (!is.null(sufficient_rna)) {
    filt_scd <- subset(scd, cells = sufficient_rna)
  }

  if (!is.null(min_gene_counts)) {
    sufficiently_observed <- rowSums(filt_scd) >= min_gene_counts
    sufficient_genes <- rownames(filt_scd)[sufficiently_observed]
    filt_scd <- subset(filt_scd, features = sufficient_genes)
  }

  sufficient_ribo <- NULL
  if (!is.null(min_pct_ribo)) {
    sufficient_ribo <- WhichCells(filt_scd, expression = percent_ribo >= min_pct_ribo)
  } else if (!is.null(max_pct_ribo)) {
    sufficient_ribo <- WhichCells(filt_scd, expression = percent_ribo <= max_pct_ribo)
  }
  if (!is.null(sufficient_ribo)) {
    filt_scd <- subset(scd, cells = sufficient_ribo)
  }

  sufficient_mito <- NULL
  if (!is.null(min_pct_mito)) {
    sufficient_mito <- WhichCells(filt_scd, expression = percent_mito >= min_pct_mito)
  } else if (!is.null(max_pct_mito)) {
    sufficient_mito <- WhichCells(filt_scd, expression = percent_mito <= max_pct_mito)
  }
  if (!is.null(sufficient_mito)) {
    filt_scd <- subset(scd, cells = sufficient_mito)
  }

  filt_scd <- record_seurat_samples(filt_scd, type = "num_cells") %>%
    record_seurat_samples(type = "nFeature_RNA") %>%
    record_seurat_samples(type = "nCount_RNA") %>%
    record_seurat_samples(type = "percent_mito", pattern = "^mt-") %>%
    record_seurat_samples(type = "percent_ribo", pattern = "^Rp[sl]")

  post_plots <- scatter_seurat(scd)

  retlist <- list(
      "pre_plots" = pre_plots,
      "scd" = filt_scd,
      "post_plots" = post_plots)
  return(retlist)
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
                                  verbose = TRUE) {
  sample_meta <- scd@misc[["sample_metadata"]]

  if (!is.null(pattern)) {
    scd[[type]] <- Seurat::PercentageFeatureSet(scd, pattern = pattern)
  }

  ## I think this should be smart enough to take a list.
  if (length(type) > 1) {
    count <- 0
    for (t in type) {
      count <- count + 1
      scd <- record_seurat_samples(scd, type = t, pattern = pattern[count],
                                   column_name = column_name, column_prefix = column_prefix,
                                   verbose = verbose)
    }
    return(scd)
  }

  scd_meta <- scd@meta.data
  if (type == "num_cells") {
    sample_meta[[type]] <- 0
    for (s in 1:nrow(sample_meta)) {
      sample_name <- rownames(sample_meta)[s]
      sample_meta[s, type] <- sum(scd_meta[["orig.ident"]] == sample_name)
    }
  } else {
    sample_meta <- skim_seurat_metadata(
        sample_meta, scd_meta, meta_query = type,
        group_column = "orig.ident", summary_query = "numeric.mean",
        column_name = column_name, column_prefix = column_prefix,
        verbose = verbose)
  }

  scd@misc[["sample_metadata"]] <- sample_meta
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
    xy <- set[p]
    if (is.null(scd[[x]])) {
      message("The x-axis data does not exist: ", x, ".")
      next
    }
    if (is.null(scd[[y]])) {
      message("The y-axis data does not exist: ", y, ".")
      next
    }
    plots[[name]] <- Seurat::FeatureScatter(scd, x, y)
  }
  return(plots)
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
                                 group_column = "orig.ident", summary_query = "numeric.mean",
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
      skim_tee(meta_query) %>%
      skim(meta_query) %>%
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
