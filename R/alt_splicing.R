## A group of functions to examine the results from alternative splicing
## methods. The hope is to remove some corner cases from tools like
## suppa/miso/etc and somewhat standardize the resulting outputs.

#' Given some psi and tpm data, make a pretty plot!
#'
#' This should take either a dataframe or filename for the psi data from suppa,
#' along with the same for the average log tpm data (acquired from suppa
#' diffSplice with --save_tpm_events)
#'
#' @param dpsi Table provided by suppa containing all the metrics.
#' @param tpm Table provided by suppa containing all the tpm values.
#' @param events List of event types to include.
#' @param psi Limit the set of included events by psi value?
#' @param sig_threshold Use this significance threshold.
#' @param label_type Choose a type of event to label.
#' @param alpha How see-through should the points be in the plot?
#' @return List containing the plot and some of the requisite data.
#' @seealso [plot_rmats()]
#' @examples
#'  \dontrun{
#'  suppa_plot <- plot_suppa(dpsi_file, tmp_file)
#' }
#' @export
plot_suppa <- function(dpsi, tpm, events = NULL, psi = NULL, sig_threshold = 0.05,
                       label_type = NULL, alpha = 0.7) {
  dpsi_data <- NULL
  if (class(psi) == "character") {
    dpsi_data <- read.table(dpsi, sep = "\t")
  } else if (class(dpsi) == "data.frame") {
    dpsi_data <- dpsi
  } else {
    stop("I only understand filenames and data frames, your psi are neither.")
  }
  rownames(dpsi_data) <- dpsi_data[[1]]
  dpsi_data[[1]] <- NULL
  colnames(dpsi_data) <- c("dpsi", "pvalue")

  tpm_data <- NULL
  if (class(tpm) == "character") {
    tpm_data <- read.table(tpm, sep = "\t")
  } else if (class(tpm) == "data.frame") {
    tpm_data <- tpm
  } else {
    stop("I only understand filenames and data frames, your tpms are neither.")
  }
  colnames(tpm_data) <- c("event", "avglogtpm")

  events_data <- NULL
  if (!is.null(events)) {
    if (class(events) == "character") {
      events_data <- read.table(events, sep = "\t", header = TRUE)
    } else if (class(events) == "data.frame") {
      events_data <- events
    } else {
      stop("I only understand filenames and data frames, your events are neither.")
    }
  }

  psi_data <- NULL
  if (!is.null(psi)) {
    if (class(psi) == "character") {
      psi_data <- read.table(psi, sep = "\t")
    } else if (class(psi) == "data.frame") {
      psi_data <- psi
    } else {
      stop("I only understand filenames and data frames, your events are neither.")
    }
  }
  colnames(psi_data) <- c(
    "denominator1", "denominator2", "denominator3",
    "numerator1", "numerator2", "numerator3", "numerator4", "numerator5", "numerator6")

  plotting_data <- merge(dpsi_data, tpm_data, by.x = "row.names", by.y = "event")
  rownames(plotting_data) <- plotting_data[["Row.names"]]
  plotting_data <- plotting_data[, -1]
  ## Now we have the basis for everything we are likely to plot.
  ## Add some categorizations of the data.
  plotting_data[["log10pval"]] <- -1.0 * log10(plotting_data[["pvalue"]])
  plotting_data[["adjp"]] <- p.adjust(plotting_data[["pvalue"]], method = "fdr")
  plotting_data[["psig"]] <- FALSE
  plotting_data[["adjpsig"]] <- FALSE
  plotting_data[["gene_name"]] <- gsub(x = rownames(plotting_data),
                                       pattern = "^(.*);.*$",
                                       replacement = "\\1")
  plotting_data[["category"]] <- gsub(x = rownames(plotting_data),
                                      pattern = "^(.*);(.{2}):.*$",
                                      replacement = "\\2")
  plotting_data[["coordinates"]] <- gsub(x = rownames(plotting_data),
                                         pattern = "^(.*);(.{2}):(.*)$",
                                         replacement = "\\3")
  plotting_data[["plot_cat"]] <- plotting_data[["category"]]
  plotting_data[["plot_cat"]] <- ifelse(
    test = plotting_data[["plot_cat"]] == "SE",
    yes = "Skipping exon",
    no = ifelse(
      test = plotting_data[["plot_cat"]] == "MX",
      yes = "Mutually exclusive exons",
      no = ifelse(
        test = plotting_data[["plot_cat"]] == "A5",
        yes = "Alternate 5 prime",
        no = ifelse(
          test = plotting_data[["plot_cat"]] == "A3",
          yes = "Alternate 3 prime",
          no = ifelse(
            test = plotting_data[["plot_cat"]] == "RI",
            yes = "Retained intron",
            no = ifelse(plotting_data[["plot_cat"]] == "AF",
                        yes = "Alternate first exon",
                        no = ifelse(plotting_data[["plot_cat"]] == "AL",
                                    yes = "Alternate last exon",
                                    no = "Unknown")))))))
  plotting_data[["category"]] <- plotting_data[["plot_cat"]]
  insig_idx <- plotting_data[["pvalue"]] > 0.05
  plotting_data[insig_idx, "plot_cat"] <- "Insignificant"
  ## If, somehow something is observed as unknown, make it insignificant.
  unknown_idx <- plotting_data[["plot_cat"]] == "Unknown"
  plotting_data[unknown_idx, "plot_cat"] <- "Insignificant"

  level_names <- c("Skipping exon", "Mutually exclusive exons", "Alternate 5 prime",
                   "Alternate 3 prime", "Retained intron", "Alternate first exon",
                   "Alternate last exon", "Insignificant")
  plotting_data[["category"]] <- factor(plotting_data[["category"]], levels = level_names,
                                        labels = level_names)
  plotting_data[["event"]] <- rownames(plotting_data)

  psig_idx <- plotting_data[["pvalue"]] <= sig_threshold
  plotting_data[psig_idx, "psig"] <- TRUE
  adjpsig_idx <- plotting_data[["adjp"]] <= sig_threshold
  plotting_data[adjpsig_idx, "adjpsig"] <- TRUE

  ## A quick volcano plot, which should be made prettier soon.
  sig_splicing_volplot <- ggplot(plotting_data,
                                 aes(x = .data[["dpsi"]], y = .data[["log10pval"]],
                                     color = .data[["psig"]])) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = c(-0.5, 0.5), size = 1) +
    ggplot2::geom_hline(yintercept = 1.3, size = 1)

  ## Now a somewhat more involved ma plot, first dropping the super-low tpm stuff
  plotting_subset_idx <- plotting_data[["avglogtpm"]] >= -1
  plotting_subset <- plotting_data[plotting_subset_idx, ]
  label_subset_idx <- plotting_subset[["psig"]] == TRUE
  label_subset <- plotting_subset[label_subset_idx, ]
  if (!is.null(label_type)) {
    type_subset_idx <- label_subset[["plot_cat"]] == label_type
    label_subset <- label_subset[type_subset_idx, ]
  }

  color_values <- c("Skipping exon" = "#92250E",
                    "Mutually exclusive exons" = "#717600",
                    "Alternate 5 prime" = "#5B095A",
                    "Alternate 3 prime" = "#1D6E72",
                    "Retained intron" = "#31326D",
                    "Alternate first exon" = "black",
                    "Alternate last exon" = "#003300",
                    "Insignificant" = "#666666")

  ## Now make a quick and dirty ma plot.
  sig_splicing_maplot <- ggplot(plotting_subset,
                                aes(x = .data[["avglogtpm"]], y = .data[["dpsi"]],
                                    color = .data[["plot_cat"]], fill = .data[["plot_cat"]])) +
    ggplot2::geom_point(alpha = alpha) +
    ggplot2::scale_shape_manual(values = 21) +
    ggplot2::scale_fill_manual(name = "Category",
                               guide = "legend", ## keep my preferred order.
                               breaks = levels(plotting_data[["category"]]),
                               values = color_values) +
    ggplot2::scale_color_manual(name = "Category",
                                values = color_values, ## keep my preferred order.
                                breaks = levels(plotting_data[["category"]]),
                                guide = ggplot2::guide_legend(override.aes = list(size = 5))) +
    ggrepel::geom_text_repel(data = label_subset,
                             show.legend = FALSE,
                             arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc")),
                             aes(x = .data[["avglogtpm"]], y = .data[["dpsi"]],
                                 label = .data[["gene_name"]])) +
    ggplot2::xlab("Average log(transcripts per million).") +
    ggplot2::ylab("Delta PSI calculated by Suppa.") +
    ggplot2::theme_bw(base_size = base_size)

  ## If available, add the event data to the table
  if (!is.null(events_data)) {
    events_data <- events_data[, c("event_id", "alternative_transcripts", "total_transcripts")]
    plotting_data <- merge(plotting_data,
                           events_data,
                           by.x = "row.names",
                           by.y = "event_id")
    rownames(plotting_data) <- plotting_data[["Row.names"]]
    plotting_data <- plotting_data[, -1]
  }
  ## If available, add all the psi data to the table
  if (!is.null(psi_data)) {
    plotting_data <- merge(plotting_data,
                           psi_data,
                           by = "row.names")
    rownames(plotting_data) <- plotting_data[["Row.names"]]
    plotting_data <- plotting_data[, -1]
  }

  retlist <- list(
    "volcano" = sig_splicing_volplot,
    "ma" = sig_splicing_maplot,
    "data" = plotting_data)
  return(retlist)
}

#' Take a set of results from suppa and attempt to write it to a pretty xlsx file.
#'
#' Suppa provides a tremendous amount of output, this attempts to standardize
#' those results and print them to an excel sheet.
#'
#' @param table Result table from suppa.
#' @param annotations Set of annotation data to include with the suppa result.
#' @param by_table Use this column to merge the annotations and data tables from
#'  the perspective of the data table.
#' @param by_annot Use this column to merge the annotations and data tables
#'  from the perspective of the annotations.
#' @param columns Choose a subset of columns to include, or leave the defaults.
#' @param excel Provide an excel file to write.
#' @return Data frame of the merged data.
#' @seealso [write_xlsx()]
#' @examples
#'  \dontrun{
#'  prettier_table <- write_suppa_table(suppa_result_file,
#'                                      annotations = gene_info,
#'                                      excel = "excel/pretty_suppa_table.xlsx")
#' }
#' @export
write_suppa_table <- function(table, annotations = NULL, by_table = "gene_name",
                              by_annot = "ensembl_gene_id",
                              columns = "default", excel = "excel/suppa_table.xlsx") {

  default_columns <- c("event", "dpsi", "pvalue", "adjp", "avglogtpm", "category",
                       "coordinates", "alternative_transcripts", "total_transcripts",
                       "denominator1", "denominator2", "denominator3", "numerator1",
                       "numerator2", "numerator3", "numerator4", "numerator5", "numerator6")
  full_table <- data.frame()
  if (is.null(annotations)) {
    full_table <- data.table::as.data.table(table)
  } else {
    annot <- data.table::as.data.table(annotations)
    annot[[by_annot]] <- make.names(annot[[by_annot]], unique = TRUE)
    tab <- data.table::as.data.table(table)
    full_table <- merge(annot, tab, by.x = by_annot, by.y = by_table, all.y = TRUE)
  }

  xls_data <- as.data.frame(full_table)
  ## Now coerce numeric columns
  xlsx_result <- write_xlsx(data = xls_data, excel = excel)
  return(xls_data)
}

#' Given some psi and tpm data from rMATS, make a pretty plot!
#'
#' This should take either a dataframe or filename for the psi data from rMATS.
#' This was mostly copy/pasted from plot_suppa().
#'
#' @param se Table of skipped exon data from rmats.
#' @param a5ss Table of alternate 5p exons.
#' @param a3ss Table of alternate 3p exons.
#' @param mxe Table of alternate exons.
#' @param ri Table of retained introns.
#' @param sig_threshold Use this significance threshold.
#' @param dpsi_threshold Use a delta threshold.
#' @param label_type Choose a type of event to label.
#' @param alpha How see-through should the points be in the plot?
#' @return List containing the plot and some of the requisite data.
#' @seealso [plot_supps()]
#' @examples
#'  \dontrun{
#'  rmats_plot <- plot_rmats(se_table, a5_table, a3_table)
#' }
#' @export
plot_rmats <- function(se = NULL, a5ss = NULL, a3ss = NULL, mxe = NULL, ri = NULL,
                       sig_threshold = 0.05, dpsi_threshold = 0.7,
                       label_type = NULL, alpha = 0.7) {
  if (is.null(se) && is.null(a5ss) && is.null(a3ss) &&
        is.null(mxe) && is.null(ri)) {
    stop("No data was provided.")
  }

  se_data <- data.frame()
  if (class(se) == "character") {
    se_data <- read.table(se, header = TRUE)
  } else if (class(se) == "data.frame" | class(se) == "NULL") {
    se_data <- se_data
  } else {
    stop("I only understand filenames and data frames, your psi are neither.")
  }

  a5ss_data <- data.frame()
  if (class(a5ss) == "character") {
    a5ss_data <- read.table(a5ss, header = TRUE)
  } else if (class(a5ss)[1] == "data.frame" || class(a5ss)[1] == "NULL") {
    a5ss_data <- a5ss_data
  } else {
    stop("I only understand filenames and data frames, your psi are neither.")
  }

  a3ss_data <- data.frame()
  if (class(a3ss) == "character") {
    a3ss_data <- read.table(a3ss, header = TRUE)
  } else if (class(a3ss)[1] == "data.frame" || class(a3ss)[1] == "NULL") {
    a3ss_data <- a3ss_data
  } else {
    stop("I only understand filenames and data frames, your psi are neither.")
  }

  mxe_data <- data.frame()
  if (class(mxe) == "character") {
    mxe_data <- read.table(mxe, header = TRUE)
  } else if (class(mxe)[1] == "data.frame" || class(mxe)[1] == "NULL") {
    mxe_data <- mxe_data
  } else {
    stop("I only understand filenames and data frames, your psi are neither.")
  }

  ri_data <- data.frame()
  if (class(ri) == "character") {
    ri_data <- read.table(ri, header = TRUE)
  } else if (class(ri)[1] == "data.frame" || class(ri)[1] == "NULL") {
    ri_data <- ri_data
  } else {
    stop("I only understand filenames and data frames, your psi are neither.")
  }

  all_data <- data.frame()
  kept_columns <- c("GeneID", "PValue", "FDR", "IncLevel1",
                    "IncLevel2", "IncLevelDifference")
  new_colnames <- c("gene_id", "pvalue", "adjp", "level1", "level2", "dpsi")

  if (nrow(se_data > 0)) {
    se_subset <- se_data[, kept_columns]
    colnames(se_subset) <- new_colnames
    se_subset[["event"]] <- "SE"
    all_data <- rbind(all_data, se_subset)
  }
  if (nrow(a5ss_data) > 0) {
    a5ss_subset <- a5ss_data[, kept_columns]
    colnames(a5ss_subset) <- new_colnames
    a5ss_subset[["event"]] <- "A5"
    all_data <- rbind(all_data, a5ss_subset)
  }
  if (nrow(a3ss_data) > 0) {
    a3ss_subset <- a3ss_data[, kept_columns]
    colnames(a3ss_subset) <- new_colnames
    a3ss_subset[["event"]] <- "A3"
    all_data <- rbind(all_data, a3ss_subset)
  }
  if (nrow(mxe_data) > 0) {
    mxe_subset <- mxe_data[, kept_columns]
    colnames(mxe_subset) <- new_colnames
    mxe_subset[["event"]] <- "MX"
    all_data <- rbind(all_data, mxe_subset)
  }
  if (nrow(ri_data) > 0) {
    ri_subset <- ri_data[, kept_columns]
    colnames(ri_subset) <- new_colnames
    ri_subset[["event"]] <- "RI"
    all_data <- rbind(all_data, ri_subset)
  }

  all_data <- data.table::as.data.table(all_data)
  ## Adding stub variables to make dplyr/tidyr NSE operations not throw warnings
  ## when running R CMD check
  l1a <- l1b <- l2a <- l2b <- id <- l1mean <- l2mean <- NULL
  plotting_data <- all_data %>%
    tidyr::separate("level1", c("l1a", "l1b"), "\\,") %>%
    tidyr::separate("level2", c("l2a", "l2b"), "\\,")
  plotting_data[["l1a"]] <- suppressWarnings(as.numeric(plotting_data[["l1a"]]))
  plotting_data[["l1b"]] <- suppressWarnings(as.numeric(plotting_data[["l1b"]]))
  plotting_data[["l2a"]] <- suppressWarnings(as.numeric(plotting_data[["l2a"]]))
  plotting_data[["l2b"]] <- suppressWarnings(as.numeric(plotting_data[["l2b"]]))
  plotting_data[["id"]] <- rownames(plotting_data)
  ## suppressWarnings(plotting_data[, `:=` (l1mean = mean(c(l1a, l1b), na.rm = TRUE)), by = id])
  plotting_data <- plotting_data %>%
    group_by(id) %>%
    dplyr::mutate(l1mean = mean(c(l1a, l1b), na.rm = TRUE),
                  l2mean = mean(c(l2a, l2b), na.rm = TRUE),
                  all_mean = mean(c(l1mean, l2mean), na.rm = TRUE))
  ## plotting_data[, `:=` (l2mean = mean(c(l2a, l2b), na.rm = TRUE)), by = id]
  ## plotting_data[, `:=` (all_mean = mean(c(l1mean, l2mean), na.rm = TRUE)), by = id]
  plotting_data[["check"]] <- plotting_data[["l1mean"]] - plotting_data[["l2mean"]]
  test <- all.equal(plotting_data[["check"]],
                    plotting_data[["dpsi"]], tolerance = 0.01)
  plotting_data[["log10pval"]] <- -1.0 * log10(plotting_data[["pvalue"]])
  plotting_data[["log10adjpval"]] <- -1.0 * log10(plotting_data[["adjp"]])

  plotting_data[["psig"]] <- FALSE
  plotting_data[["adjpsig"]] <- FALSE

  plotting_data[["plot_cat"]] <- plotting_data[["event"]]
  plotting_data[["plot_cat"]] <- ifelse(
    test = plotting_data[["plot_cat"]] == "SE",
    yes = "Skipping exon",
    no = ifelse(
      test = plotting_data[["plot_cat"]] == "MX",
      yes = "Mutually exclusive exons",
      no = ifelse(
        test = plotting_data[["plot_cat"]] == "A5",
        yes = "Alternate 5 prime",
        no = ifelse(
          test = plotting_data[["plot_cat"]] == "A3",
          yes = "Alternate 3 prime",
          no = ifelse(
            test = plotting_data[["plot_cat"]] == "RI",
            yes = "Retained intron",
            no = ifelse(plotting_data[["plot_cat"]] == "AF",
                        yes = "Alternate first exon",
                        no = ifelse(plotting_data[["plot_cat"]] == "AL",
                                    yes = "Alternate last exon",
                                    no = "Unknown")))))))
  plotting_data[["category"]] <- plotting_data[["plot_cat"]]
  insig_idx <- plotting_data[["adjp"]] > sig_threshold
  plotting_data[insig_idx, "plot_cat"] <- "Insignificant"
  ## If, somehow something is observed as unknown, make it insignificant.
  unknown_idx <- plotting_data[["plot_cat"]] == "Unknown"
  plotting_data[unknown_idx, "plot_cat"] <- "Insignificant"

  level_names <- c("Skipping exon", "Mutually exclusive exons", "Alternate 5 prime",
                   "Alternate 3 prime", "Retained intron", "Alternate first exon",
                   "Alternate last exon", "Insignificant")
  plotting_data[["category"]] <- factor(plotting_data[["category"]], levels = level_names,
                                        labels = level_names)
  plotting_data[["event"]] <- rownames(plotting_data)

  psig_idx <- plotting_data[["adjp"]] <= sig_threshold
  plotting_data[psig_idx, "psig"] <- TRUE
  adjpsig_idx <- plotting_data[["adjp"]] <= sig_threshold
  plotting_data[adjpsig_idx, "adjpsig"] <- TRUE

  ## A quick volcano plot, which should be made prettier soon.
  sig_splicing_volplot <- ggplot(plotting_data,
                                 aes(x = .data[["dpsi"]], y = .data[["log10pval"]],
                                     color = .data[["psig"]])) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = c(-0.5, 0.5), linewidth = 1) +
    ggplot2::geom_hline(yintercept = 1.3, linewidth = 1)

  ## Now a somewhat more involved ma plot, first dropping the super-low tpm stuff
  label_subset_idx <- plotting_data[["psig"]] == TRUE
  label_subset <- plotting_data[label_subset_idx, ]
  label_subset_idx <- abs(label_subset[["dpsi"]]) > dpsi_threshold ## |
  ##    abs(label_subset[["all_mean"]] - 0.5) > (dpsi_threshold / 2)
  label_subset <- label_subset[label_subset_idx, ]
  ##label_subset_idx <- label_subset[["plot_cat"]] != "Insignificant"
  ##label_subset <- label_subset[label_subset_idx, ]
  if (!is.null(label_type)) {
    type_subset_idx <- label_subset[["plot_cat"]] == label_type
    label_subset <- label_subset[type_subset_idx, ]
  }

  color_values <- c("Skipping exon" = "#92250E",
                    "Mutually exclusive exons" = "#717600",
                    "Alternate 5 prime" = "#5B095A",
                    "Alternate 3 prime" = "#1D6E72",
                    "Retained intron" = "#31326D",
                    "Alternate first exon" = "black",
                    "Alternate last exon" = "#003300",
                    "Insignificant" = "#666666")

  ## Now make a quick and dirty ma plot.
  sig_splicing_maplot <- ggplot(
    plotting_data,
    aes(x = .data[["all_mean"]], y = .data[["dpsi"]],
        color = .data[["plot_cat"]], fill = .data[["plot_cat"]])) +
    ggplot2::geom_point(alpha = alpha) +
    ggplot2::scale_shape_manual(values = 21) +
    ggplot2::scale_fill_manual(name = "Category",
                               guide = "legend", ## keep my preferred order.
                               breaks = levels(plotting_data[["category"]]),
                               values = color_values) +
    ggplot2::scale_color_manual(name = "Category",
                                values = color_values, ## keep my preferred order.
                                breaks = levels(plotting_data[["category"]]),
                                guide = ggplot2::guide_legend(override.aes = list(size = 5))) +
    ggrepel::geom_text_repel(
      data = label_subset,
      show.legend = FALSE,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc")),
      aes(x = .data[["all_mean"]], y = .data[["dpsi"]], label = .data[["gene_id"]])) +
    ggplot2::xlab("Average Inclusion.") +
    ggplot2::ylab("Delta PSI calculated by rMATS.") +
    ggplot2::theme_bw(base_size = base_size)

  retlist <- list(
    "volcano" = sig_splicing_volplot,
    "ma" = sig_splicing_maplot,
    "data" = plotting_data)
  return(retlist)
}

## EOF
