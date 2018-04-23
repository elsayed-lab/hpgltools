#' Given some psi and tpm data, make a pretty plot!
#'
#' This should take either a dataframe or filename for the psi data from suppa,
#' along with the same for the average log tpm data (acquired from suppa
#' diffSplice with --save_tpm_events)
#' @export
plot_suppa <- function(dpsi, tpm, events=NULL, psi=NULL, sig_threshold=0.05,
                       label_type=NULL, alpha=0.7) {
  dpsi_data <- NULL
  if (class(psi) == "character") {
    dpsi_data <- read.table(dpsi, sep="\t")
  } else if (class(dpsi) == "data.frame") {
    dpsi_data <- dpsi
  } else {
    stop("I only understand filenames and data frames, your psi are neither.")
  }
  colnames(dpsi_data) <- c("dpsi", "pvalue")

  tpm_data <- NULL
  if (class(tpm) == "character") {
    tpm_data <- read.table(tpm, sep="\t")
  } else if (class(tpm) == "data.frame") {
    tpm_data <- tpm
  } else {
    stop("I only understand filenames and data frames, your tpms are neither.")
  }
  colnames(tpm_data) <- c("event", "avglogtpm")

  events_data <- NULL
  if (!is.null(events)) {
    if (class(events) == "character") {
      events_data <- read.table(events, sep="\t", header=TRUE)
    } else if (class(events) == "data.frame") {
      events_data <- events
    } else {
      stop("I only understand filenames and data frames, your events are neither.")
    }
  }

  psi_data <- NULL
  if (!is.null(psi)) {
    if (class(psi) == "character") {
      psi_data <- read.table(psi, sep="\t")
    } else if (class(psi) == "data.frame") {
      psi_data <- psi
    } else {
      stop("I only understand filenames and data frames, your events are neither.")
    }
  }
  colnames(psi_data) <- c(
    "denominator1", "denominator2", "denominator3",
    "numerator1", "numerator2", "numerator3", "numerator4", "numerator5", "numerator6")

  plotting_data <- merge(dpsi_data, tpm_data, by.x="row.names", by.y="event")
  rownames(plotting_data) <- plotting_data[["Row.names"]]
  plotting_data <- plotting_data[, -1]
  ## Now we have the basis for everything we are likely to plot.
  ## Add some categorizations of the data.
  plotting_data[["log10pval"]] <- -1.0 * log10(plotting_data[["pvalue"]])
  plotting_data[["adjp"]] <- p.adjust(plotting_data[["pvalue"]], method="fdr")
  plotting_data[["psig"]] <- FALSE
  plotting_data[["adjpsig"]] <- FALSE
  plotting_data[["gene_name"]] <- gsub(x=rownames(plotting_data),
                                  pattern="^(.*);.*$",
                                  replacement="\\1")
  plotting_data[["category"]] <- gsub(x=rownames(plotting_data),
                                      pattern="^(.*);(.{2}):.*$",
                                      replacement="\\2")
  plotting_data[["coordinates"]] <- gsub(x=rownames(plotting_data),
                                         pattern="^(.*);(.{2}):(.*)$",
                                         replacement="\\3")
  plotting_data[["plot_cat"]] <- plotting_data[["category"]]
  plotting_data[["plot_cat"]] <- ifelse(
    test=plotting_data[["plot_cat"]] == "SE",
    yes="Skipping exon",
    no=ifelse(
      test=plotting_data[["plot_cat"]] == "MX",
      yes="Mutually exclusive exons",
      no=ifelse(
        test=plotting_data[["plot_cat"]] == "A5",
        yes="Alternate 5 prime",
        no=ifelse(
          test=plotting_data[["plot_cat"]] == "A3",
          yes="Alternate 3 prime",
          no=ifelse(
            test=plotting_data[["plot_cat"]] == "RI",
            yes="Retained intron",
            no=ifelse(plotting_data[["plot_cat"]] == "AF",
                      yes="Alternate first exon",
                      no=ifelse(plotting_data[["plot_cat"]] == "AL",
                                yes="Alternate last exon",
                                no="Unknown")))))))
  plotting_data[["category"]] <- plotting_data[["plot_cat"]]
  insig_idx <- plotting_data[["pvalue"]] > 0.05
  plotting_data[insig_idx, "plot_cat"] <- "Insignificant"
  ## If, somehow something is observed as unknown, make it insignificant.
  unknown_idx <- plotting_data[["plot_cat"]] == "Unknown"
  plotting_data[unknown_idx, "plot_cat"] <- "Insignificant"

  level_names <- c("Skipping exon", "Mutually exclusive exons", "Alternate 5 prime",
                   "Alternate 3 prime", "Retained intron", "Alternate first exon",
                   "Alternate last exon", "Insignificant")
  plotting_data[["category"]] <- factor(plotting_data[["category"]], levels=level_names,
                                        labels=level_names)
  plotting_data[["event"]] <- rownames(plotting_data)

  psig_idx <- plotting_data[["pvalue"]] <= sig_threshold
  plotting_data[psig_idx, "psig"] <- TRUE
  adjpsig_idx <- plotting_data[["adjp"]] <= sig_threshold
  plotting_data[adjpsig_idx, "adjpsig"] <- TRUE

  ## A quick volcano plot, which should be made prettier soon.
  sig_splicing_volplot <- ggplot(plotting_data, aes_string(x="dpsi", y="log10pval", color="psig")) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept=c(-0.5, 0.5), size=1) +
    ggplot2::geom_hline(yintercept=1.3, size=1)

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
                                aes_string(x="avglogtpm", y="dpsi",
                                           color="plot_cat", fill="plot_cat")) +
    ggplot2::geom_point(alpha=alpha) +
    ggplot2::scale_shape_manual(values=21) +
    ggplot2::scale_fill_manual(name="Category",
                               guide="legend",
                               breaks=levels(plotting_data[["category"]]), ## keep my preferred order.
                               values=color_values) +
    ggplot2::scale_color_manual(name="Category",
                                values=color_values,
                                breaks=levels(plotting_data[["category"]]), ## keep my preferred order.
                                guide=ggplot2::guide_legend(override.aes=list(size=5))) +
    ggrepel::geom_text_repel(data=label_subset,
                             show.legend=FALSE,
                             arrow=ggplot2::arrow(length=ggplot2::unit(0.01, "npc")),
                             aes_string(x="avglogtpm", y="dpsi", label="gene_name")) +
    ggplot2::xlab("Average log(transcripts per million).") +
    ggplot2::ylab("Delta PSI calculated by Suppa.") +
    ggplot2::theme_bw(base_size=base_size)

  ## If available, add the event data to the table
  if (!is.null(events_data)) {
    events_data <- events_data[, c("event_id", "alternative_transcripts", "total_transcripts")]
    plotting_data <- merge(plotting_data,
                           events_data,
                           by.x="row.names",
                           by.y="event_id")
    rownames(plotting_data) <- plotting_data[["Row.names"]]
    plotting_data <- plotting_data[, -1]
  }
  ## If available, add all the psi data to the table
  if (!is.null(psi_data)) {
    plotting_data <- merge(plotting_data,
                           psi_data,
                           by="row.names")
    rownames(plotting_data) <- plotting_data[["Row.names"]]
    plotting_data <- plotting_data[, -1]
  }

  retlist <- list(
    "volcano" = sig_splicing_volplot,
    "ma" = sig_splicing_maplot,
    "data" = plotting_data)
  return(retlist)
}

write_suppa_table <- function(table, annotations=NULL, by_table="gene_name", by_annot="ensembl_gene_id",
                              columns="default", excel="excel/suppa_table.xlsx") {
  default_columns <- c("ensembl_gene_id", "version", "hgnc_symbol", "description",
                       "gene_biotype", "cds_length", "chromosome_name", "strand",
                       "start_position", "end_position", "event", "dpsi", "pvalue", "adjp",
                       "avglogtpm", "category", "coordinates", "alternative_transcripts",
                       "total_transcripts", "denominator1", "denominator2", "denominator3",
                       "numerator1", "numerator2", "numerator3", "numerator4", "numerator5",
                       "numerator6")
  full_table <- data.frame()
  if (is.null(annotations)) {
    full_table <- data.table::as.data.table(table)
  } else {
    annot <- data.table::as.data.table(annotations)
    annot[[by_annot]] <- make.names(annot[[by_annot]], unique=TRUE)
    tab <- data.table::as.data.table(table)
    full_table <- merge(annot, tab, by.x=by_annot, by.y=by_table, all.y=TRUE)
  }

  chosen_columns <- default_columns
  if (columns != "default") {
    chosen_columns <- columns
  }
  full_table <- full_table[, chosen_columns, with=FALSE]
  xls_data <- as.data.frame(full_table)
  ## Now coerce numeric columns
  xlsx_result <- write_xls(data=xls_data, excel=excel)
}
