#' Run simple_gprofiler on every table from extract_significant_genes()
#'
#' @param sig Result from extract_significant_genes
#' @param according_to Use this result type for the gprofiler searches.
#' @param together Concatenate the up/down genes into one set?
#' @param ... Arguments to pass to simple_gprofiler().
#' @export
all_gprofiler <- function(sig, according_to = "deseq", together = FALSE, ...) {
  ret <- list()
  input_up <- list()
  input_down <- list()
  source <- "significant"
  ## Check if this came from extract_significant_genes or extract_abundant_genes.
  if (!is.null(sig[[according_to]][["ups"]])) {
    input_up <- sig[[according_to]][["ups"]]
    input_down <- sig[[according_to]][["downs"]]
  } else if (!is.null(sig[["abundances"]])) {
    source <- "abundance"
    input_up <- sig[["abundances"]][[according_to]][["high"]]
    input_down <- sig[["abundances"]][[according_to]][["low"]]
  } else {
    stop("I do not understand this input.")
  }

  sig_names <- names(input_up)
  for (i in seq_along(sig_names)) {
    slept <- Sys.sleep(1)
    name <- sig_names[i]
    retname_up <- paste0(name, "_up")
    retname_down <- paste0(name, "_down")
    up <- input_up[[name]]
    down <- input_down[[name]]
    up_elements <- 0
    down_elements <- 0
    if (source == "abundance") {
      up <- names(up)
      down <- names(down)
      up_elements <- length(up)
      down_elements <- length(down)
    } else {
      up_elements <- nrow(up)
      down_elements <- nrow(down)
    }
    if (isTRUE(together)) {
      if (source == "abundance") {
        up <- c(up, down)
        up_elements <- up_elements + down_elements
        down <- c()
        down_elements <- 0
      } else {
        up <- rbind(up, down)
        up_elements <- nrow(up)
        down <- data.frame()
        down_elements <- 0
      }
    }
    if (up_elements > 0) {
      ret[[retname_up]] <- simple_gprofiler(up, ...)
    } else {
      ret[[retname_up]] <- NULL
    }
    if (down_elements > 0) {
      ret[[retname_down]] <- simple_gprofiler(down, ...)
    } else {
      ret[[retname_down]] <- NULL
    }
  }
  return(ret)
}

#' Run searches against the web service g:Profiler.
#'
#' This is the beginning of a reimplementation to use gprofiler2.  However,
#' AFAICT gprofiler2 does not yet actually work for anything other than their GO
#' data.
#'
#' @param sig_genes Guess!  The set of differentially expressed/interesting
#'  genes.
#' @param species Organism supported by gprofiler.
#' @param convert Use gProfileR's conversion utility?
#' @param first_col First place used to define the order of 'significant'.
#' @param second_col If that fails, try a second column.
#' @param do_go Perform GO search?
#' @param do_kegg Perform KEGG search?
#' @param do_reactome Perform reactome search?
#' @param do_mi Do miRNA search?
#' @param do_tf Search for transcription factors?
#' @param do_corum Do corum search?
#' @param do_hp Do the hp search?
#' @param do_hpa Do the hpa search?
#' @param do_wp Do the wp search?
#' @param significant Only return the statistically significant hits?
#' @param exclude_iea Passed directly to gprofiler2.
#' @param do_under Perform under-representation search?
#' @param evcodes Get the set of evcodes in the data?  This makes it take
#'  longer.
#' @param threshold p-value 'significance' threshold.
#' @param adjp Method to adjust p-values.
#' @param domain_scope Passed to gprofiler2.
#' @param bg Background genes.
#' @param pseudo_gsea Is the data in a ranked order by significance?
#' @param id_col Which column in the table should be used for gene ID
#'  crossreferencing?  gProfiler uses Ensembl ids.  So if you have a table of
#'  entrez or whatever, translate it!
#' @param excel Print the results to an excel file?
#' @return a list of results for go, kegg, reactome, and a few more.
#' @seealso [gProfiler]
#' @examples
#' \dontrun{
#'  gprofiler_is_nice_and_easy <- simple_gprofiler(genes, species='mmusculus')
#' }
#' @export
simple_gprofiler2 <- function(sig_genes, species = "hsapiens", convert = TRUE,
                              first_col = "deseq_logfc", second_col = "logfc", do_go = TRUE,
                              do_kegg = TRUE, do_reactome = TRUE, do_mi = TRUE, do_tf = TRUE,
                              do_corum = TRUE, do_hp = TRUE, do_hpa = TRUE, do_wp = TRUE,
                              significant = TRUE, exclude_iea = FALSE, do_under = FALSE,
                              evcodes = TRUE, threshold = 0.05, adjp = "g_SCS",
                              domain_scope = "annotated", bg = NULL,
                              pseudo_gsea = TRUE, id_col = "row.names", excel = NULL) {
  gene_list <- NULL
  if ("character" %in% class(sig_genes)) {
    gene_ids <- sig_genes
  } else {
    if (!is.null(sig_genes[[first_col]])) {
      gene_list <- sig_genes[order(-sig_genes[[first_col]]), ]
      pseudo_gsea <- TRUE
    } else if (!is.null(sig_genes[[second_col]])) {
      gene_list <- sig_genes[order(-sig_genes[[second_col]]), ]
      pseudo_gsea <- TRUE
    }
    gene_ids <- NULL
    if (is.null(id_col)) {
      id_col <- "ID"
    }
    if (id_col == "row.names") {
      gene_ids <- rownames(gene_list)
    } else {
      gene_ids <- gene_list[[id_col]]
    }
  }

  retlst <- list()
  if (isTRUE(do_go)) {
    retlst[["GO"]] <- data.frame()
  }
  if (isTRUE(do_kegg)) {
    retlst[["KEGG"]] <- data.frame()
  }
  if (isTRUE(do_reactome)) {
    retlst[["REAC"]] <- data.frame()
  }
  if (isTRUE(do_wp)) {
    retlst[["WP"]] <- data.frame()
  }
  if (isTRUE(do_tf)) {
    retlst[["TF"]] <- data.frame()
  }
  if (isTRUE(do_mi)) {
    retlst[["MIRNA"]] <- data.frame()
  }
  if (isTRUE(do_hpa)) {
    retlst[["HPA"]] <- data.frame()
  }
  if (isTRUE(do_corum)) {
    retlst[["CORUM"]] <- data.frame()
  }
  if (isTRUE(do_hp)) {
    retlst[["HP"]] <- data.frame()
  }

  if (sum(grepl(pattern = "gene:", x = gene_ids)) > 0) {
    warning("Hey, it looks like you forgot to strip off the htseq prefix for the gene IDs.")
    gene_ids <- gsub(x = gene_ids, pattern = "gene:", replacement = "")
  }

  type_names <- names(retlst)
  retlst[["input"]] <- sig_genes

  interactive_plots <- list()
  gost_plots <- list()
  gost_links <- list()
  sig_tables <- list()
  types <- c("GO", "KEGG", "REAC", "WP", "TF", "MIRNA", "HPA", "CORUM", "HP")
  num_hits <- rep(0, length(types))
  names(num_hits) <- types
  for (t in seq_along(type_names)) {
    type <- type_names[t]
    mesg("Performing gProfiler ", type, " search of ",
         length(gene_ids), " genes against ", species, ".")
    Sys.sleep(1)
    ## To avoid the error: "'names' attribute [14] must be the same length as
    ## the vector [1]"
    gene_ids <- as.vector(gene_ids)
    a_result <- try(gprofiler2::gost(
                                    query = gene_ids,
                                    organism = species,
                                    evcodes = evcodes,
                                    significant = significant,
                                    ordered_query = pseudo_gsea,
                                    user_threshold = threshold,
                                    correction_method = adjp,
                                    domain_scope = domain_scope,
                                    custom_bg = bg,
                                    sources = type))
    a_df <- data.frame(stringsAsFactors = FALSE)
    if ("try-error" %in% class(a_result)) {
      mesg("The ", type, " method failed for this organism.")
    } else if (is.null(a_result)) {
      mesg("There was no result for the ", type, " search.")
    } else {
      a_df <- a_result[["result"]]
      sig_idx <- a_df[["p_value"]] <= threshold
      sig_df <- a_df[sig_idx, ]
      mesg(type, " search found ", nrow(sig_df), " hits.")
      num_hits[[type]] <- nrow(sig_df)
      sig_tables[[type]] <- sig_df
      gost_links[[type]] <- gprofiler2::gost(query = gene_ids, organism = species,
                                             evcodes = evcodes, significant = significant,
                                             ordered_query = pseudo_gsea,
                                             user_threshold = threshold,
                                             correction_method = adjp,
                                             domain_scope = domain_scope,
                                             custom_bg = bg, sources = type,
                                             as_short_link = TRUE)
      interactive_plots[[type]] <- try(
          gprofiler2::gostplot(a_result, capped = TRUE, interactive = TRUE), silent = TRUE)
      gost_plots[[type]] <- try(
          gprofiler2::gostplot(a_result, capped = FALSE, interactive = FALSE), silent = TRUE)
    }
    retlst[[type]] <- a_df
  } ## End iterating over the set of default sources.

  retlst[["interactive_plots"]] <- interactive_plots
  retlst[["num_hits"]] <- num_hits
  retlst[["gost_plots"]] <- gost_plots
  retlst[["gost_links"]] <- gost_links
  retlst[["significant"]] <- sig_tables
  retlst[["pvalue_plots"]] <- try(plot_gprofiler2_pval(retlst))
  if (!is.null(excel)) {
    mesg("Writing data to: ", excel, ".")
    excel_ret <- sm(try(write_gprofiler_data(retlst, excel = excel)))
    retlst[["excel"]] <- excel_ret
    mesg("Finished writing data.")
  }
  for (t in seq_along(type_names)) {
    type <- type_names[t]
    type_name <- paste0(type, "_enrich")
    ## Note to self, now that I think about it I think gprofiler2 provides its own p-adjustment.
    retlst[[type_name]] <- gprofiler2enrich(retlst, ontology = type,
                                             cutoff = threshold)
  }

  class(retlst) <- c("gprofiler_result", "list")
  return(retlst)
}

#' Redirect users to simple_gprofiler2
#'
#' @param ... Arguments passed to simple_gprofiler2()
#' @export
simple_gprofiler <- function(...) {
  simple_gprofiler2(...)
}

#' Run searches against the web service g:Profiler.
#'
#' Thank you Ginger for showing me your thesis, gProfiler is pretty cool!
#'
#' @param sig_genes Guess!  The set of differentially expressed/interesting
#'  genes.
#' @param species Organism supported by gprofiler.
#' @param convert Use gProfileR's conversion utility?
#' @param first_col First place used to define the order of 'significant'.
#' @param second_col If that fails, try a second column.
#' @param do_go Perform GO search?
#' @param do_kegg Perform KEGG search?
#' @param do_reactome Perform reactome search?
#' @param do_mi Do miRNA search?
#' @param do_tf Search for transcription factors?
#' @param do_corum Do corum search?
#' @param do_hp Do the hp search?
#' @param significant Only return the statistically significant hits?
#' @param pseudo_gsea Is the data in a ranked order by significance?
#' @param id_col Which column in the table should be used for gene ID
#'  crossreferencing?  gProfiler uses Ensembl ids.  So if you have a table of
#'  entrez or whatever, translate it!
#' @param excel Print the results to an excel file?
#' @return List of results for go, kegg, reactome, and a few more.
#' @seealso [gProfiler]
#' @examples
#' \dontrun{
#'  gprofiler_is_nice_and_easy <- simple_gprofiler(genes, species='mmusculus')
#' }
#' @export
simple_gprofiler_old <- function(sig_genes, species = "hsapiens", convert = TRUE,
                             first_col = "logFC", second_col = "limma_logfc", do_go = TRUE,
                             do_kegg = TRUE, do_reactome = TRUE, do_mi = TRUE, do_tf = TRUE,
                             do_corum = TRUE, do_hp = TRUE, significant = TRUE,
                             pseudo_gsea = TRUE, id_col = "row.names", excel = NULL) {
  gene_list <- NULL
  if (class(sig_genes)[1] == "AsIs") {
    gene_ids <- as.character(sig_genes)
  } else if (class(sig_genes)[1] == "character") {
    gene_ids <- sig_genes
  } else if (class(sig_genes)[1] == "numeric") {
    gene_ids <- names(sig_genes)
  } else {
    if (!is.null(sig_genes[[first_col]])) {
      gene_list <- sig_genes[order(-sig_genes[[first_col]]), ]
      pseudo_gsea <- TRUE
    } else if (!is.null(sig_genes[[second_col]])) {
      gene_list <- sig_genes[order(-sig_genes[[second_col]]), ]
      pseudo_gsea <- TRUE
    }
    gene_ids <- NULL
    if (is.null(id_col)) {
      id_col <- "ID"
    }
    if (id_col == "row.names") {
      gene_ids <- rownames(gene_list)
    } else {
      gene_ids <- gene_list[[id_col]]
    }
  }

  retlst <- list(
      "GO" = do_go,
      "KEGG" = do_kegg,
      "REAC" = do_reactome,
      "MI" = do_mi,
      "TF" = do_tf,
      "CORUM" = do_corum,
      "HP" = do_hp)
  type_names <- names(retlst)
  gene_ids <- as.vector(gene_ids)
  for (t in seq_along(type_names)) {
    type <- type_names[t]
    if (isTRUE(retlst[[type]])) {
      message("Performing gProfiler ", type, " search of ",
              length(gene_ids), " genes against ", species, ".")
      Sys.sleep(3)
      a_result <- suppressWarnings(
          try(gProfileR::gprofiler(
                             query = gene_ids,
                             organism = species,
                             significant = significant,
                             ordered_query = pseudo_gsea,
                             src_filter = type), silent = TRUE))
      if (class(a_result)[1] != "try-error") {
        retlst[[type]] <- a_result
        message(type, " search found ", nrow(a_result), " hits.")
      }
    } else {
      retlst[[type]] <- data.frame()
    }
  }
  names(retlst) <- tolower(names(retlst))
  retlst[["input"]] <- sig_genes
  retlst[["pvalue_plots"]] <- try(plot_gprofiler_pval(retlst), silent = TRUE)
  if (!is.null(excel)) {
    message("Writing data to: ", excel, ".")
    excel_ret <- sm(try(write_gprofiler_data(retlst, excel = excel)))
    retlst[["excel"]] <- excel_ret
    message("Finished writing data.")
  }
  class(retlst) <- c("gprofiler_result", "list")
  return(retlst)
}

#' Recast gProfiler data to the output class produced by clusterProfiler.
#'
#' I would like to use the various clusterProfiler plots more easily.
#' Therefore I figured it would be advantageous to coerce the various
#' outputs from gprofiler and friends into the datastructure produced by
#' clusterProfiler.
#'
#' @param retlst Output from simple_gprofiler()
#' @param ontology Category type to extract, currently only GO?
#' @param cutoff Use a p-value cutoff to get only the significant
#'  categories?
#' @param organism Set the orgdb organism name?
#' @param padjust_method what it says on the tin.
gprofiler2enrich <- function(retlst, ontology = "MF", cutoff = 1,
                              organism = NULL, padjust_method = "BH") {
  interesting <- retlst[[ontology]]

  sig_genes <- c()
  sig_genes_input <- retlst[["input"]]
  if (class(sig_genes_input)[1] == "character") {
    sig_genes <- sig_genes_input
  } else if ("data.frame" %in% class(sig_genes_input)) {
    sig_genes <- rownames(retlst[["input"]])
  } else {
    stop("I do not know this input data type when extracting the input genes.")
  }
  if (is.null(interesting)) {
    return(NULL)
  }
  if (nrow(interesting) == 0) {
    return(NULL)
  }

  bg_genes <- sum(!duplicated(sort(interesting[["term_id"]])))
  interesting[["tmp"]] <- bg_genes
  interesting[["adjusted"]] <- p.adjust(interesting[["p_value"]], method = padjust_method)

  genes_per_category <- interesting[, c("term_id", "intersection")]
  category_genes <- gsub(pattern=",\\s*", replacement="/", x = genes_per_category[["intersection"]])

  ## Right now the cutoff is 1.0, which is not particularly interesting/useful.
  interesting_cutoff_idx <- interesting[["p_value"]] <= cutoff
  interesting_cutoff <- interesting[interesting_cutoff_idx, ]

  ## Note that for the moment I am repeating the pvalue/p.adjust/qvalue because
  ## I am reasonably certain that gprofiler2 does its own adjustment.
  representation_df <- data.frame(
      "ID" = interesting[["term_id"]],
      "Description" = interesting[["term_name"]],
      ## The following two lines are ridiculous, but required for the enrichplots to work.
      "GeneRatio" = paste0(interesting[["intersection_size"]], "/", interesting[["term_size"]]),
      "BgRatio" = paste0(interesting[["term_size"]], "/", interesting[["query_size"]]),
      "pvalue" = interesting[["p_value"]],
      "p.adjust" = interesting[["p_value"]],
      "qvalue" = interesting[["p_value"]],
      "geneID" = category_genes,
      "Count" = interesting[["intersection_size"]],
      stringsAsFactors = FALSE)
  rownames(representation_df) <- representation_df[["ID"]]
  if (is.null(organism)) {
    organism <- "UNKNOWN"
  }
  ret <- new("enrichResult",
             result = representation_df,
             pvalueCutoff = cutoff,
             pAdjustMethod = padjust_method,
             qvalueCutoff = cutoff,
             gene = sig_genes,
             ## universe = extID,
             geneSets = list(up=sig_genes),
             ## geneSets = geneSets,
             organism = organism,
             keytype = "UNKNOWN",
             ontology = ontology,
             readable = FALSE)
  return(ret)
}

## EOF
