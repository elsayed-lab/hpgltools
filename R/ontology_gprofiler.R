#' Run simple_gprofiler on every table from extract_significant_genes()
#'
#' @param sig
all_gprofiler <- function(sig, according_to = "deseq", together = FALSE, ...) {
  ret <- list()
  input_up <- sig[[according_to]][["ups"]]
  input_down <- sig[[according_to]][["downs"]]
  sig_names <- names(input_up)
  for (i in seq_along(sig_names)) {
    Sys.sleep(1)
    name <- sig_names[i]
    retname_up <- paste0(name, "_up")
    retname_down <- paste0(name, "_down")
    up <- input_up[[name]]
    down <- input_down[[name]]
    if (isTRUE(together)) {
      up <- rbind(up, down)
      down <- data.frame()
    }
    if (nrow(up) > 0) {
      ret[[retname_up]] <- simple_gprofiler(up, ...)
    } else {
      ret[[retname_up]] <- NULL
    }
    if (nrow(down) > 0) {
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
                              first_col = "logFC", second_col = "deseq_logfc", do_go = TRUE,
                              do_kegg = TRUE, do_reactome = TRUE, do_mi = TRUE, do_tf = TRUE,
                              do_corum = TRUE, do_hp = TRUE, do_hpa = TRUE, do_wp = TRUE,
                              significant = FALSE, exclude_iea = FALSE, do_under = FALSE,
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

  ## Check that the provided species is in gprofiler2
  species_info <- try(gprofiler2::get_version_info(organism = species))
  if ("try-error" %in% class(species_info)) {
    stop("The organism ", species, " is not supported by gprofiler2.")
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
  interactive_plots <- list()
  gost_plots <- list()
  gost_links <- list()
  for (t in seq_along(type_names)) {
    type <- type_names[t]
    message("Performing gProfiler ", type, " search of ",
            length(gene_ids), " against ", species, ".")
    Sys.sleep(3)
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
                                    sources = type), silent = TRUE)

    a_df <- data.frame(stringsAsFactors = FALSE)
    if ("try-error" %in% class(a_result)) {
      message("The ", type, " method failed for this organism.")
    } else {
      a_df <- a_result[["result"]]
      gost_links[[type]] <- gprofiler2::gost(
                                            query = gene_ids,
                                            organism = species,
                                            evcodes = evcodes,
                                            significant = significant,
                                            ordered_query = pseudo_gsea,
                                            user_threshold = threshold,
                                            correction_method = adjp,
                                            domain_scope = domain_scope,
                                            custom_bg = bg,
                                            sources = type,
                                            as_short_link = TRUE)
      interactive_plots[[type]] <- gprofiler2::gostplot(a_result, capped = TRUE, interactive = TRUE)
      gost_plots[[type]] <- gprofiler2::gostplot(a_result, capped = FALSE, interactive = FALSE)
    }
    message(type, " search found ", nrow(a_result), " hits.")
    retlst[[type]] <- a_df
  } ## End iterating over the set of default sources.

  retlst[["interactive_plots"]] <- interactive_plots
  retlst[["gost_plots"]] <- gost_plots
  retlst[["gost_links"]] <- gost_links
  retlst[["pvalue_plots"]] <- try(plot_gprofiler_pval(retlst))
  retlst[["species_info"]] <- species_info
  if (!is.null(excel)) {
    message("Writing data to: ", excel, ".")
    excel_ret <- sm(try(write_gprofiler_data(retlst, excel = excel)))
    retlst[["excel"]] <- excel_ret
    message("Finished writing data.")
  }
  class(retlst) <- c("gprofiler_result", "list")
  return(retlst)
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
simple_gprofiler <- function(sig_genes, species = "hsapiens", convert = TRUE,
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

## EOF
