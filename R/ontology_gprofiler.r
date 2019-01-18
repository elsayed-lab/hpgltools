#' Run searches against the web service g:Profiler.
#'
#' Thank you Ginger for showing me your thesis, gProfiler is pretty cool!
#'
#' @param sig_genes Guess!  The set of differentially expressed/interesting
#'   genes.
#' @param species  Organism supported by gprofiler.
#' @param convert  Use gProfileR's conversion utility?
#' @param first_col  First place used to define the order of 'significant'.
#' @param second_col  If that fails, try a second column.
#' @param do_go  Perform GO search?
#' @param do_kegg  Perform KEGG search?
#' @param do_reactome  Perform reactome search?
#' @param do_mi  Do miRNA search?
#' @param do_tf  Search for transcription factors?
#' @param do_corum  Do corum search?
#' @param do_hp  Do the hp search?
#' @param significant  Only return the statistically significant hits?
#' @param pseudo_gsea  Is the data in a ranked order by significance?
#' @param id_col  Which column in the table should be used for gene ID
#'   crossreferencing?  gProfiler uses Ensembl ids.  So if you have a table of
#'   entrez or whatever, translate it!
#' @param excel  Print the results to an excel file?
#' @return a list of results for go, kegg, reactome, and a few more.
#' @seealso \pkg{gProfiler}
#' @examples
#' \dontrun{
#'  gprofiler_is_nice_and_easy <- simple_gprofiler(genes, species='mmusculus')
#' }
#' @export
simple_gprofiler <- function(sig_genes, species="hsapiens", convert=TRUE,
                             first_col="logFC", second_col="limma_logfc", do_go=TRUE,
                             do_kegg=TRUE, do_reactome=TRUE, do_mi=TRUE, do_tf=TRUE,
                             do_corum=TRUE, do_hp=TRUE, significant=TRUE,
                             pseudo_gsea=TRUE, id_col="row.names", excel=NULL) {
  ## Assume for the moment a limma-ish data frame
  ## An idea from Dr. Mount: Add the enrichment number of genes as (overlap /
  ## #term) * (total genes / #query) However, the total number is a constant, so
  ## we can likely get the same information from the overlap.size
  gene_list <- NULL
  if (class(sig_genes) == "character") {
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

  do_lst <- list(
    "GO" = do_go,
    "KEGG" = do_kegg,
    "REAC" = do_reactome,
    "MI" = do_mi,
    "TF" = do_tf,
    "CORUM" = do_corum,
    "HP" = do_hp
    )
  result_lst <- list()

  message_string <- "Performing gProfiler %s search of %s genes against %s."
  for (type in names(do_lst)) {
    message(sprintf(message_string, type, length(gene_ids), species))
    Sys.sleep(5)
    ## To avoid the error: "'names' attribute [14] must be the same length as
    ## the vector [1]"
    gene_ids <- as.vector(gene_ids)
    a_result <- try(gProfileR::gprofiler(
                                 query=gene_ids,
                                 organism=species,
                                 significant=significant,
                                 ordered_query=pseudo_gsea,
                                 src_filter=type))
    if (class(a_result) == "try-error") {
      a_result <- data.frame(stringsAsFactors=FALSE)
    }
    a_result[["domain"]] <- as.factor(a_result[["domain"]])
    message(type, " search found ", nrow(a_result), " hits.")
    result_lst[[type]] <- a_result
  }

  retlist <- list(
    "go" = result_lst[["GO"]],
    "kegg" = result_lst[["KEGG"]],
    "reactome" = result_lst[["REAC"]],
    "mi" = result_lst[["MI"]],
    "tf" = result_lst[["TF"]],
    "corum" = result_lst[["CORUM"]],
    "hp" = result_lst[["HP"]],
    "input" = sig_genes)
  retlist[["pvalue_plots"]] <- try(plot_gprofiler_pval(retlist))

  if (!is.null(excel)) {
    message("Writing data to: ", excel, ".")
    excel_ret <- sm(try(write_gprofiler_data(retlist, excel=excel)))
    retlist[["excel"]] <- excel_ret
    message("Finished writing data.")
  }
  class(retlist) <- c("gprofiler_result", "list")
  return(retlist)
}

## EOF
