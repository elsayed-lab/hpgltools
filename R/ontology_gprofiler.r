#' Run searches against the web service g:Profiler.
#'
#' Thank you Ginger for showing me your thesis, gProfiler is pretty cool!
#'
#' @param sig_genes Guess!  The set of differentially expressed/interesting genes.
#' @param species  Organism supported by gprofiler.
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
#' @return a list of results for go, kegg, reactome, and a few more.
#' @export
simple_gprofiler <- function(sig_genes, species="hsapiens", first_col="logFC",
                             second_col="limma_logfc", do_go=TRUE, do_kegg=TRUE,
                             do_reactome=TRUE, do_mi=TRUE, do_tf=TRUE,
                             do_corum=TRUE, do_hp=TRUE, significant=TRUE,
                             pseudo_gsea=TRUE) {
    ## Assume for the moment a limma-ish data frame
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
        if (!is.null(gene_list[["ID"]])) {
            gene_ids <- as.vector(gene_list[["ID"]])
        } else if (!is.null(rownames(gene_list))) {
            gene_ids <- rownames(gene_list)
        } else {
            stop("Unable to get the set of gene ids.")
        }
    }

    ## Setting 'ordered_query' to TRUE, so rank these by p-value or FC or something
    go_result <- data.frame()
    if (isTRUE(do_go)) {
        message("Performing g:Profiler GO search.")
        Sys.sleep(3)
        go_result <- try(gProfileR::gprofiler(query=gene_ids,
                                              organism=species,
                                              significant=significant,
                                              ordered_query=pseudo_gsea,
                                              src_filter="GO"))
        if (class(go_result) == "try-error") {
            go_result <- data.frame()
        }
        message(paste0("GO search found ", nrow(go_result), " hits."))
    }

    kegg_result <- data.frame()
    if (isTRUE(do_kegg)) {
        message("Performing g:Profiler KEGG search.")
        Sys.sleep(3)
        kegg_result <- try(gProfileR::gprofiler(query=gene_ids,
                                                organism=species,
                                                significant=significant,
                                                ordered_query=pseudo_gsea,
                                                src_filter="KEGG"))
        if (class(kegg_result) == "try-error") {
            kegg_result <- data.frame()
        }
        message(paste0("KEGG search found ", nrow(kegg_result), " hits."))
    }

    reactome_result <- data.frame()
    if (isTRUE(do_reactome)) {
        message("Performing g:Profiler reactome.db search.")
        Sys.sleep(3)
        reactome_result <- try(gProfileR::gprofiler(query=gene_ids,
                                                    organism=species,
                                                    significant=significant,
                                                    ordered_query=pseudo_gsea,
                                                    src_filter="REAC"))
        if (class(reactome_result) == "try-error") {
            reactome_result <- data.frame()
        }
        message(paste0("Reactome search found ", nrow(reactome_result), " hits."))
    }

    mi_result <- data.frame()
    if (isTRUE(do_mi)) {
        message("Performing g:Profiler miRNA search.")
        Sys.sleep(3)
        mi_result <- try(gProfileR::gprofiler(query=gene_ids,
                                              organism=species,
                                              significant=significant,
                                              ordered_query=pseudo_gsea,
                                              src_filter="MI"))
        if (class(mi_result) == "try-error") {
            mi_result <- data.frame()
        }
        message(paste0("miRNA search found ", nrow(mi_result), " hits."))
    }

    tf_result <- data.frame()
    if (isTRUE(do_tf)) {
        message("Performing g:Profiler transcription-factor search.")
        Sys.sleep(3)
        tf_result <- try(gProfileR::gprofiler(query=gene_ids,
                                              organism=species,
                                              significant=significant,
                                              ordered_query=pseudo_gsea,
                                              src_filter="TF"))
        if (class(tf_result) == "try-error") {
            tf_result <- data.frame()
        }
        message(paste0("transcription-factor search found ", nrow(tf_result), " hits."))
    }

    corum_result <- data.frame()
    if (isTRUE(do_corum)) {
        message("Performing g:Profiler corum search.")
        Sys.sleep(3)
        corum_result <- try(gProfileR::gprofiler(query=gene_ids,
                                                 organism=species,
                                                 significant=significant,
                                                 ordered_query=pseudo_gsea,
                                                 src_filter="CORUM"))
        if (class(corum_result) == "try-error") {
            corum_result <- data.frame()
        }
        message(paste0("corum search found ", nrow(corum_result), " hits."))
    }

    hp_result <- data.frame()
    if (isTRUE(do_hp)) {
        message("Performing g:Profiler hp search.")
        Sys.sleep(3)
        hp_result <- try(gProfileR::gprofiler(query=gene_ids,
                                              organism=species,
                                              significant=significant,
                                              ordered_query=pseudo_gsea,
                                              src_filter="HP"))
        if (class(hp_result) == "try-error") {
            hp_result <- data.frame()
        }
        message(paste0("hp search found ", nrow(hp_result), " hits."))
    }

    retlist <- list(
        "go" = go_result,
        "kegg" = kegg_result,
        "reactome" = reactome_result,
        "mi" = mi_result,
        "tf" = tf_result,
        "corum" = corum_result,
        "hp" = hp_result)
    retlist[["plots"]] <- try(plot_gprofiler_pval(retlist))
    return(retlist)
}

## EOF
