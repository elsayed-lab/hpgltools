## Time-stamp: <Wed Apr 27 17:24:12 2016 Ashton Trey Belew (abelew@gmail.com)>

#' Run searches against the web service g:Profiler
#'
#' @param gene_list guess!
#' @param organism an organism supported by gprofiler
#' @param first_col where to search for the order of 'significant' first
#' @param second_col if that fails, try some where else.
#' @return a list of results for go, kegg, reactome, and a few more.
#' @export
simple_gprofiler <- function(gene_list, organism="hsapiens", first_col="logFC", second_col="limma_logfc") {
    ## Assume for the moment a limma-ish data frame
    if (!is.null(gene_list[[first_col]])) {
        gene_list <- gene_list[order(-gene_list[[first_col]]), ]
    } else if (!is.null(gene_list[[second_col]])) {
        gene_list <- gene_list[order(-gene_list[[second_col]]), ]
    }

    gene_ids <- NULL
    if (!is.null(gene_list[["ID"]])) {
        gene_ids <- as.vector(gene_list[["ID"]])
    } else if (!is.null(rownames(gene_list))) {
        gene_ids <- rownames(gene_list)
    } else {
        stop("Unable to get the set of gene ids.")
    }
    ## Setting 'ordered_query' to TRUE, so rank these by p-value or FC or something
    message("Performing g:Profiler GO search.")
    go_result <- try(gProfileR::gprofiler(query=gene_ids, organism=organism, significant=TRUE,
                                          ordered_query=TRUE, src_filter="GO"))
    message("Performing g:Profiler KEGG search.")
    kegg_result <-  try(gProfileR::gprofiler(query=gene_ids, organism=organism, significant=TRUE,
                                             ordered_query=TRUE, src_filter="KEGG"))
    message("Performing g:Profiler reactome.db search.")
    reactome_result <- try(gProfileR::gprofiler(query=gene_ids, organism=organism, significant=TRUE,
                                                ordered_query=TRUE, src_filter="REAC"))
    message("Performing g:Profiler mi search.")
    mi_result <- try(gProfileR::gprofiler(query=gene_ids, organism=organism, significant=TRUE,
                                          ordered_query=TRUE, src_filter="MI"))
    message("Performing g:Profiler tf search.")
    tf_result <- try(gProfileR::gprofiler(query=gene_ids, organism=organism, significant=TRUE,
                                          ordered_query=TRUE, src_filter="TF"))
    message("Performing g:Profiler corum search.")
    corum_result <- try(gProfileR::gprofiler(query=gene_ids, organism=organism, significant=TRUE,
                                             ordered_query=TRUE, src_filter="CORUM"))
    message("Performing g:Profiler hp search.")
    hp_result <- try(gProfileR::gprofiler(query=gene_ids, organism=organism, significant=TRUE,
                                          ordered_query=TRUE, src_filter="HP"))
    retlist <- list(
        "go" = go_result,
        "kegg" = kegg_result,
        "reactome" = reactome_result,
        "mi" = mi_result,
        "tf" = tf_result,
        "corum" = corum_result,
        "hp" = hp_result)
    return(retlist)
}
