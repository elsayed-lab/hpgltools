## Time-stamp: <Fri Apr 15 10:38:25 2016 Ashton Trey Belew (abelew@gmail.com)>

simple_gprofiler <- function(gene_list, organism="hsapiens", direction="descending") {
    ## Assume for the moment a limma-ish data frame
    gene_list <- gene_list[order(-gene_list[["logFC"]]), ]
    gene_ids <- as.vector(gene_list[["ID"]])
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
