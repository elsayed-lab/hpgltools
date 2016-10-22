#' Make a pvalue plot from gprofiler data.
#'
#' The p-value plots from clusterProfiler are pretty, this sets the gprofiler data into a format
#' suitable for plotting in that fashion and returns the resulting plots of significant ontologies.
#'
#' @param gp_result Some data from gProfiler.
#' @param wrapped_width  Maximum width of the text names.
#' @param cutoff P-value cutoff for the plots.
#' @param n Maximum number of ontologies to include.
#' @param group_minsize Minimum ontology group size to include.
#' @param ... Options I might pass from other functions are dropped into arglist.
#' @return List of MF/BP/CC pvalue plots.
#' @seealso \pkg{topgo} \code{clusterProfiler}
#' @export
plot_gprofiler_pval <- function(gp_result, wrapped_width=30, cutoff=0.1, n=30, group_minsize=5, scorer="recall", ...) {
    go_result <- gp_result[["go"]]
    kegg_result <- gp_result[["kegg"]]
    reactome_result <- gp_result[["reactome"]]
    mi_result <- gp_result[["mi"]]
    tf_result <- gp_result[["tf"]]
    corum_result <- gp_result[["corum"]]
    hp_result <- gp_result[["hp"]]

    kept_columns <- c("p.value", "term.size", "query.size", "overlap.size", "recall", "precision", "term.id", "term.name", "relative.depth")
    old_options <- options(scipen=4)
    mf_over[["p.value"]] <- as.numeric(format(x=mf_over[["p.value"]], digits=3, scientific=TRUE))
    bp_over[["p.value"]] <- as.numeric(format(x=bp_over[["p.value"]], digits=3, scientific=TRUE))
    cc_over[["p.value"]] <- as.numeric(format(x=cc_over[["p.value"]], digits=3, scientific=TRUE))
    mf_over <- go_result[go_result[["domain"]] == "MF", ]
    mf_over <- mf_over[, kept_columns]
    bp_over <- go_result[go_result[["domain"]] == "BP", ]
    bp_over <- bp_over[, kept_columns]
    cc_over <- go_result[go_result[["domain"]] == "CC", ]
    cc_over <- cc_over[, kept_columns]

    plotting_mf_over <- mf_over
    mf_pval_plot_over <- NULL
    if (is.null(mf_over) | nrow(mf_over) == 0) {
        plotting_mf_over <- NULL
    } else {
        ## First set the order of the table to be something most descriptive.
        ## For the moment, we want that to be the score.
        plotting_mf_over[["score"]] <- plotting_mf_over[[scorer]]
        new_order <- order(plotting_mf_over[["score"]], decreasing=FALSE)
        plotting_mf_over <- plotting_mf_over[new_order, ]
        ## Drop anything with no term name
        plotting_mf_over <- plotting_mf_over[plotting_mf_over[["term.name"]] != "NULL", ]
        ## Drop anything outside of our pvalue cutoff
        plotting_mf_over <- plotting_mf_over[plotting_mf_over[["p.value"]] <= cutoff, ]
        ## Drop anything with fewer than x genes in the group
        plotting_mf_over <- plotting_mf_over[plotting_mf_over[["query.size"]] >= group_minsize, ]
        ## Because of the way ggplot wants to order the bars, we need to go from the bottom up, ergo tail here.
        ## This ordering will be maintained in the plot by setting the levels of the factor in plot_ontpval, which should have a note.
        plotting_mf_over <- tail(plotting_mf_over, n=n)
        plotting_mf_over <- plotting_mf_over[, c("term.name", "p.value", "recall")]
        colnames(plotting_mf_over) <- c("term", "pvalue", "score")
        plotting_mf_over[["term"]] <- as.character(lapply(strwrap(plotting_mf_over[["term"]], wrapped_width, simplify=FALSE),
                                                          paste, collapse="\n"))
        mf_pval_plot_over <- try(plot_ontpval(plotting_mf_over, ontology="MF"), silent=TRUE)
    }
    if (class(mf_pval_plot_over)[[1]] == "try-error") {
        mf_pval_plot_over <- NULL
    }

    plotting_bp_over <- bp_over
    bp_pval_plot_over <- NULL
    if (is.null(bp_over) | nrow(bp_over) == 0) {
        plotting_bp_over <- NULL
    } else {
        ## First set the order of the table to be something most descriptive.
        ## For the moment, we want that to be the score.
        plotting_bp_over[["score"]] <- plotting_bp_over[[scorer]]
        new_order <- order(plotting_bp_over[["score"]], decreasing=FALSE)
        plotting_bp_over <- plotting_bp_over[new_order, ]
        ## Drop anything with no term name
        plotting_bp_over <- plotting_bp_over[plotting_bp_over[["term.name"]] != "NULL", ]
        ## Drop anything outside of our pvalue cutoff
        plotting_bp_over <- plotting_bp_over[plotting_bp_over[["p.value"]] <= cutoff, ]
        ## Drop anything with fewer than x genes in the group
        plotting_bp_over <- plotting_bp_over[plotting_bp_over[["query.size"]] >= group_minsize, ]
        ## Because of the way ggplot wants to order the bars, we need to go from the bottom up, ergo tail here.
        ## This ordering will be maintained in the plot by setting the levels of the factor in plot_ontpval, which should have a note.
        plotting_bp_over <- tail(plotting_bp_over, n=n)
        plotting_bp_over <- plotting_bp_over[, c("term.name", "p.value", "recall")]
        colnames(plotting_bp_over) <- c("term", "pvalue", "score")
        plotting_bp_over[["term"]] <- as.character(lapply(strwrap(plotting_bp_over[["term"]], wrapped_width, simplify=FALSE),
                                                          paste, collapse="\n"))
        bp_pval_plot_over <- try(plot_ontpval(plotting_bp_over, ontology="BP"), silent=TRUE)
    }
    if (class(bp_pval_plot_over)[[1]] == "try-error") {
        bp_pval_plot_over <- NULL
    }

    plotting_cc_over <- cc_over
    cc_pval_plot_over <- NULL
    if (is.null(cc_over) | nrow(cc_over) == 0) {
        plotting_cc_over <- NULL
    } else {
        ## First set the order of the table to be something most descriptive.
        ## For the moment, we want that to be the score.
        plotting_cc_over[["score"]] <- plotting_cc_over[[scorer]]
        new_order <- order(plotting_cc_over[["score"]], decreasing=FALSE)
        plotting_cc_over <- plotting_cc_over[new_order, ]
        ## Drop anything with no term name
        plotting_cc_over <- plotting_cc_over[plotting_cc_over[["term.name"]] != "NULL", ]
        ## Drop anything outside of our pvalue cutoff
        plotting_cc_over <- plotting_cc_over[plotting_cc_over[["p.value"]] <= cutoff, ]
        ## Drop anything with fewer than x genes in the group
        plotting_cc_over <- plotting_cc_over[plotting_cc_over[["query.size"]] >= group_minsize, ]
        ## Because of the way ggplot wants to order the bars, we need to go from the bottom up, ergo tail here.
        ## This ordering will be maintained in the plot by setting the levels of the factor in plot_ontpval, which should have a note.
        plotting_cc_over <- tail(plotting_cc_over, n=n)
        plotting_cc_over <- plotting_cc_over[, c("term.name", "p.value", "recall")]
        colnames(plotting_cc_over) <- c("term", "pvalue", "score")
        plotting_cc_over[["term"]] <- as.character(lapply(strwrap(plotting_cc_over[["term"]], wrapped_width, simplify=FALSE),
                                                          paste, collapse="\n"))
        cc_pval_plot_over <- try(plot_ontpval(plotting_cc_over, ontology="CC"), silent=TRUE)
    }
    if (class(cc_pval_plot_over)[[1]] == "try-error") {
        cc_pval_plot_over <- NULL
    }

    plotting_kegg_over <- kegg_result
    kegg_pval_plot_over <- NULL
    if (is.null(kegg_result) | nrow(kegg_result) == 0) {
        plotting_kegg_over <- NULL
    } else {
        ## First set the order of the table to be something most descriptive.
        ## For the moment, we want that to be the score.
        plotting_kegg_over[["score"]] <- plotting_kegg_over[[scorer]]
        new_order <- order(plotting_kegg_over[["score"]], decreasing=FALSE)
        plotting_kegg_over <- plotting_kegg_over[new_order, ]
        ## Drop anything with no term name
        plotting_kegg_over <- plotting_kegg_over[plotting_kegg_over[["term.name"]] != "NULL", ]
        ## Drop anything outside of our pvalue cutoff
        plotting_kegg_over <- plotting_kegg_over[plotting_kegg_over[["p.value"]] <= cutoff, ]
        ## Drop anything with fewer than x genes in the group
        plotting_kegg_over <- plotting_kegg_over[plotting_kegg_over[["query.size"]] >= group_minsize, ]
        ## Because of the way ggplot wants to order the bars, we need to go from the bottom up, ergo tail here.
        ## This ordering will be maintained in the plot by setting the levels of the factor in plot_ontpval, which should have a note.
        plotting_kegg_over <- tail(plotting_kegg_over, n=n)
        plotting_kegg_over <- plotting_kegg_over[, c("term.name", "p.value", "recall")]
        colnames(plotting_kegg_over) <- c("term", "pvalue", "score")
        plotting_kegg_over[["term"]] <- as.character(lapply(strwrap(plotting_kegg_over[["term"]], wrapped_width, simplify=FALSE),
                                                          paste, collapse="\n"))
        kegg_pval_plot_over <- try(plot_ontpval(plotting_kegg_over, ontology="KEGG"), silent=TRUE)
    }
    if (class(kegg_pval_plot_over)[[1]] == "try-error") {
        kegg_pval_plot_over <- NULL
    }

    plotting_reactome_over <- reactome_result
    reactome_pval_plot_over <- NULL
    if (is.null(reactome_result) | nrow(reactome_result) == 0) {
        plotting_reactome_over <- NULL
    } else {
        ## First set the order of the table to be something most descriptive.
        ## For the moment, we want that to be the score.
        plotting_reactome_over[["score"]] <- plotting_reactome_over[[scorer]]
        new_order <- order(plotting_reactome_over[["score"]], decreasing=FALSE)
        plotting_reactome_over <- plotting_reactome_over[new_order, ]
        ## Drop anything with no term name
        plotting_reactome_over <- plotting_reactome_over[plotting_reactome_over[["term.name"]] != "NULL", ]
        ## Drop anything outside of our pvalue cutoff
        plotting_reactome_over <- plotting_reactome_over[plotting_reactome_over[["p.value"]] <= cutoff, ]
        ## Drop anything with fewer than x genes in the group
        plotting_reactome_over <- plotting_reactome_over[plotting_reactome_over[["query.size"]] >= group_minsize, ]
        ## Because of the way ggplot wants to order the bars, we need to go from the bottom up, ergo tail here.
        ## This ordering will be maintained in the plot by setting the levels of the factor in plot_ontpval, which should have a note.
        plotting_reactome_over <- tail(plotting_reactome_over, n=n)
        plotting_reactome_over <- plotting_reactome_over[, c("term.name", "p.value", "recall")]
        colnames(plotting_reactome_over) <- c("term", "pvalue", "score")
        plotting_reactome_over[["term"]] <- as.character(lapply(strwrap(plotting_reactome_over[["term"]], wrapped_width, simplify=FALSE),
                                                          paste, collapse="\n"))
        reactome_pval_plot_over <- try(plot_ontpval(plotting_reactome_over, ontology="Reactome"), silent=TRUE)
    }
    if (class(reactome_pval_plot_over)[[1]] == "try-error") {
        reactome_pval_plot_over <- NULL
    }

    plotting_mi_over <- mi_result
    mi_pval_plot_over <- NULL
    if (is.null(mi_result) | nrow(mi_result) == 0) {
        plotting_mi_over <- NULL
    } else {
        ## First set the order of the table to be something most descriptive.
        ## For the moment, we want that to be the score.
        plotting_mi_over[["score"]] <- plotting_mi_over[[scorer]]
        new_order <- order(plotting_mi_over[["score"]], decreasing=FALSE)
        plotting_mi_over <- plotting_mi_over[new_order, ]
        ## Drop anything with no term name
        plotting_mi_over <- plotting_mi_over[plotting_mi_over[["term.name"]] != "NULL", ]
        ## Drop anything outside of our pvalue cutoff
        plotting_mi_over <- plotting_mi_over[plotting_mi_over[["p.value"]] <= cutoff, ]
        ## Drop anything with fewer than x genes in the group
        plotting_mi_over <- plotting_mi_over[plotting_mi_over[["query.size"]] >= group_minsize, ]
        ## Because of the way ggplot wants to order the bars, we need to go from the bottom up, ergo tail here.
        ## This ordering will be maintained in the plot by setting the levels of the factor in plot_ontpval, which should have a note.
        plotting_mi_over <- tail(plotting_mi_over, n=n)
        plotting_mi_over <- plotting_mi_over[, c("term.name", "p.value", "recall")]
        colnames(plotting_mi_over) <- c("term", "pvalue", "score")
        plotting_mi_over[["term"]] <- as.character(lapply(strwrap(plotting_mi_over[["term"]], wrapped_width, simplify=FALSE),
                                                          paste, collapse="\n"))
        mi_pval_plot_over <- try(plot_ontpval(plotting_mi_over, ontology="miRNA"), silent=TRUE)
    }
    if (class(mi_pval_plot_over)[[1]] == "try-error") {
        mi_pval_plot_over <- NULL
    }

    plotting_tf_over <- tf_result
    tf_pval_plot_over <- NULL
    if (is.null(tf_result) | nrow(tf_result) == 0) {
        plotting_tf_over <- NULL
    } else {
        ## First set the order of the table to be something most descriptive.
        ## For the moment, we want that to be the score.
        plotting_tf_over[["score"]] <- plotting_tf_over[[scorer]]
        new_order <- order(plotting_tf_over[["score"]], decreasing=FALSE)
        plotting_tf_over <- plotting_tf_over[new_order, ]
        ## Drop anything with no term name
        plotting_tf_over <- plotting_tf_over[plotting_tf_over[["term.name"]] != "NULL", ]
        ## Drop anything outside of our pvalue cutoff
        plotting_tf_over <- plotting_tf_over[plotting_tf_over[["p.value"]] <= cutoff, ]
        ## Drop anything with fewer than x genes in the group
        plotting_tf_over <- plotting_tf_over[plotting_tf_over[["query.size"]] >= group_minsize, ]
        ## Because of the way ggplot wants to order the bars, we need to go from the bottom up, ergo tail here.
        ## This ordering will be maintained in the plot by setting the levels of the factor in plot_ontpval, which should have a note.
        plotting_tf_over <- tail(plotting_tf_over, n=n)
        plotting_tf_over <- plotting_tf_over[, c("term.name", "p.value", "recall")]
        colnames(plotting_tf_over) <- c("term", "pvalue", "score")
        plotting_tf_over[["term"]] <- as.character(lapply(strwrap(plotting_tf_over[["term"]], wrapped_width, simplify=FALSE),
                                                          paste, collapse="\n"))
        tf_pval_plot_over <- try(plot_ontpval(plotting_tf_over, ontology="Transcription factors"), silent=TRUE)
    }
    if (class(tf_pval_plot_over)[[1]] == "try-error") {
        tf_pval_plot_over <- NULL
    }

    plotting_corum_over <- corum_result
    corum_pval_plot_over <- NULL
    if (is.null(corum_result) | nrow(corum_result) == 0) {
        plotting_corum_over <- NULL
    } else {
        ## First set the order of the table to be something most descriptive.
        ## For the moment, we want that to be the score.
        plotting_corum_over[["score"]] <- plotting_corum_over[[scorer]]
        new_order <- order(plotting_corum_over[["score"]], decreasing=FALSE)
        plotting_corum_over <- plotting_corum_over[new_order, ]
        ## Drop anything with no term name
        plotting_corum_over <- plotting_corum_over[plotting_corum_over[["term.name"]] != "NULL", ]
        ## Drop anything outside of our pvalue cutoff
        plotting_corum_over <- plotting_corum_over[plotting_corum_over[["p.value"]] <= cutoff, ]
        ## Drop anything with fewer than x genes in the group
        plotting_corum_over <- plotting_corum_over[plotting_corum_over[["query.size"]] >= group_minsize, ]
        ## Because of the way ggplot wants to order the bars, we need to go from the bottom up, ergo tail here.
        ## This ordering will be maintained in the plot by setting the levels of the factor in plot_ontpval, which should have a note.
        plotting_corum_over <- tail(plotting_corum_over, n=n)
        plotting_corum_over <- plotting_corum_over[, c("term.name", "p.value", "recall")]
        colnames(plotting_corum_over) <- c("term", "pvalue", "score")
        plotting_corum_over[["term"]] <- as.character(lapply(strwrap(plotting_corum_over[["term"]], wrapped_width, simplify=FALSE),
                                                          paste, collapse="\n"))
        corum_pval_plot_over <- try(plot_ontpval(plotting_corum_over, ontology="Corum"), silent=TRUE)
    }
    if (class(corum_pval_plot_over)[[1]] == "try-error") {
        corum_pval_plot_over <- NULL
    }

    plotting_hp_over <- hp_result
    hp_pval_plot_over <- NULL
    if (is.null(hp_result) | nrow(hp_result) == 0) {
        plotting_hp_over <- NULL
    } else {
        ## First set the order of the table to be something most descriptive.
        ## For the moment, we want that to be the score.
        plotting_hp_over[["score"]] <- plotting_hp_over[[scorer]]
        new_order <- order(plotting_hp_over[["score"]], decreasing=FALSE)
        plotting_hp_over <- plotting_hp_over[new_order, ]
        ## Drop anything with no term name
        plotting_hp_over <- plotting_hp_over[plotting_hp_over[["term.name"]] != "NULL", ]
        ## Drop anything outside of our pvalue cutoff
        plotting_hp_over <- plotting_hp_over[plotting_hp_over[["p.value"]] <= cutoff, ]
        ## Drop anything with fewer than x genes in the group
        plotting_hp_over <- plotting_hp_over[plotting_hp_over[["query.size"]] >= group_minsize, ]
        ## Because of the way ggplot wants to order the bars, we need to go from the bottom up, ergo tail here.
        ## This ordering will be maintained in the plot by setting the levels of the factor in plot_ontpval, which should have a note.
        plotting_hp_over <- tail(plotting_hp_over, n=n)
        plotting_hp_over <- plotting_hp_over[, c("term.name", "p.value", "recall")]
        colnames(plotting_hp_over) <- c("term", "pvalue", "score")
        plotting_hp_over[["term"]] <- as.character(lapply(strwrap(plotting_hp_over[["term"]], wrapped_width, simplify=FALSE),
                                                          paste, collapse="\n"))
        hp_pval_plot_over <- try(plot_ontpval(plotting_hp_over, ontology="Human pathology"), silent=TRUE)
    }
    if (class(hp_pval_plot_over)[[1]] == "try-error") {
        hp_pval_plot_over <- NULL
    }

    pval_plots <- list(
        "mfp_plot_over" = mf_pval_plot_over,
        "bpp_plot_over" = bp_pval_plot_over,
        "ccp_plot_over" = cc_pval_plot_over,
        "kegg_plot_over" = kegg_pval_plot_over,
        "reactome_plot_over" = reactome_pval_plot_over,
        "mi_plot_over" = mi_pval_plot_over,
        "tf_plot_over" = tf_pval_plot_over,
        "corum_plot_over" = corum_pval_plot_over,
        "hp_plot_over" = hp_pval_plot_over,
        "mf_subset_over" = plotting_mf_over,
        "bp_subset_over" = plotting_bp_over,
        "cc_subset_over" = plotting_cc_over,
        "kegg_subset" = plotting_kegg_over,
        "reactome_subset" = plotting_reactome_over,
        "mi_subset" = plotting_mi_over,
        "tf_subset" = plotting_tf_over,
        "corum_subset" = plotting_corum_over,
        "hp_subset" = plotting_hp_over
    )
    new_options <- options(old_options)
    return(pval_plots)
}

#' Run searches against the web service g:Profiler.
#'
#' Thank you Ginger for showing me your thesis, gProfiler is pretty cool!
#'
#' @param de_genes guess!
#' @param species an organism supported by gprofiler
#' @param first_col where to search for the order of 'significant' first
#' @param second_col if that fails, try some where else.
#' @param do_go Perform GO search?
#' @param do_kegg Perform KEGG search?
#' @param do_reactome Perform reactome search?
#' @param do_mi Do miRNA search?
#' @param do_tf Search for transcription factors?
#' @param do_corum Do corum search?
#' @param do_hp Do the hp search?
#' @return a list of results for go, kegg, reactome, and a few more.
#' @export
simple_gprofiler <- function(de_genes, species="hsapiens", first_col="logFC",
                             second_col="limma_logfc", do_go=TRUE, do_kegg=TRUE,
                             do_reactome=TRUE, do_mi=TRUE, do_tf=TRUE,
                             do_corum=TRUE, do_hp=TRUE) {
    ## Assume for the moment a limma-ish data frame
    gene_list <- NULL
    if (!is.null(de_genes[[first_col]])) {
        gene_list <- de_genes[order(-de_genes[[first_col]]), ]
    } else if (!is.null(de_genes[[second_col]])) {
        gene_list <- de_genes[order(-de_genes[[second_col]]), ]
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
    go_result <- data.frame()
    if (isTRUE(do_go)) {
        message("Performing g:Profiler GO search.")
        Sys.sleep(3)
        go_result <- try(gProfileR::gprofiler(query=gene_ids, organism=species, significant=TRUE,
                                              ordered_query=TRUE, src_filter="GO"))
        if (class(go_result) == "try-error") {
            go_result <- data.frame()
        }
        message(paste0("GO search found ", nrow(go_result), " hits."))
    }

    kegg_result <- data.frame()
    if (isTRUE(do_kegg)) {
        message("Performing g:Profiler KEGG search.")
        Sys.sleep(3)
        kegg_result <- try(gProfileR::gprofiler(query=gene_ids, organism=species, significant=TRUE,
                                                ordered_query=TRUE, src_filter="KEGG"))
        if (class(kegg_result) == "try-error") {
            kegg_result <- data.frame()
        }
        message(paste0("KEGG search found ", nrow(kegg_result), " hits."))
    }

    reactome_result <- data.frame()
    if (isTRUE(do_reactome)) {
        message("Performing g:Profiler reactome.db search.")
        Sys.sleep(3)
        reactome_result <- try(gProfileR::gprofiler(query=gene_ids, organism=species, significant=TRUE,
                                                    ordered_query=TRUE, src_filter="REAC"))
        if (class(reactome_result) == "try-error") {
            reactome_result <- data.frame()
        }
        message(paste0("Reactome search found ", nrow(reactome_result), " hits."))
    }

    mi_result <- data.frame()
    if (isTRUE(do_mi)) {
        message("Performing g:Profiler miRNA search.")
        Sys.sleep(3)
        mi_result <- try(gProfileR::gprofiler(query=gene_ids, organism=species, significant=TRUE,
                                              ordered_query=TRUE, src_filter="MI"))
        if (class(mi_result) == "try-error") {
            mi_result <- data.frame()
        }
        message(paste0("miRNA search found ", nrow(mi_result), " hits."))
    }

    tf_result <- data.frame()
    if (isTRUE(do_tf)) {
        message("Performing g:Profiler transcription-factor search.")
        Sys.sleep(3)
        tf_result <- try(gProfileR::gprofiler(query=gene_ids, organism=species, significant=TRUE,
                                              ordered_query=TRUE, src_filter="TF"))
        if (class(tf_result) == "try-error") {
            tf_result <- data.frame()
        }
        message(paste0("transcription-factor search found ", nrow(tf_result), " hits."))
    }

    corum_result <- data.frame()
    if (isTRUE(do_corum)) {
        message("Performing g:Profiler corum search.")
        Sys.sleep(3)
        corum_result <- try(gProfileR::gprofiler(query=gene_ids, organism=species, significant=TRUE,
                                                 ordered_query=TRUE, src_filter="CORUM"))
        if (class(corum_result) == "try-error") {
            corum_result <- data.frame()
        }
        message(paste0("corum search found ", nrow(corum_result), " hits."))
    }

    hp_result <- data.frame()
    if (isTRUE(do_hp)) {
        message("Performing g:Profiler hp search.")
        Sys.sleep(3)
        hp_result <- try(gProfileR::gprofiler(query=gene_ids, organism=species, significant=TRUE,
                                              ordered_query=TRUE, src_filter="HP"))
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

#' Write some excel results from a gprofiler search.
#'
#' Gprofiler is pretty awesome.  This function will attempt to write its results to an excel file.
#'
#' @param gprofiler_result  The result from simple_gprofiler().
#' @param wb  Optional workbook object, if you wish to append to an existing workbook.
#' @param excel  Excel file to which to write.
#' @param add_plots  Add some pvalue plots?
#' @param ... More options, not currently used I think.
#' @export
write_gprofiler_data <- function(gprofiler_result, wb=NULL, excel="excel/gprofiler_result.xlsx",
                                 add_plots=TRUE, height=15, width=10, ...) {
    arglist <- list(...)
    table_style <- "TableStyleMedium9"
    if (!is.null(arglist[["table_style"]])) {
        table_style <- arglist[["TableStyleMedium9"]]
    }

    if (is.null(wb)) {
        wb <- openxlsx::createWorkbook(creator="atb")
        hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold", border="Bottom", fontSize="30")
    }

    bp_data <- NULL
    mf_data <- NULL
    cc_data <- NULL
    kegg_data <- NULL
    tf_data <- NULL
    react_data <- NULL
    mi_data <- NULL
    hp_data <- NULL
    corum_data <- NULL

    if (!is.null(gprofiler_result[["go"]]) & nrow(gprofiler_result[["go"]]) > 0) {
        new_row <- 1
        sheet <- "GO"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        go_data <- gprofiler_result[["go"]]
        bp_data <- go_data[go_data[["domain"]]=="BP", ]
        mf_data <- go_data[go_data[["domain"]]=="MF", ]
        cc_data <- go_data[go_data[["domain"]]=="CC", ]

        openxlsx::writeData(wb, sheet, paste0("BP Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=bp_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["bpp_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(bp_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(bp_data) + 2

        openxlsx::writeData(wb, sheet, paste0("MF Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=mf_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["mfp_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(mf_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(mf_data) + 2

        openxlsx::writeData(wb, sheet, paste0("CC Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=cc_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["ccp_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(cc_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(cc_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking if go data is null

    if (!is.null(gprofiler_result[["kegg"]]) & nrow(gprofiler_result[["kegg"]]) > 0) {
        new_row <- 1
        sheet <- "KEGG"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        kegg_data <- gprofiler_result[["kegg"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=kegg_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["kegg_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(kegg_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(kegg_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking KEGG data


    if (!is.null(gprofiler_result[["tf"]]) & nrow(gprofiler_result[["tf"]]) > 0) {
        new_row <- 1
        sheet <- "transcription_factor"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        tf_data <- gprofiler_result[["tf"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=tf_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["tf_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(tf_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(tf_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking tf data

    if (!is.null(gprofiler_result[["reactome"]]) & nrow(gprofiler_result[["reactome"]]) > 0) {
        new_row <- 1
        sheet <- "reactome"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        react_data <- gprofiler_result[["reactome"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=react_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["reactome_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(react_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(react_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking reactome data


    if (!is.null(gprofiler_result[["mi"]]) & nrow(gprofiler_result[["mi"]]) > 0) {
        new_row <- 1
        sheet <- "mirna"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        mi_data <- gprofiler_result[["mi"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=mi_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["mi_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(mi_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(mi_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking mirna data

    if (!is.null(gprofiler_result[["hp"]]) & nrow(gprofiler_result[["hp"]]) > 0) {
        new_row <- 1
        sheet <- "hp"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        hp_data <- gprofiler_result[["hp"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=hp_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["hp_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(hp_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(hp_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking corum data

    if (!is.null(gprofiler_result[["corum"]]) & nrow(gprofiler_result[["corum"]]) > 0) {
        new_row <- 1
        sheet <- "corum"
        openxlsx::addWorksheet(wb, sheetName=sheet)
        corum_data <- gprofiler_result[["corum"]]
        openxlsx::writeData(wb, sheet, paste0("Results from ", sheet, "."), startRow=new_row)
        openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
        new_row <- new_row + 1
        dfwrite <- try(openxlsx::writeDataTable(wb, sheet, x=corum_data, tableStyle=table_style, startRow=new_row))
        ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
        if (isTRUE(add_plots)) {
            a_plot <- gprofiler_result[["plots"]][["corum_plot_over"]]
            try(print(a_plot), silent=TRUE)
            ins <- try(openxlsx::insertPlot(wb, sheet, width=width, height=height,
                                            startCol=ncol(corum_data) + 2, startRow=new_row, fileType="png", units="in"))
        }
        new_row <- new_row + nrow(corum_data) + 2
        openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    } ## End checking corum data

    excel_ret <- NULL
    if (!is.null(excel)) {
        excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite=TRUE))
    }
    if (class(excel_ret) == "try-error") {
        ## In case of an error with zip(1) in making an excel file.
        csvfile <- paste0(excel, "_gomf.csv")
        write.csv(mf_data, file=csvfile)
        csvfile <- paste0(excel, "_gobp.csv")
        write.csv(bp_data, file=csvfile)
        csvfile <- paste0(excel, "_gocc.csv")
        write.csv(cc_data, file=csvfile)
        csvfile <- paste0(excel, "_kegg.csv")
        write.csv(kegg_data, file=csvfile)
        csvfile <- paste0(excel, "_tf.csv")
        write.csv(tf_data, file=csvfile)
        csvfile <- paste0(excel, "_reactome.csv")
        write.csv(react_data, file=csvfile)
        csvfile <- paste0(excel, "_mi.csv")
        write.csv(mi_data, file=csvfile)
        csvfile <- paste0(excel, "_hp.csv")
        write.csv(hp_data, file=csvfile)
        csvfile <- paste0(excel, "_corum.csv")
        write.csv(corum_data, file=csvfile)
        excel_ret <- "csv"
    }
    return(excel_ret)
}

## EOF
