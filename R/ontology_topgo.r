#' Perform a simplified topgo analysis.
#'
#' This will attempt to make it easier to run topgo on a set of genes.
#'
#' @param de_genes Data frame of differentially expressed genes, containing IDs any other columns.
#' @param goid_map File containing mappings of genes to goids in the format expected by topgo.
#' @param goids_df Data frame of the goids which may be used to make the goid_map.
#' @param pvals Set of pvalues in the DE data which may be used to improve the topgo results.
#' @param limitby Test to index the results by.
#' @param limit Ontology pvalue to use as the lower limit.
#' @param signodes I don't remember right now.
#' @param sigforall Provide the significance for all nodes?
#' @param numchar Character limit for the table of results.
#' @param selector Function name for choosing genes to include.
#' @param pval_column Column from which to acquire scores.
#' @param overwrite Yeah I do not remember this one either.
#' @param densities Densities, yeah, the densities...
#' @param pval_plots Include pvalue plots of the results a la clusterprofiler?
#' @param parallel  Perform some operations in parallel to speed this up?
#' @param ... Other options which I do not remember right now!
#' @return Big list including the various outputs from topgo
#' @export
simple_topgo <- function(de_genes, goid_map="id2go.map", goids_df=NULL,
                         pvals=NULL, limitby="fisher", limit=0.1, signodes=100,
                         sigforall=TRUE, numchar=300, selector="topDiffGenes",
                         pval_column="adj.P.Val", overwrite=FALSE, densities=FALSE,
                         pval_plots=TRUE, parallel=FALSE, ...) {
### Some neat ideas from the topGO documentation:
### geneList <- getPvalues(exprs(eset), classlabel = y, alternative = "greater")
### A variant of these operations make it possible to give topGO scores so that
### a larger array of tests may be performed
### x <- topDiffGenes(geneList)
### sum(x) ## the number of selected genes
### If we do something like above to give scores to all the 'DEgenes', then we set up the GOdata object like this:
### mf_GOdata = new("topGOdata", description="something", ontology="BP", allGenes = entire_geneList, geneSel=topDiffGenes, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=2)
    ## The following library invocation is in case it was unloaded for pathview
    requireNamespace("topGO")
    require.auto("Hmisc")
    requireNamespace("Hmisc")
    gomap_info <- make_id2gomap(goid_map=goid_map, goids_df=goids_df, overwrite=overwrite)
    geneID2GO <- topGO::readMappings(file=goid_map)
    annotated_genes <- names(geneID2GO)
    if (is.null(de_genes[["ID"]])) {
        de_genes[["ID"]] <- make.names(rownames(de_genes), unique=TRUE)
    }
    fisher_interesting_genes <- as.factor(as.integer(annotated_genes %in% de_genes[["ID"]]))
    names(fisher_interesting_genes) <- annotated_genes
    ks_interesting_genes <- as.integer(!annotated_genes %in% de_genes[["ID"]])
    if (!is.null(de_genes[[pval_column]])) {
        ## I think this needs to include the entire gene universe, not only the set of x differentially expressed genes
        ## Making this an explicit as.vector(as.numeric()) because it turns out the values from DESeq are characters.
        pvals <- as.vector(as.numeric(de_genes[[pval_column]]))
        names(pvals) <- rownames(de_genes)
        for (p in 1:length(pvals)) {
            name <- names(pvals)[p]
            ks_interesting_genes[[name]] <- pvals[p]
        }
    }
    ks_interesting_genes <- as.vector(ks_interesting_genes)
    names(ks_interesting_genes) <- annotated_genes

    tt <- try(sm(requireNamespace("topGO")), silent=TRUE)
    tt <- try(sm(attachNamespace("topGO")), silent=TRUE)
    ## Instead of invoking library(topGO), I can requireNamespace && attachNamespace.
    ## "GOMFTerm not found"
    ## Ergo, requireNamespace() is insufficient!
    godata_fisher_result <- list(
        "MF" = NULL,
        "BP" = NULL,
        "CC" = NULL)
    godata_ks_result <- list(
        "MF" = NULL,
        "BP" = NULL,
        "CC" = NULL)
    sig_fisher_result <- list(
        "MF" = NULL,
        "BP" = NULL,
        "CC" = NULL)
    sig_ks_result <- list(
        "MF" = NULL,
        "BP" = NULL,
        "CC" = NULL)
    sig_weight_result <- list(
        "MF" = NULL,
        "BP" = NULL,
        "CC" = NULL)
    sig_el_result <- list(
        "MF" = NULL,
        "BP" = NULL,
        "CC" = NULL)
    godata_fisher_res <- NULL
    godata_ks_res <- NULL
    sig_fisher_res <- NULL
    sig_ks_res <- NULL
    sig_weight_res <- NULL
    sig_el_res <- NULL
    results <- list()
    if (isTRUE(parallel)) {
        ontologies <- c("MF","BP","CC")
        cl <- parallel::makeCluster(3)  ## 1 for each ontology
        doParallel::registerDoParallel(cl)
        requireNamespace("parallel")
        requireNamespace("doParallel")
        requireNamespace("iterators")
        requireNamespace("foreach")

        ## In each of the following steps, we will make a list of 3 items in an unknown order
        ## But they will have an ontology slot which we can query to figure out if they are
        ## MF/BP/CC.  Ergo, let foreach() %dopar% {}  create the lists in whatever order they want.
        ## Then step through the lists and re-order them appropriately.

        ## Step 1a:  Iterate through the 3 ontology groups and create Fisher-testable data sets
        godata_fisher_res <- foreach(c=1:length(ontologies), .packages=c("hpgltools","topGO")) %dopar% {
            ont <- ontologies[[c]]
            results[[ont]] <- new("topGOdata",
                                  ontology=ont,
                                  allGenes=fisher_interesting_genes,
                                  annot=topGO::annFUN.gene2GO,
                                  gene2GO=geneID2GO)
        } ## End foreach %dopar% for the godata fisher data set
        for (r in 1:length(godata_fisher_res)) {
            type <- godata_fisher_res[[r]]@ontology
            a_result <- godata_fisher_res[[r]]
            godata_fisher_result[[type]] <- a_result
        }
        rm(godata_fisher_res)

        ## Step 1b: As above, but make them suitable for KS tests.
        godata_ks_res <- foreach(c=1:length(ontologies), .packages=c("hpgltools","topGO")) %dopar% {
            ont <- ontologies[[c]]
            results[[ont]] <- try(new("topGOdata",
                                      description=ont,
                                      ontology=ont,
                                      allGenes=ks_interesting_genes,
                                      geneSel=get(selector),
                                      annot=topGO::annFUN.gene2GO,
                                      gene2GO=geneID2GO))
        } ## End the foreach %dopar% for the ks godata set
        for (r in 1:length(godata_ks_res)) {
            type <- godata_ks_res[[r]]@ontology
            a_result <- godata_ks_res[[r]]
            godata_ks_result[[type]] <- a_result
        }
        rm(godata_ks_res)

        ## Step 2:  Perform a fisher test using the fisher-testable data
        test_stat <- new("classicCount", testStatistic=topGO::GOFisherTest, name="Fisher test")
        sig_fisher_res <- foreach(c=1:length(ontologies), .packages=c("hpgltools","topGO")) %dopar% {
            ont <- ontologies[[c]]
            results[[ont]] <- try(topGO::getSigGroups(godata_fisher_result[[ont]], test_stat))
        } ## End the foreach %dopar% to get a significant fisher result
        for (r in 1:length(sig_fisher_res)) {
            description_string <- sig_fisher_res[[r]]@description
            type <- strsplit(x=description_string, split=": ", perl=TRUE)[[1]][[2]]
            a_result <- sig_fisher_res[[r]]
            sig_fisher_result[[type]] <- a_result
        }
        rm(sig_fisher_res)

        ## Step 3:  Perform a KS test using the appropriate data set
        test_stat <- new("classicScore", testStatistic=topGO::GOKSTest, name="KS tests")
        sig_ks_res <- foreach(c=1:length(ontologies), .packages=c("hpgltools","topGO")) %dopar% {
            ont <- ontologies[[c]]
            results[[ont]] <- try(topGO::getSigGroups(godata_ks_result[[ont]], test_stat))
        } ## End the foreach %dopar% to get a significant fisher result
        for (r in 1:length(sig_ks_res)) {
            description_string <- sig_ks_res[[r]]@description
            type <- strsplit(x=description_string, split=": ", perl=TRUE)[[1]][[2]]
            a_result <- sig_ks_res[[r]]
            sig_ks_result[[type]] <- a_result
        }
        rm(sig_ks_res)

        ## Step 4:  Use the KS-testable data to do an elimination score
        test_stat <- new("elimScore", testStatistic=topGO::GOKSTest, name="KS test", cutOff=0.05)
        sig_el_res <- foreach(c=1:length(ontologies), .packages=c("hpgltools","topGO")) %dopar% {
            ont <- ontologies[[c]]
            results[[ont]] <- try(topGO::getSigGroups(godata_ks_result[[ont]], test_stat))
        } ## End the foreach %dopar% to get a significant KS result
        for (r in 1:length(sig_el_res)) {
            description_string <- sig_el_res[[r]]@description
            type <- strsplit(x=description_string, split=": ", perl=TRUE)[[1]][[2]]
            a_result <- sig_el_res[[r]]
            sig_el_result[[type]] <- a_result
        }
        rm(sig_el_res)

        ## Step 5: Finally, use the weighted test on the Fisher-data
        test_stat <- new("weightCount", testStatistic=topGO::GOFisherTest, name="Fisher test", cutOff=0.05)
        sig_weight_res <- foreach(c=1:length(ontologies), .packages=c("hpgltools","topGO")) %dopar% {
            ont <- ontologies[[c]]
            results[[ont]] <- try(topGO::getSigGroups(godata_fisher_result[[ont]], test_stat))
        } ## End the foreach %dopar% to get a significant fisher result
        for (r in 1:length(sig_weight_res)) {
            description_string <- sig_weight_res[[r]]@description
            type <- strsplit(x=description_string, split=": ", perl=TRUE)[[1]][[2]]
            a_result <- sig_weight_res[[r]]
            sig_weight_result[[type]] <- a_result
        }
        rm(sig_weight_res)

        ## We have collected all the data sets, close the 3-cpu cluster.
        parallel::stopCluster(cl)
    } else {
        godata_fisher_result[["MF"]] <- new("topGOdata", ontology="MF",
                                            allGenes=fisher_interesting_genes,
                                            annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        godata_fisher_result[["BP"]] <- new("topGOdata", ontology="BP",
                                            allGenes=fisher_interesting_genes,
                                            annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        godata_fisher_result[["CC"]] <- new("topGOdata", ontology="CC",
                                            allGenes=fisher_interesting_genes,
                                            annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        godata_ks_result[["MF"]] <- new("topGOdata", description="MF", ontology="MF",
                                        allGenes=ks_interesting_genes, geneSel=get(selector),
                                        annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        godata_ks_result[["BP"]] <- new("topGOdata", description="BP", ontology="BP",
                                        allGenes=ks_interesting_genes, geneSel=get(selector),
                                        annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        godata_ks_result[["CC"]] <- new("topGOdata", description="CC", ontology="CC",
                                        allGenes=ks_interesting_genes, geneSel=get(selector),
                                        annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)

        test_stat <- new("classicCount", testStatistic=topGO::GOFisherTest, name="Fisher test")
        sig_fisher_result[["MF"]] <- topGO::getSigGroups(godata_fisher_result[["MF"]], test_stat)
        sig_fisher_result[["BP"]] <- topGO::getSigGroups(godata_fisher_result[["BP"]], test_stat)
        sig_fisher_result[["CC"]] <- topGO::getSigGroups(godata_fisher_result[["CC"]], test_stat)

        test_stat <- new("classicScore", testStatistic=topGO::GOKSTest, name="KS tests")
        sig_ks_result[["MF"]] <- topGO::getSigGroups(godata_ks_result[["MF"]], test_stat)
        sig_ks_result[["BP"]] <- topGO::getSigGroups(godata_ks_result[["BP"]], test_stat)
        sig_ks_result[["CC"]] <- topGO::getSigGroups(godata_ks_result[["CC"]], test_stat)

        ## test_stat <- new("elimScore", testStatistic=topGO::GOKSTest, name="Fisher test", cutOff=0.05)
        test_stat <- new("elimScore", testStatistic=topGO::GOKSTest, name="KS test", cutOff=0.05)
        sig_el_result[["MF"]] <- topGO::getSigGroups(godata_ks_result[["MF"]], test_stat)
        sig_el_result[["BP"]] <- topGO::getSigGroups(godata_ks_result[["BP"]], test_stat)
        sig_el_result[["CC"]] <- topGO::getSigGroups(godata_ks_result[["CC"]], test_stat)
        ## I think the following lines were in error, and they should be using the ks test data.
###        sig_el_result[["mf"]] <- topGO::getSigGroups(fisher_mf_GOdata, test_stat)
###        sig_el_result[["bp"]] <- topGO::getSigGroups(fisher_bp_GOdata, test_stat)
###        sig_el_result[["cc"]] <- topGO::getSigGroups(fisher_cc_GOdata, test_stat)

        test_stat <- new("weightCount", testStatistic=topGO::GOFisherTest, name="Fisher test", sigRatio="ratio")
        sig_weight_result[["MF"]] <- topGO::getSigGroups(godata_fisher_result[["MF"]], test_stat)
        sig_weight_result[["BP"]] <- topGO::getSigGroups(godata_fisher_result[["BP"]], test_stat)
        sig_weight_result[["CC"]] <- topGO::getSigGroups(godata_fisher_result[["CC"]], test_stat)
    }

    mf_fisher_pdist <- try(plot_histogram(sig_fisher_result[["MF"]]@score, bins=20))
    mf_ks_pdist <- try(plot_histogram(sig_ks_result[["MF"]]@score, bins=20))
    mf_el_pdist <- try(plot_histogram(sig_el_result[["MF"]]@score, bins=20))
    mf_weight_pdist <- try(plot_histogram(sig_weight_result[["MF"]]@score, bins=20))
    bp_fisher_pdist <- try(plot_histogram(sig_fisher_result[["BP"]]@score, bins=20))
    bp_ks_pdist <- try(plot_histogram(sig_ks_result[["BP"]]@score, bins=20))
    bp_el_pdist <- try(plot_histogram(sig_el_result[["BP"]]@score, bins=20))
    bp_weight_pdist <- try(plot_histogram(sig_weight_result[["BP"]]@score, bins=20))
    cc_fisher_pdist <- try(plot_histogram(sig_fisher_result[["CC"]]@score, bins=20))
    cc_ks_pdist <- try(plot_histogram(sig_ks_result[["CC"]]@score, bins=20))
    cc_el_pdist <- try(plot_histogram(sig_el_result[["CC"]]@score, bins=20))
    cc_weight_pdist <- try(plot_histogram(sig_weight_result[["CC"]]@score, bins=20))

    p_dists <- list(
        "mf_fisher" = mf_fisher_pdist,
        "bp_fisher" = bp_fisher_pdist,
        "cc_fisher" = cc_fisher_pdist,
        "mf_ks" = mf_ks_pdist,
        "bp_ks" = bp_ks_pdist,
        "cc_ks" = cc_ks_pdist,
        "mf_el" = mf_el_pdist,
        "bp_el" = bp_el_pdist,
        "cc_el" = cc_el_pdist,
        "mf_weight" = mf_weight_pdist,
        "bp_weight" = bp_weight_pdist,
        "cc_weight" = cc_weight_pdist)

    results <- list(
        "fmf_godata" = godata_fisher_result[["MF"]],
        "fbp_godata" = godata_fisher_result[["BP"]],
        "fcc_godata" = godata_fisher_result[["CC"]],
        "kmf_godata" = godata_ks_result[["MF"]],
        "kbp_godata" = godata_ks_result[["BP"]],
        "kcc_godata" = godata_ks_result[["CC"]],
        "mf_fisher" = sig_fisher_result[["MF"]],
        "bp_fisher" = sig_fisher_result[["BP"]],
        "cc_fisher" = sig_fisher_result[["CC"]],
        "mf_ks" = sig_ks_result[["MF"]],
        "bp_ks" = sig_ks_result[["BP"]],
        "cc_ks" = sig_ks_result[["CC"]],
        "mf_el" = sig_el_result[["MF"]],
        "bp_el" = sig_el_result[["BP"]],
        "cc_el" = sig_el_result[["CC"]],
        "mf_weight" = sig_weight_result[["MF"]],
        "bp_weight" = sig_weight_result[["BP"]],
        "cc_weight" = sig_weight_result[["CC"]])

    tables <- try(topgo_tables(results, limitby=limitby, limit=limit), silent=TRUE)
    if (class(tables)[1] == 'try-error') {
        tables <- NULL
    }

    mf_densities <- bp_densities <- cc_densities <- list()
    if (isTRUE(densities)) {
        mf_densities <- suppressMessages(plot_topgo_densities(results[["fmf_godata"]], tables[["mf"]]))
        bp_densities <- suppressMessages(plot_topgo_densities(results[["fbp_godata"]], tables[["bp"]]))
        cc_densities <- suppressMessages(plot_topgo_densities(results[["fcc_godata"]], tables[["cc"]]))
    } else {
        message("simple_topgo(): Set densities=TRUE for ontology density plots.")
    }

    information <- list(
        "results" = results,
        "tables" = tables,
        "mf_densities" = mf_densities,
        "bp_densities" = bp_densities,
        "cc_densities" = cc_densities,
        "pdists" = p_dists)
    pval_plots <- plot_topgo_pval(information)
    information[["pvalue_plots"]] <- pval_plots
    return(information)
}

#' Make pretty tables out of topGO data
#'
#' The topgo function GenTable is neat, but it needs some simplification to not be obnoxious.
#'
#' @param result Topgo result.
#' @param limit Pvalue limit defining 'significant'.
#' @param limitby Type of test to perform.
#' @param numchar How many characters to allow in the description?
#' @param orderby Which of the available columns to order the table by?
#' @param ranksof Which of the available columns are used to rank the data?
#' @export
topgo_tables <- function(result, limit=0.1, limitby="fisher",
                         numchar=300, orderby="classic", ranksof="classic") {
    ## The following if statement could be replaced by get(limitby)
    ## But I am leaving it as a way to ensure that no shenanigans ensue
    mf_allRes <- bp_allRes <- cc_allRes <- mf_interesting <- bp_interesting <- cc_interesting <- NULL
    if (limitby == "fisher") {
        mf_siglist <- names(which(result$mf_fisher@score <= limit))
        bp_siglist <- names(which(result$bp_fisher@score <= limit))
        cc_siglist <- names(which(result$bp_fisher@score <= limit))
    } else if (limitby == "KS") {
        mf_siglist <- names(which(result$mf_ks@score <= limit))
        bp_siglist <- names(which(result$bp_ks@score <= limit))
        cc_siglist <- names(which(result$bp_ks@score <= limit))
    } else if (limitby == "EL") {
        mf_siglist <- names(which(result$mf_el@score <= limit))
        bp_siglist <- names(which(result$bp_el@score <= limit))
        cc_siglist <- names(which(result$bp_el@score <= limit))
    } else if (limitby == "weight") {
        mf_siglist <- names(which(result$mf_weight@score <= limit))
        bp_siglist <- names(which(result$bp_weight@score <= limit))
        cc_siglist <- names(which(result$bp_weight@score <= limit))
    } else {
        stop("I can only limit by: fisher, KS, EL, or weight.")
    }
    mf_topnodes <- length(mf_siglist)
    if (mf_topnodes > 0) {
        mf_allRes <- try(topGO::GenTable(result$fmf_godata, classic=result$mf_fisher, KS=result$mf_ks,
                                         EL=result$mf_el, weight=result$mf_weight, orderBy=orderby,
                                         ranksOf=ranksof, topNodes=mf_topnodes, numChar=numchar))
        if (class(mf_allRes) != 'try-error') {
            mf_qvalues <- as.data.frame(qvalue::qvalue(topGO::score(result$mf_fisher))$qvalues)
            mf_allRes <- merge(mf_allRes, mf_qvalues, by.x="GO.ID", by.y="row.names")
            mf_allRes$classic <- as.numeric(mf_allRes$classic)
            mf_allRes <- mf_allRes[with(mf_allRes, order(classic)), ]
            colnames(mf_allRes) <- c("GO.ID","Term","Annotated","Significant","Expected",
                                     "fisher","KS","EL","weight","qvalue")
            mf_interesting <- subset(mf_allRes, get(limitby) <= limit)
            rownames(mf_interesting) <- NULL
            mf_interesting$ont <- "MF"
            mf_interesting <- mf_interesting[,c("GO.ID","ont","Annotated","Significant","Expected",
                                                "fisher","qvalue","KS","EL","weight","Term")]
        }
    }
    bp_topnodes <- length(bp_siglist)
    if (bp_topnodes > 0) {
        bp_allRes <- try(topGO::GenTable(result$fbp_godata, classic=result$bp_fisher, KS=result$bp_ks,
                                         EL=result$bp_el, weight=result$bp_weight, orderBy=orderby,
                                         ranksOf=ranksof, topNodes=bp_topnodes, numChar=numchar))
        if (class(bp_allRes) != 'try-error') {
            bp_qvalues <- as.data.frame(qvalue::qvalue(topGO::score(result$bp_fisher))$qvalues)
            bp_allRes <- merge(bp_allRes, bp_qvalues, by.x="GO.ID", by.y="row.names", all.x=TRUE)
            bp_allRes$classic <- as.numeric(bp_allRes$classic)
            bp_allRes <- bp_allRes[with(bp_allRes, order(classic)), ]
            colnames(bp_allRes) <- c("GO.ID","Term","Annotated","Significant","Expected",
                                     "fisher","KS","EL","weight","qvalue")
            bp_interesting <- subset(bp_allRes, get(limitby) <= limit)
            rownames(bp_interesting) <- NULL
            bp_interesting$ont <- "BP"
            bp_interesting = bp_interesting[,c("GO.ID","ont","Annotated","Significant","Expected",
                                               "fisher","qvalue","KS","EL","weight","Term")]
        }
    }
    cc_topnodes <- length(cc_siglist)
    if (cc_topnodes > 0) {
        cc_allRes <- try(topGO::GenTable(result$fcc_godata, classic=result$cc_fisher, KS=result$cc_ks,
                                         EL=result$cc_el, weight=result$cc_weight, orderBy=orderby,
                                         ranksOf=ranksof, topNodes=cc_topnodes, numChar=numchar))
        if (class(cc_allRes) != 'try-error') {
            cc_qvalues <- as.data.frame(qvalue::qvalue(topGO::score(result$cc_fisher))$qvalues)
            cc_allRes <- merge(cc_allRes, cc_qvalues, by.x="GO.ID", by.y="row.names")
            cc_allRes$classic <- as.numeric(cc_allRes$classic)
            cc_allRes <- cc_allRes[with(cc_allRes, order(classic)), ]
            colnames(cc_allRes) <- c("GO.ID","Term","Annotated","Significant","Expected",
                                     "fisher","KS","EL","weight","qvalue")
            cc_interesting <- subset(cc_allRes, get(limitby) <= limit)
            rownames(cc_interesting) <- NULL
            cc_interesting$ont <- "CC"
            cc_interesting <- cc_interesting[,c("GO.ID","ont","Annotated","Significant","Expected",
                                                "fisher","qvalue","KS","EL","weight","Term")]
        }
    }
    tables <- list(
        "mf" = mf_allRes,
        "bp" = bp_allRes,
        "cc" = cc_allRes,
        "mf_interesting" = mf_interesting,
        "bp_interesting" = bp_interesting,
        "cc_interesting" = cc_interesting)
    return(tables)
}

#' Make a go mapping from IDs in a format suitable for topGO.
#'
#' When using a non-supported organism, one must write out mappings in the format expected by
#' topgo.  This handles that process and gives a summary of the new table.
#'
#' @param goid_map TopGO mapping file.
#' @param goids_df If there is no goid_map, create it with this data frame.
#' @param overwrite Rewrite the mapping file?
#' @return Summary of the new goid table.
#' @export
make_id2gomap <- function(goid_map="reference/go/id2go.map", goids_df=NULL, overwrite=FALSE) {
    id2go_test <- file.info(goid_map)
    goids_dir <- dirname(goid_map)
    new_go <- NULL
    if (!file.exists(goids_dir)) {
        dir.create(goids_dir, recursive=TRUE)
    }
    if (isTRUE(overwrite)) {
        if (is.null(goids_df)) {
            stop("There is neither a id2go file nor a data frame of goids.")
        } else {
            message("Attempting to generate a id2go file in the format expected by topGO.")

            new_go = reshape2::dcast(goids_df, ID~., value.var="GO",
                            fun.aggregate=paste, collapse = ",")

            ##new_go <- dplyr::ddply(goids_df, plyr::.("ID"), "summarise", GO=paste(unique("GO"), collapse=','))
            write.table(new_go, file=goid_map, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
            rm(id2go_test)
        }
    } else { ## overwrite is not true
        if (is.na(id2go_test[["size"]])) {
            if (is.null(goids_df)) {
                stop("There is neither a id2go file nor a data frame of goids.")
            } else {
                message("Attempting to generate a id2go file in the format expected by topGO.")
                new_go = reshape2::dcast(goids_df, ID~., value.var="GO",
                                         fun.aggregate=paste, collapse = ",")
                ##new_go <- plyr::ddply(goids_df, plyr::.("ID"), "summarise", GO=paste(unique("GO"), collapse=','))
                write.table(new_go, file=goid_map, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
                id2go_test <- file.info(goid_map)
            }
        } else { ## There already exists a file, so return its stats
            ## new_go = id2go_test
            new_go <- goid_map
        }
    }
    return(new_go)
}

hpgl_topdiffgenes <- function(scores, df=get0("de_genes"), direction="up") {
    ## Testing parameters
    ##scores = pvals
    ##df = epi_cl14clbr_high
    ## Here is the original topDiffGenes
    ## topDiffGenes <- function(allScore) {
    ##   return(allScore < 0.01)
    ##}
    ## my version of this will expect a limma result table from which I will extract the entries with low p-values
    ## and logFCs which are high or low
    quartiles <- summary(df)
}

#' A very simple selector of strong scoring genes (by p-value)
#'
#' This function was provided in the topGO documentation, but not defined.
#' It was copied/pasted here.  I have ideas for including up/down expression
#' but have so far deemed them not needed because I am feeding topGO
#' already explicit lists of genes which are up/down/whatever.
#' But it still is likely to be useful to be able to further subset the data.
#'
#' @param allScore The scores of the genes
#' @export
topDiffGenes <- function(allScore) { return(allScore < 0.01) }

#' Plot the ontology DAG.
#'
#' This function was stolen from topgo in order to figure out where it was failing.
#'
#' @param graph Graph from topGO
#' @return Weights!
#' @export
getEdgeWeights <- function(graph) {
    weightsList <- graph::edgeWeights(graph)
    to <- lapply(weightsList, names)
    from <- graph::nodes(graph)

    if (any(is.na(unlist(to))) || any(is.na(from))) {
        stop("Edge names do not match node names.")
    }

    edge.names <- paste(rep(from, Biobase::listLen(to)), unlist(to), sep = "~")
    edge.weights <- unlist(weightsList)
    names(edge.weights) <- edge.names
    return(edge.weights)
}

#' A minor hack of the topGO GOplot function.
#'
#' This allows me to change the line widths from the default.
#'
#' @param dag DAG tree of ontologies.
#' @param sigNodes Set of significant ontologies (with p-values).
#' @param dag.name Name for the graph.
#' @param edgeTypes Types of the edges for graphviz.
#' @param nodeShape.type Shapes on the tree.
#' @param genNodes Generate the nodes?
#' @param wantedNodes Subset of the ontologies to plot.
#' @param showEdges Show the arrows?
#' @param useFullNames Full names of the ontologies (they can get long).
#' @param oldSigNodes I dunno.
#' @param nodeInfo Hmm.
#' @param maxchars Maximum characters per line inside the shapes.
#' @return Topgo plot!
#' @export
hpgl_GOplot <- function(dag, sigNodes, dag.name='GO terms', edgeTypes=TRUE,
                        nodeShape.type=c('box','circle','ellipse','plaintext')[3],
                        genNodes=NULL, wantedNodes=NULL, showEdges=TRUE, useFullNames=TRUE,
                        oldSigNodes=NULL, nodeInfo=NULL, maxchars=30) {
    ## Original function definition had nodeInfo=nodeInfo
    if(!missing(sigNodes)) {
        sigNodeInd = TRUE
    } else {
        sigNodeInd = FALSE
    }

    ## we set the global Graphviz attributes
    ## graphAttrs <- getDefaultAttrs(layoutType = 'dot')
    graphAttrs <- Rgraphviz::getDefaultAttrs(layoutType = 'dot')
    graphAttrs$cluster <- NULL
    graphAttrs$edge$arrowsize = "0.4"
    graphAttrs$edge$weight = "0.01"

    ##graphAttrs$graph$splines <- FALSE
    graphAttrs$graph$size = "12.0,12.0"
    graphAttrs$graph$margin = "0.0,0.0"
    ##  graphAttrs$graph$ranksep = "0.02"
    ##  graphAttrs$graph$nodesep = "0.30"

    ## set the node shape
    graphAttrs$node$shape <- nodeShape.type
    ##graphAttrs$node$fixedsize <- FALSE
    ## set the fontsize for the nodes labels
    graphAttrs$node$fontsize <- '20.0'
    graphAttrs$node$height <- '2.0'
    graphAttrs$node$width <- '3.0'
    graphAttrs$graph$size = "12,12"
    graphAttrs$node$color = "lightblue"
    graphAttrs$node$fontname = "arial"
    graphAttrs$node$style = "invis"

    ## set the local attributes lists
    nodeAttrs <- list()
    edgeAttrs <- list()

    ## try to use adaptive node size
    ##nodeAttrs$fixedsize[nodes(dag)] <- rep(FALSE, numNodes(dag))

    if(is.null(nodeInfo)) {
        nodeInfo <- character(graph::numNodes(dag))
        names(nodeInfo) <- graph::nodes(dag)
    } else {
##        print(class(nodeInfo))
##        nodeInfo <- paste('\\\n', nodeInfo, sep = '')
        nodeInfo = gsub("(\\w.{18}).*(\\\\\\n)","\\1\\2", nodeInfo, perl=TRUE)
        nodeInfo <- paste('\\\n', nodeInfo, sep = '')
    }
    ##teststring = paste("test:", nodeInfo)
    ##print(teststring)

  ## a good idea is to use xxxxxxx instead of GO:xxxxxxx as node labes
    node.names <- graph::nodes(dag)
    if(!useFullNames) {
        nodeAttrs$label <- sapply(node.names,
                                  function(x) {
                                      return(paste(substr(x, 4, nchar(node.names[1])),
                                                   nodeInfo[x], sep = ''))
                                  })
    } else {
        nodeAttrs$label <- paste(node.names, nodeInfo, sep = '')
        names(nodeAttrs$label) <- node.names
    }

  ## we will change the shape and the color of the nodes that generated the dag
    if(!is.null(wantedNodes)) {
        diffNodes <- setdiff(wantedNodes, genNodes)
        if(length(diffNodes) > 0) {
            nodeAttrs$color[diffNodes] <- rep('lightblue', .ln <- length(diffNodes))
            nodeAttrs$shape[diffNodes] <- rep('circle', .ln)
            nodeAttrs$height[diffNodes] <- rep('0.45', .ln)
            ## nodeAttrs$width[diffNodes] <- rep('0.6', .ln)
            ## nodeAttrs$fixedsize[wantedNodes] <- rep(TRUE, .ln)
        }
    }

    ## we will change the shape and the color of the nodes we want back
    if(!is.null(genNodes)) {
        nodeAttrs$color[genNodes] <- rep('lightblue', .ln <- length(genNodes))
        nodeAttrs$shape[genNodes] <- rep('box', .ln)
        ## nodeAttrs$fixedsize[genNodes] <- rep(FALSE, .ln)
    }

    ## we will use different fillcolors for the nodes
    if(sigNodeInd) {
        if(!is.null(oldSigNodes)) {
            old.logSigNodes <- log10(sort(oldSigNodes[graph::nodes(dag)]))
            old.range <- range(old.logSigNodes)
            logSigNodes <- log10(sort(sigNodes))
            logSigNodes[logSigNodes < old.range[1]] <- old.range[1]
            logSigNodes[logSigNodes > old.range[2]] <- old.range[2]
            ## debug:  old.range == range(logSigNodes)
            ## if(!identical(all.equal(old.range, range(logSigNodes)), TRUE)) {
            ## print(old.range)
            ## print(range(logSigNodes))
            ## stop('some stupid error here :)')
            ## }
        } else {
            old.logSigNodes <- logSigNodes <- log10(sort(sigNodes))
        }

        sigColor <- round(logSigNodes - range(logSigNodes)[1] + 1)
        old.sigColor <- round(old.logSigNodes - range(old.logSigNodes)[1] + 1)

        mm <- max(sigColor, old.sigColor)
        sigColor <- sigColor + (mm - max(sigColor))

        colorMap <- heat.colors(mm)
        nodeAttrs$fillcolor <- unlist(lapply(sigColor, function(x) return(colorMap[x])))
    }

    if(!showEdges) {
        graphAttrs$edge$color <- 'white'
    } else {
        ## if we want to differentiate between 'part-of' and 'is-a' edges
        if (edgeTypes) {
            ## 0 for a is_a relation,  1 for a part_of relation
            ## edgeAttrs$color <- ifelse(getEdgeWeights(dag) == 0, 'black', 'red')
            edgeAttrs$color <- ifelse(hpgltools::getEdgeWeights(dag) == 0, 'black', 'black')
        }
        ## plot(dag, attrs = graphAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs)
    }

    final_dag <- Rgraphviz::agopen(graph=dag, name=dag.name, attrs=graphAttrs,
                                   nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs)
    return(final_dag)
}

#' A hack of topGO's groupDensity()
#'
#' This just adds a couple wrappers to avoid errors in groupDensity.
#'
#' @param object TopGO enrichment object.
#' @param whichGO Individual ontology group to compare against.
#' @param ranks Rank order the set of ontologies?
#' @param rm.one Remove pvalue=1 groups?
#' @return plot of group densities.
hpgl_GroupDensity = function(object, whichGO, ranks=TRUE, rm.one=FALSE) {
    groupMembers <- try(topGO::genesInTerm(object, whichGO)[[1]])
    if (class(groupMembers)[1] == 'try-error') {
        return(NULL)
    }
    allS <- topGO::geneScore(object, use.names = TRUE)
    if (rm.one) {
        allS <- allS[allS < 0.99]
    }
    xlab <- "Gene' score"
    if (ranks) {
        allS <- BiocGenerics::rank(allS, ties.method = "random")
        xlab <- "Gene's rank"
    }
    group <- as.integer(names(allS) %in% groupMembers)
    xx <- data.frame(score=allS, group = factor(group,
                                                labels=paste(c("complementary", whichGO),
                                                             "  (", table(group), ")", sep="")))
    plot = lattice::densityplot( ~ score | group, data=xx, layout=c(1,2), xlab=xlab)
    return(plot)
}

## EOF
