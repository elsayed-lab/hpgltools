## Time-stamp: <Sat Mar  5 00:50:55 2016 Ashton Trey Belew (abelew@gmail.com)>

#' Enhance the goseq table of gene ontology information.
#'
#' @param df a dataframe of ontology information.  This is intended to
#' be the output from goseq including information like
#' numbers/category, GOids, etc.  It requires a column 'category' which contains: GO:000001 and such.
#' @param file a csv file to which to write the table
#' @return the ontology table with annotation information included
#' @seealso \pkg{goseq}
#' @examples
#' \dontrun{
#'  annotated_go = goseq_table(go_ids)
#'  head(annotated_go, n=1)
#' ## >        category numDEInCat numInCat over_represented_pvalue
#' ## > 571  GO:0006364          9       26            4.655108e-08
#' ## >      under_represented_pvalue       qvalue ontology
#' ## > 571                 1.0000000 6.731286e-05       BP
#' ## >                                term
#' ## > 571                 rRNA processing
#' ## >                               synonym
#' ## > 571        "35S primary transcript processing, GO:0006365"
#' ## >        secondary    definition
#' ## > 571    GO:0006365   Any process involved in the conversion of a primary ribosomal
#' ##          RNA (rRNA) transcript into one or more mature rRNA molecules.
#' }
#' @export
goseq_table <- function(df, file=NULL) {
    if (is.null(df$term)) {
        df$term <- goterm(df$category)
    }
    if (is.null(df$ontology)) {
        df$ontology <- goont(df$category)
    }
    ## df = subset(df, !is.null(term))
    ## Something about this is a disaster FIXME
    ## df = df[ which(!is.null(df$term)), ]
    message("Testing that go categories are defined.")
    df$good <- gotest(df$category)
    message("Removing undefined categories.")
    ## df = subset(df, good == 1)
    df <- df[ which(df$good == 1), ]
    message("Gathering synonyms.")
    df$synonym <- gosyn(df$category)
    ##message("Gathering secondary ids.")
    ##secondary <- try(gosec(df$category), silent=TRUE)
    ##if (class(secondary) != 'try-error') {
    ##    df$secondary <- secondary
    ##}
##    print("Gather approximate go levels.")  ## This function is too slow, commented it out.
##    df$level = golevel(df$categoy)
    message("Gathering category definitions.")
    df$definition <- godef(df$category)
    df <- df[,c("category","numDEInCat","numInCat","over_represented_pvalue",
                "under_represented_pvalue","qvalue","ontology","term",
                "synonym","definition")]
    if (!is.null(file)) {
        write.csv(df, file=file)
    }
    return(df)
}

#' simple_goseq() Perform a simplified goseq analysis
#'
#' @param de_genes a data frame of differentially expressed genes, containing IDs and whatever other columns
#' @param all_genes the universe of possible genes
#' @param lengths the length of each gene with an ID in de_genes
#' @param goids a list of ontology accessions to gene accessions
#' @param doplot   include pwf plots
#' @param adjust   minimum adjusted pvalue
#' @param pvalue   minimum pvalue
#' @param qvalue   minimum qvalue
#' @param goseq_method   testing used by goseq
#' @param padjust_method   which method to adjust the pvalues
#' @param species   optionally choose a species from supportedOrganisms()
#' @param length_db   Source of gene lengths
#' @param gff   gff file source of gene lengths
#' @param ... extra parameters which I do not recall
#' @return a big list including:
#'   the pwd:pwf function,
#'   alldata:the godata dataframe,
#'   pvalue_histogram:p-value histograms,
#'   godata_interesting:the ontology information of the enhanced groups,
#'   term_table:the goterms with some information about them,
#'   mf_subset:a plot of the MF enhanced groups,
#'   mfp_plot:the pvalues of the MF group,
#'   bp_subset:a plot of the BP enhanced groups,
#'   bpp_plot,
#'   cc_subset,
#'   and ccp_plot
#' @seealso \pkg{goseq} \link[goseq]{goseq} \link[goseq]{nullp}
#' @export
simple_goseq <- function(de_genes, all_genes=NULL, lengths=NULL, goids=NULL, doplot=TRUE,
                         adjust=0.1, pvalue=0.1, qvalue=0.1, goseq_method="Wallenius",
                         padjust_method="BH", species=NULL, length_db="ensGene", gff=NULL, ...) {
    message("simple_goseq() makes some pretty hard assumptions about the data it is fed:")
    message("It requires 2 tables, one of GOids which must have columns (gene)ID and GO(category)")
    message("The other table is of gene lengths with columns (gene)ID and (gene)width.")
    message("Other columns are fine, but ignored.")
    if (is.null(lengths) & is.null(gff)) {
        stop("simple_goseq(): Need a length dataframe or gff file for gene lengths.")
    } else if (!is.null(lengths)) {
        message("simple_goseq(): Using the explicit lengths df for gene lengths.")
    } else {
        ## This is probably hopelessly fragile and requires further thought
        length_df <- gff2df(gff)
        lengths <- length_df[length_df$type == 'CDS', ]
        lengths <- length_df[,c("gene_id", "width")]
        colnames(lengths) <- c("ID","width")
    }
    if (is.null(de_genes$ID)) {
        de_genes$ID <- make.names(rownames(de_genes), unique=TRUE)
    }
    if (is.null(de_genes$DE)) {
        de_genes$DE <- 1
    }
    de_vector <- NULL
    de_table <- de_genes[,c("ID","DE")]
    if (is.null(lengths) & is.null(all_genes) & is.null(species)) {
        stop("simple_goseq(): Need either a set of all genes or gene lengths")
    } else if (!is.null(lengths)) {
        message("simple_goseq(): Using the length data to fill in the de vector.")
        de_table <- merge(de_table, lengths, by.x="ID", by.y="ID", all.y=TRUE)
        de_table[is.na(de_table)] <- 0  ## Set the new entries DE status to 0
        rownames(de_table) <- make.names(de_table$ID, unique=TRUE)
        de_vector <- as.vector(de_table$DE)
        names(de_vector) <- rownames(de_table)
    } else if (!is.null(species)) {
        message("simple_goseq(): Using species and length_db to get metadata.")
        gene_names <- as.data.frame(get(paste(species, length_db, "LENGTH", sep = "."))$Gene)
        colnames(gene_names) <- c("ID")
        de_table <- merge(de_table, gene_names, by.x="ID", by.y="ID", all.y=TRUE)
        de_table[is.na(de_table)] <- 0
        rownames(de_table) <- make.names(de_table$ID, unique=TRUE)
        de_vector <- as.vector(de_table$DE)
        names(de_vector) <- rownames(de_table)
    } else { ## If both lengths and all_genes are defined, use all_genes.
        message("simple_goseq(): Using all genes to fill in the de vector.")
        de_table <- merge(de_table, all_genes, by.x="ID", by.y="row.names", all.y=TRUE)
        ##de_table[is.na(de_table)] <- 0  ## Set the new entries DE status to 0
        de_table <- merge(de_table, lengths, by.x="ID", by.y="ID", all.x=TRUE)
        de_table$DE[is.na(de_table$DE)] <- 0  ## Set the new entries DE status to 0
        rownames(de_table) <- make.names(de_table$ID, unique=TRUE)
        de_vector <- as.vector(de_table$DE)
        names(de_vector) <- rownames(de_table)
    }
    pwf <- NULL
    if (is.null(species)) {
        ## length_table = lengths[,c("ID","width")]
        width_vector <- as.vector(de_table$width)
        names(width_vector) <- de_table$ID
        if (is.null(goids)) {
            stop("simple_goseq(): The goids are not defined.")
        }
        goids <- goids[, c("ID","GO")]
        ##colnames(goids) <- c("ID", "GO")
        pwf <- goseq::nullp(DEgenes=de_vector, bias.data=width_vector, plot.fit=doplot)
    } else {
        pwf <- goseq::nullp(de_vector, species, length_db, plot.fit=doplot) ## Taken from the goseq() reference manual
    }
    pwf_plot <- NULL
    if (isTRUE(doplot)) {
        pwf_plot <- recordPlot()
    }
##    godata = goseq(pwf, gene2cat=goids, method='Wallenius')
    godata <- NULL
    if (is.null(species)) {
        godata <- goseq::goseq(pwf, gene2cat=goids, use_genes_without_cat=TRUE, method=goseq_method)
    } else {
        godata <- goseq::goseq(pwf, species, length_db, use_genes_without_cat=TRUE, method=goseq_method)
    }
    goseq_p <- try(hpgltools::hpgl_histogram(godata$over_represented_pvalue, bins=20))
    goseq_p_second <- sort(unique(table(goseq_p$data)), decreasing=TRUE)[2]
    ## Set the y scale to 2x the second highest number
    ## (assuming always that the highest is a p-value of 1)
    goseq_y_limit <- goseq_p_second * 2
    goseq_p <- goseq_p + ggplot2::scale_y_continuous(limits=c(0, goseq_y_limit))
    message("simple_goseq(): Calculating q-values")
    qdata <- godata$over_represented_pvalue
    qdata[qdata > 1] <- 1 ## For scientific numbers which are 1.0000E+00 it might evaluate to 1.0000000000000001
    qdata <- qvalue::qvalue(qdata)
    godata$term <- goterm(godata$category)
    godata$ontology <- goont(godata$category)
    godata <- cbind(godata, qdata$qvalues)
    colnames(godata) <- c("category","over_represented_pvalue","under_represented_pvalue",
                          "numDEInCat","numInCat","term","ontology","qvalue")
    if (is.null(adjust)) {
        godata_interesting <- subset(godata, godata$over_represented_pvalue < pvalue)
        padjust_method <- "none"
    } else {  ## There is a requested pvalue adjustment
        godata_interesting <- subset(godata, p.adjust(godata$over_represented_pvalue, method=padjust_method) <= adjust)
        if (dim(godata_interesting)[1] == 0) {
            message(paste("simple_goseq(): There are no genes with an adj.p<", adjust, " using: ", padjust_method, ".", sep=""))
            message(sprintf("simple_goseq(): Providing genes with raw pvalue<%s", pvalue))
            godata_interesting <- subset(godata, godata$over_represented_pvalue <= pvalue)
            padjust_method <- "none"
        }
    }
    message("simple_goseq(): Filling godata with terms, this is slow.")
    godata_interesting <- goseq_table(godata_interesting)
    message("simple_goseq(): Making pvalue plots for the ontologies.")
    pvalue_plots <- goseq_pval_plots(godata)
    ## mf_subset <- subset(godata, ontology == "MF")
    mf_subset <- godata[godata$ontology == "MF", ]
    ##bp_subset <- subset(godata, ontology == "BP")
    bp_subset <- godata[godata$ontology == "BP", ]
    ## cc_subset <- subset(godata, ontology == "CC")
    cc_subset <- godata[godata$ontology == "CC", ]
    ## mf_interesting <- subset(godata_interesting, ontology == "MF")
    mf_interesting <- godata_interesting[godata_interesting$ontology == "MF", ]
    rownames(mf_interesting) <- mf_interesting$category
    mf_interesting <- mf_interesting[,c("ontology","numDEInCat","numInCat","over_represented_pvalue","qvalue","term")]
    ##bp_interesting <- subset(godata_interesting, ontology == "BP")
    bp_interesting <- godata_interesting[godata_interesting$ontology == "BP", ]
    rownames(bp_interesting) <- bp_interesting$category
    bp_interesting <- bp_interesting[,c("ontology","numDEInCat","numInCat","over_represented_pvalue","qvalue","term")]
    ##cc_interesting <- subset(godata_interesting, ontology == "CC")
    cc_interesting <- godata_interesting[godata_interesting$ontology == "CC", ]
    cc_interesting <- cc_interesting[,c("ontology","numDEInCat","numInCat","over_represented_pvalue","qvalue","term")]
    return_list <- list(input=de_genes, pwf=pwf, pwf_plot=pwf_plot,
                        alldata=godata, pvalue_histogram=goseq_p,
                        godata_interesting=godata_interesting,
                        mf_interesting=mf_interesting, bp_interesting=bp_interesting,
                        cc_interesting=cc_interesting, goadjust_method=goseq_method,
                        adjust_method=padjust_method, mf_subset=mf_subset,
                        mfp_plot=pvalue_plots$mfp_plot, bp_subset=bp_subset,
                        bpp_plot=pvalue_plots$bpp_plot, cc_subset=cc_subset,
                        ccp_plot=pvalue_plots$ccp_plot, qdata=qdata)
    return(return_list)
}

#' Make a pvalue plot from goseq data
#'
#' @param goterms some data from goseq!
#' @param wrapped_width the number of characters before wrapping to help legibility
#' @param cutoff   pvalue cutoff for the plot
#' @param n    how many groups to include
#' @param mincat   minimum size of the category
#' @param level   levels of the ontology tree to use
#' @return plots!
#' @seealso \link[goseq]{goseq} \pkg{clusterProfiler} \code{\link{pval_plot}}
#' @export
goseq_pval_plots <- function(goterms, wrapped_width=20, cutoff=0.1, n=10, mincat=10, level=NULL) {
    ## The following supports stuff like level='level > 3 & level < 6'
        if (!is.null(level)) {
        keepers <- data.frame()
        message("Getting all go levels.  This takes a moment.")
        mf_go <- golevel_df(ont="MF")
        bp_go <- golevel_df(ont="BP")
        cc_go <- golevel_df(ont="CC")
        message("Finished getting go levels.")
        if (class(level) == 'numeric') {
            stmt <- paste0("subset(mf_go, level == ", level, ")")
            mf_go <- eval(parse(text=stmt))
            stmt <- paste0("subset(bp_go, level == ", level, ")")
            bp_go <- eval(parse(text=stmt))
            stmt <- paste0("subset(cc_go, level == ", level, ")")
            cc_go <- eval(parse(text=stmt))
        } else {
            stmt <- paste0("subset(mf_go, ", level, ")")
            mf_go <- eval(parse(text=stmt))
            stmt <- paste0("subset(bp_go, ", level, ")")
            bp_go <- eval(parse(text=stmt))
            stmt <- paste0("subset(cc_go, ", level, ")")
            cc_go <- eval(parse(text=stmt))
        }
        keepers <- rbind(keepers, mf_go)
        keepers <- rbind(keepers, bp_go)
        keepers <- rbind(keepers, cc_go)
        message("Extracting the goterms in your chosen level.")
        goterms <- merge(goterms, keepers, by.x="category", by.y="GO")
    }
    ## TODO: Replace the subset calls with the less noxious which calls.
    plotting_mf <- subset(goterms, complete.cases(goterms))
    plotting_mf$score <- plotting_mf$numDEInCat / plotting_mf$numInCat
    ## plotting_mf <- subset(plotting_mf, ontology == "MF")
    plotting_mf <- plotting_mf[ plotting_mf$ontology == "MF", ]
    ## plotting_mf <- subset(plotting_mf, term != "NULL")
    plotting_mf <- plotting_mf[ plotting_mf$term != "NULL", ]
    ## plotting_mf <- subset(plotting_mf, over_represented_pvalue <= cutoff)
    plotting_mf <- plotting_mf[ plotting_mf$over_represented_pvalue <= cutoff, ]
    ## plotting_mf <- subset(plotting_mf, numInCat > mincat)
    plotting_mf <- plotting_mf[ plotting_mf$numInCat > mincat, ]
    plotting_mf <- plotting_mf[order(plotting_mf$over_represented_pvalue),]
    plotting_mf <- head(plotting_mf, n=n)
    plotting_mf <- plotting_mf[,c("term","over_represented_pvalue","score")]
    plotting_mf$term <- as.character(lapply(strwrap(plotting_mf$term, wrapped_width, simplify=FALSE), paste, collapse="\n"))
    colnames(plotting_mf) <- c("term","pvalue","score")
    mf_pval_plot <- pval_plot(plotting_mf, ontology="MF")

    plotting_bp <- subset(goterms, complete.cases(goterms))
    plotting_bp$score <- plotting_bp$numDEInCat / plotting_bp$numInCat
    ## plotting_bp <- subset(plotting_bp, ontology == "BP")
    plotting_bp <- plotting_bp[ plotting_bp$ontology == "BP", ]
    ## plotting_bp <- subset(plotting_bp, term != "NULL")
    plotting_bp <- plotting_bp[ plotting_bp$term != "NULL", ]
    ## plotting_bp <- subset(plotting_bp, over_represented_pvalue <= cutoff)
    plotting_bp <- plotting_bp[ plotting_bp$over_represented_pvalue <= cutoff, ]
    ## plotting_bp <- subset(plotting_bp, numInCat > mincat)
    plotting_bp <- plotting_bp[ plotting_bp$numInCat > mincat, ]
    plotting_bp <- plotting_bp[order(plotting_bp$over_represented_pvalue),]
    plotting_bp <- head(plotting_bp, n=n)
    plotting_bp <- plotting_bp[,c("term","over_represented_pvalue","score")]
    colnames(plotting_bp) <- c("term","pvalue","score")
    plotting_bp$term <- as.character(lapply(strwrap(plotting_bp$term, wrapped_width, simplify=FALSE), paste, collapse="\n"))
    bp_pval_plot <- pval_plot(plotting_bp, ontology="BP")

    plotting_cc <- subset(goterms, complete.cases(goterms))
    plotting_cc$score <- plotting_cc$numDEInCat / plotting_cc$numInCat
    ## plotting_cc <- subset(plotting_cc, ontology == "CC")
    plotting_cc <- plotting_cc[ plotting_cc$ontology == "CC", ]
    ## plotting_cc <- subset(plotting_cc, term != "NULL")
    plotting_cc <- plotting_cc[ plotting_cc$term != "NULL", ]
    ## plotting_cc <- subset(plotting_cc, over_represented_pvalue <= cutoff)
    plotting_cc <- plotting_cc[ plotting_cc$over_represented_pvalue <= cutoff, ]
    ## plotting_cc <- subset(plotting_cc, numInCat > mincat)
    plotting_cc <- plotting_cc[ plotting_cc$numInCat > mincat, ]
    plotting_cc <- plotting_cc[order(plotting_cc$over_represented_pvalue),]
    plotting_cc <- head(plotting_cc, n=n)
    plotting_cc <- plotting_cc[,c("term","over_represented_pvalue","score")]
    colnames(plotting_cc) <- c("term","pvalue","score")
    plotting_cc$term <- as.character(lapply(strwrap(plotting_cc$term, wrapped_width, simplify=FALSE), paste, collapse="\n"))
    cc_pval_plot <- pval_plot(plotting_cc, ontology="CC")

    pval_plots <- list(mfp_plot=mf_pval_plot, bpp_plot=bp_pval_plot, ccp_plot=cc_pval_plot,
                       mf_subset=plotting_mf, bp_subset=plotting_bp, cc_subset=plotting_cc)
    return(pval_plots)
}

#' Make fun trees a la topgo from goseq data.
#'
#' @param de_genes some differentially expressed genes
#' @param godata data from goseq
#' @param goid_map   file to save go id mapping
#' @param score_limit   score limit for the coloring
#' @param goids_df   a mapping of IDs to GO in the Ramigo expected format
#' @param overwrite   overwrite the trees
#' @param selector   a function for choosing genes
#' @param pval_column  column to acquire pvalues
#' @return a plot!
#' @seealso \pkg{Ramigo}
#' @export
goseq_trees <- function(de_genes, godata, goid_map="reference/go/id2go.map",
                        score_limit=0.01, goids_df=NULL, overwrite=FALSE,
                        selector="topDiffGenes", pval_column="adj.P.Val") {
    make_id2gomap(goid_map=goid_map, goids_df=goids_df, overwrite=overwrite)
    geneID2GO <- topGO::readMappings(file=goid_map)
    annotated_genes <- names(geneID2GO)
    if (is.null(de_genes$ID)) {
        de_genes$ID <- make.names(rownames(de_genes), unique=TRUE)
    }
    interesting_genes <- factor(annotated_genes %in% de_genes$ID)
    names(interesting_genes) <- annotated_genes

    if (is.null(de_genes[[pval_column]])) {
        mf_GOdata <- new("topGOdata", ontology="MF", allGenes=interesting_genes, annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata <- new("topGOdata", ontology="BP", allGenes=interesting_genes, annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata <- new("topGOdata", ontology="CC", allGenes=interesting_genes, annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
    } else {
        pvals <- as.vector(as.numeric(de_genes[[pval_column]]))
        names(pvals) <- rownames(de_genes)
        mf_GOdata <- new("topGOdata", description="MF", ontology="MF", allGenes=pvals,
                         geneSel=get(selector), annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata <- new("topGOdata", description="BP", ontology="BP", allGenes=pvals,
                         geneSel=get(selector), annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata <- new("topGOdata", description="CC", ontology="CC", allGenes=pvals,
                         geneSel=get(selector), annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
    }

    enriched_ids <- godata$alldata$category
    enriched_scores <- godata$alldata$over_represented_pvalue
    names(enriched_scores) <- enriched_ids

    mf_avail_nodes <- as.list(mf_GOdata@graph@nodes)
    names(mf_avail_nodes) <- mf_GOdata@graph@nodes
    mf_nodes <- enriched_scores[names(enriched_scores) %in% names(mf_avail_nodes)]
    mf_included <- length(which(mf_nodes <= score_limit))
    mf_tree_data <- try(suppressWarnings(topGO::showSigOfNodes(mf_GOdata, mf_nodes, useInfo="all",
                                                               sigForAll=TRUE, firstSigNodes=mf_included,
                                                               useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(mf_tree_data) == 'try-error') {
        message("There was an error generating the MF tree.")
        mf_tree <- NULL
    } else {
        mf_tree <- recordPlot()
    }

    bp_avail_nodes <- as.list(bp_GOdata@graph@nodes)
    names(bp_avail_nodes) <- bp_GOdata@graph@nodes
    bp_nodes <- enriched_scores[names(enriched_scores) %in% names(bp_avail_nodes)]
    bp_included <- length(which(bp_nodes <= score_limit))
    bp_tree_data <- try(suppressWarnings(topGO::showSigOfNodes(bp_GOdata, bp_nodes, useInfo="all",
                                                               sigForAll=TRUE, firstSigNodes=bp_included,
                                                               useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(bp_tree_data) == 'try-error') {
        message("There was an error generating the BP tree.")
        bp_tree <- NULL
    } else {
        bp_tree <- recordPlot()
    }

    cc_avail_nodes <- as.list(cc_GOdata@graph@nodes)
    names(cc_avail_nodes) <- cc_GOdata@graph@nodes
    cc_nodes <- enriched_scores[names(enriched_scores) %in% names(cc_avail_nodes)]
    cc_included <- length(which(cc_nodes <= score_limit))
    cc_tree_data <- try(suppressWarnings(topGO::showSigOfNodes(cc_GOdata, cc_nodes, useInfo="all",
                                                               sigForAll=TRUE, firstSigNodes=cc_included,
                                                               useFullNames=TRUE, plotFunction=hpgl_GOplot)))
    if (class(cc_tree_data) == 'try-error') {
        message("There was an error generating the CC tree.")
        cc_tree <- NULL
    } else {
        cc_tree <- recordPlot()
    }
    trees <- list(MF=mf_tree, BP=bp_tree, CC=cc_tree,
                  MFdata=mf_tree_data, BPdata=bp_tree_data, CCdata=cc_tree_data)
    return(trees)
}

# EOF
