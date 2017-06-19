#' Enhance the goseq table of gene ontology information.
#'
#' While goseq has some nice functionality, the table of outputs it provides is somewhat lacking.
#' This attempts to increase that with some extra helpful data like ontology categories,
#' definitions, etc.
#'
#' @param df Dataframe of ontology information.  This is intended to be the output from goseq
#'  including information like numbers/category, GOids, etc.  It requires a column 'category'
#'  which contains: GO:000001 and such.
#' @param file Csv file to which to write the table.
#' @return Ontology table with annotation information included.
#' @seealso \pkg{goseq}
#' @examples
#' \dontrun{
#'  annotated_go = goseq_table(go_ids)
#'  head(annotated_go, n=1)
#'  ## >        category numDEInCat numInCat over_represented_pvalue
#'  ## > 571  GO:0006364          9       26            4.655108e-08
#'  ## >      under_represented_pvalue       qvalue ontology
#'  ## > 571                 1.0000000 6.731286e-05       BP
#'  ## >                                term
#'  ## > 571                 rRNA processing
#'  ## >                               synonym
#'  ## > 571        "35S primary transcript processing, GO:0006365"
#'  ## >        secondary    definition
#'  ## > 571    GO:0006365   Any process involved in the conversion of a primary ribosomal
#'  ##          RNA (rRNA) transcript into one or more mature rRNA molecules.
#' }
#' @export
goseq_table <- function(df, file=NULL) {
    if (is.null(df[["term"]])) {
        df[["term"]] <- goterm(df[["category"]])
    }
    if (is.null(df[["ontology"]])) {
        df[["ontology"]] <- goont(df[["category"]])
    }
    ## df = subset(df, !is.null(term))
    ## Something about this is a disaster FIXME
    ## df = df[ which(!is.null(df$term)), ]
    message("Testing that go categories are defined.")
    df[["good"]] <- gotest(df[["category"]])
    message("Removing undefined categories.")
    ## df = subset(df, good == 1)
    df <- df[ which(df[["good"]] == 1), ]
    message("Gathering synonyms.")
    df[["synonym"]] <- gosyn(df[["category"]])
    ##message("Gathering secondary ids.")
    ##secondary <- try(gosec(df$category), silent=TRUE)
    ##if (class(secondary) != 'try-error') {
    ##    df$secondary <- secondary
    ##}
    ## print("Gather approximate go levels.")  ## This function is too slow, commented it out.
    ## df[["level"]] = golevel(df[["categoy"]])
    message("Gathering category definitions.")
    df[["definition"]] <- godef(df[["category"]])
    df <- df[, c("category", "numDEInCat", "numInCat", "over_represented_pvalue",
                 "under_represented_pvalue", "qvalue", "ontology", "term",
                 "synonym", "definition")]
    if (!is.null(file)) {
        write.csv(df, file=file)
    }
    return(df)
}

#' Perform a simplified goseq analysis.
#'
#' goseq can be pretty difficult to get set up for non-supported organisms.  This attempts to make
#' that process a bit simpler as well as give some standard outputs which should be similar to those
#' returned by clusterprofiler/topgo/gostats/gprofiler.
#'
#' @param sig_genes Data frame of differentially expressed genes, containing IDs etc.
#' @param go_db Database of go to gene mappings (OrgDb/OrganismDb)
#' @param length_db Database of gene lengths (gff/TxDb)
#' @param doplot Include pwf plots?
#' @param adjust Minimum adjusted pvalue for 'significant.'
#' @param pvalue Minimum pvalue for 'significant.'
#' @param qvalue Minimum qvalue for 'significant.'
#' @param length_keytype Keytype to provide to extract lengths
#' @param go_keytype Keytype to provide to extract go IDs
#' @param goseq_method Statistical test for goseq to use.
#' @param padjust_method Which method to use to adjust the pvalues.
#' @param bioc_length_db Source of gene lengths?
#' @param ... Extra parameters which I do not recall
#' @return Big list including:
#'  the pwd:pwf function,
#'  alldata:the godata dataframe,
#'  pvalue_histogram:p-value histograms,
#'  godata_interesting:the ontology information of the enhanced groups,
#'  term_table:the goterms with some information about them,
#'  mf_subset:a plot of the MF enhanced groups,
#'  mfp_plot:the pvalues of the MF group,
#'  bp_subset:a plot of the BP enhanced groups,
#'  bpp_plot,
#'  cc_subset,
#'  and ccp_plot
#' @seealso \pkg{goseq} \pkg{GO.db}
#' @examples
#' \dontrun{
#'  lotsotables <- simple_goseq(gene_list, godb, lengthdb)
#' }
#' @export
simple_goseq <- function(sig_genes, go_db=NULL, length_db=NULL, doplot=TRUE,
                         adjust=0.1, pvalue=0.1, qvalue=0.1,
                         length_keytype="transcripts", go_keytype="ENTREZID",
                         goseq_method="Wallenius", padjust_method="BH",
                         bioc_length_db="ensGene",
                         ...) {
    arglist <- list(...)

    minimum_interesting <- 1
    if (!is.null(arglist[["minimum_interesting"]])) {
        minimum_interesting <- arglist[["minimum_interesting"]]
    }
    length_df <- data.frame()
    length_vector <- vector()
    de_vector <- vector()
    gene_ids <- NULL
    final_keytype <- NULL
    goids_df <- NULL
    ## sig_genes may be a list, character list, or data frame.
    gene_list <- NULL
    if (class(sig_genes) == "character") {
        ## Then this is a character list of gene ids
        gene_list <- sig_genes
    } else if (class(sig_genes) == "list") {
        gene_list <- names(sig_genes)
    } else if (class(sig_genes) == "data.frame") {
        if (is.null(rownames(sig_genes)) & is.null(sig_genes[["ID"]])) {
            stop("This requires a set of gene IDs either from the rownames or a column named 'ID'.")
        } else if (!is.null(sig_genes[["ID"]])) {
            ## Use a column named 'ID' first because a bunch of annotation databases use ENTREZ
            ## IDs which are just integers, which of course is not allowed by data frame row names.
            message("Using the ID column from your table rather than the row names.")
            gene_list <- sig_genes[["ID"]]
        } else if (!is.null(rownames(sig_genes))) {
            message("Using the row names of your table.")
            gene_list <- rownames(sig_genes)
        } else {
            gene_list <- sig_genes[["ID"]]
        }
    } else {
        stop("Not sure how to handle your set of significant gene ids.")
    }
    ## At this point I should have a character list of gene ids named 'gene_list'
    de_genelist <- as.data.frame(gene_list)
    de_genelist[["DE"]] <- 1
    colnames(de_genelist) <- c("ID", "DE")

    ## Database of lengths may be a gff file, TxDb, or OrganismDb
    metadf <- NULL
    if (class(length_db)[[1]] == "character")  {
        ## Then this should be either a gff file or species name.
        if (grepl(pattern="\\.gff", x=length_db, perl=TRUE) |
            grepl(pattern="\\.gtf", x=length_db, perl=TRUE)) {
            ## gff file
            txdb <- GenomicFeatures::makeTxDbFromGFF(length_db)
            metadf <- extract_lengths(db=txdb, gene_list=gene_list)
        } else {
            ## Then species name
            message("A species name.")
        }
    } else if (class(length_db)[[1]] == "AnnotationDbi") {
        stop("This currently requires an actual OrganismDb, not AnnotationDbi.")
    } else if (class(length_db)[[1]] == "OrgDb") {
        stop("OrgDb objects contain links to other databases, but sadly are missing gene lengths.")
    } else if (class(length_db)[[1]] == "OrganismDb" | class(length_db)[[1]] == "AnnotationDbi") {
        ## metadf <- extract_lengths(db=length_db, gene_list=gene_list)
        metadf <- sm(extract_lengths(db=length_db, gene_list=gene_list, ...))
    } else if (class(length_db)[[1]] == "TxDb") {
        metadf <- sm(extract_lengths(db=length_db, gene_list=gene_list, ...))
    } else if (class(length_db)[[1]] == "data.frame") {
        metadf <- length_db
    } else {
        stop("This requires either the name of a goseq supported species or an orgdb instance.")
    }
    ## Sometimes the column with gene lengths is named 'width'
    ## In that case, fix it.
    if (is.null(metadf[["width"]]) & is.null(metadf[["length"]])) {
        stop("The length db needs to have a length or width column.")
    } else if (is.null(metadf[["length"]])) {
        ## Then it is named 'width' and I want to rename it to length
        colnames(metadf) <- gsub(x=colnames(metadf), pattern="width", replacement="length")
    }
    ## Now I should have the gene list and gene lengths

    godf <- data.frame()
    if (class(go_db) == "character") {
        ## A text table or species name
        if (grepl(pattern="\\.csv", x=go_db, perl=TRUE) |
            grepl(pattern="\\.tab", x=go_db, perl=TRUE)) {
            ## table
            godf <- read.table(go_db, ...)
            colnames(godf) <- c("ID", "GO")
        } else {
            ## Assume species name
            supported <- TRUE
            species <- go_db
        }
    } else if (class(go_db)[[1]] == "OrganismDb") {
        godf <- extract_go(go_db)
    } else if (class(go_db)[[1]] == "OrgDb") {
        godf <- extract_go(go_db)
    } else if (class(go_db)[[1]] == "data.frame") {
        godf <- go_db
        godf <- godf[, c("ID", "GO")]
    } else {
        message("Not sure what to do here.")
    }
    ## entrez IDs are numeric.  This is a problem when doing the pwf function because it sets
    ## the rownames to the IDs.  As a result, we need to call make.names() on them.
    godf[["ID"]] <- make.names(godf[["ID"]])
    metadf[["ID"]] <- make.names(metadf[["ID"]])
    de_genelist[["ID"]] <- make.names(de_genelist[["ID"]])
    ## Ok, now I have a df of GOids, all gene lengths, and DE gene list. That is everything
    ## I am supposed to need for goseq.

    ## See how many entries from the godb are in the list of genes.
    id_xref <- de_genelist[["ID"]] %in% godf[["ID"]]
    message(paste0("Found ", sum(id_xref), " genes from the sig_genes in the go_db."))

    ## So lets merge the de genes and gene lengths to ensure that they are consistent.
    ## Then make the vectors expected by goseq
    merged_ids_lengths <- metadf
    ## The following line was done in order to avoid
    ## "'unimplemented type 'list' in 'orderVector1'"
    merged_ids_lengths[["ID"]] <- as.character(merged_ids_lengths[["ID"]])
    merged_ids_lengths <- merge(merged_ids_lengths, de_genelist, by.x="ID", by.y="ID", all.x=TRUE)
    merged_ids_lengths[is.na(merged_ids_lengths)] <- 0
    ## Not casing the next lines as character/numeric causes weird errors like 'names' attribute
    ## must be the same length as the vector
    de_vector <- as.vector(as.numeric(merged_ids_lengths[["DE"]]))
    names(de_vector) <- make.names(as.character(merged_ids_lengths[["ID"]]), unique=TRUE)
    length_vector <- as.vector(as.numeric(merged_ids_lengths[["length"]]))
    names(length_vector) <- make.names(as.character(merged_ids_lengths[["ID"]]), unique=TRUE)

    pwf_plot <- NULL
    pwf <- goseq::nullp(DEgenes=de_vector, bias.data=length_vector, plot.fit=doplot)
    if (isTRUE(doplot)) {
        pwf_plot <- recordPlot()
    }
    godata <- goseq::goseq(pwf, gene2cat=godf, use_genes_without_cat=TRUE, method=goseq_method)

    goseq_p <- try(plot_histogram(godata[["over_represented_pvalue"]], bins=50))
    ## goseq_p_second <- sort(unique(table(goseq_p[["data"]])), decreasing=TRUE)[2]
    goseq_p_nearzero <- table(goseq_p[["data"]])[[1]]
    ## Set the y scale to 2x the second highest number
    ## (assuming always that the highest is a p-value of 1)
    goseq_y_limit <- goseq_p_nearzero * 2
    goseq_p <- goseq_p + ggplot2::scale_y_continuous(limits=c(0, goseq_y_limit))
    message("simple_goseq(): Calculating q-values")
    qdata <- godata[["over_represented_pvalue"]]
    qdata[qdata > 1] <- 1 ## For scientific numbers which are 1.0000E+00 it might evaluate to 1.0000000000000001
    qvalues <- tryCatch({
        ttmp <- as.numeric(qdata)
        ttmp <- qvalue::qvalue(ttmp)[["qvalues"]]
    },
    error=function(cond) {
        message(paste0("The qvalue estimate failed."))
        return(1)
    },
    finally={
    })
    message("Using GO.db to extract terms and categories.")
    godata[["term"]] <- goterm(godata[["category"]])
    godata[["ontology"]] <- goont(godata[["category"]])
    godata <- cbind(godata, qvalues)
    colnames(godata) <- c("category", "over_represented_pvalue", "under_represented_pvalue",
                          "numDEInCat", "numInCat", "term", "ontology", "qvalue")
    if (is.null(adjust)) {
        godata_interesting <- subset(godata, godata[["over_represented_pvalue"]] <= pvalue)
        padjust_method <- "none"
    } else {
        ## There is a requested pvalue adjustment
        godata_interesting <- subset(godata, stats::p.adjust(godata[["over_represented_pvalue"]],
                                                             method=padjust_method) <= adjust)
        if (dim(godata_interesting)[1] < minimum_interesting) {
            message(paste("simple_goseq(): There are no genes with an adj.p<", adjust, " using: ",
                          padjust_method, ".", sep=""))
            message(sprintf("simple_goseq(): Providing genes with raw pvalue<%s", pvalue))
            godata_interesting <- subset(godata, godata[["over_represented_pvalue"]] <= pvalue)
            padjust_method <- "none"
        }
    }
    message("simple_goseq(): Filling godata with terms, this is slow.")
    godata_interesting <- goseq_table(godata_interesting)
    message("simple_goseq(): Making pvalue plots for the ontologies.")
    pvalue_plots <- plot_goseq_pval(godata)
    ## mf_subset <- subset(godata, ontology == "MF")

    mf_subset <- godata[godata[["ontology"]] == "MF", ]
    rownames(mf_subset) <- mf_subset[["category"]]
    ##bp_subset <- subset(godata, ontology == "BP")
    bp_subset <- godata[godata[["ontology"]] == "BP", ]
    rownames(bp_subset) <- bp_subset[["category"]]
    ## cc_subset <- subset(godata, ontology == "CC")
    cc_subset <- godata[godata[["ontology"]] == "CC", ]
    rownames(cc_subset) <- cc_subset[["category"]]
    ## mf_interesting <- subset(godata_interesting, ontology == "MF")

    mf_interesting <- godata_interesting[godata_interesting[["ontology"]] == "MF", ]
    rownames(mf_interesting) <- mf_interesting[["category"]]
    mf_interesting <- mf_interesting[, c("ontology", "numDEInCat", "numInCat",
                                         "over_represented_pvalue", "qvalue", "term")]
    ##bp_interesting <- subset(godata_interesting, ontology == "BP")
    bp_interesting <- godata_interesting[godata_interesting[["ontology"]] == "BP", ]
    rownames(bp_interesting) <- bp_interesting[["category"]]
    bp_interesting <- bp_interesting[, c("ontology", "numDEInCat", "numInCat",
                                         "over_represented_pvalue", "qvalue", "term")]
    ##cc_interesting <- subset(godata_interesting, ontology == "CC")
    cc_interesting <- godata_interesting[godata_interesting[["ontology"]] == "CC", ]
    rownames(cc_interesting) <- cc_interesting[["category"]]
    cc_interesting <- cc_interesting[, c("ontology", "numDEInCat", "numInCat",
                                         "over_represented_pvalue", "qvalue", "term")]

    pval_plots <- list(
        "bpp_plot_over" = pvalue_plots[["bpp_plot_over"]],
        "mfp_plot_over" = pvalue_plots[["mfp_plot_over"]],
        "ccp_plot_over" = pvalue_plots[["ccp_plot_over"]])

    return_list <- list("input" = sig_genes,
                        "pwf" = pwf,
                        "pwf_plot" = pwf_plot,
                        "alldata" = godata,
                        "godf" = godf,
                        "pvalue_histogram" = goseq_p,
                        "godata_interesting" = godata_interesting,
                        "mf_interesting" = mf_interesting,
                        "bp_interesting" = bp_interesting,
                        "cc_interesting" = cc_interesting,
                        "goadjust_method" = goseq_method,
                        "adjust_method" = padjust_method,
                        "mf_subset" = mf_subset,
                        "pvalue_plots" = pval_plots,
                        "bp_subset" = bp_subset,
                        "cc_subset" = cc_subset,
                        "qdata" = qdata)
    return(return_list)
}

#' Given a set of goseq data from simple_goseq(),
#' make a list of genes represented in each ontology.
#'
#' This function uses the GO2ALLEG data structure to reverse map ontology
#' categories to a list of genes represented. It therefore assumes that the
#' GO2ALLEG.rda data structure has been deposited in pwd().  This in turn
#' may be generated by clusterProfilers buildGOmap() function if it doesn't
#' exist.  For some species it may also be auto-generated.
#' With little work this can be made much more generic, and it probably should.
#'
#' @param goseq List of goseq specific results as generated by simple_goseq().
#' @param ontology Ontology to search (MF/BP/CC).
#' @param pval Maximum accepted pvalue to include in the list of categories to cross reference.
#' @param include_all Include all genes in the ontology search?
#' @param ... Extra options without a purpose just yet.
#' @return Data frame of categories/genes.
#' @seealso \pkg{goseq} \pkg{clusterProfiler}
#'  \code{\link{simple_goseq}}
#' @examples
#' \dontrun{
#'  data <- simple_goseq(sig_genes=limma_output, lengths=annotation_df, goids=goids_df)
#'  genes_in_cats <- gather_genes(data, ont='BP')
#' }
#' @export
gather_goseq_genes <- function(goseq, ontology=NULL, pval=0.1, include_all=FALSE, ...) {
    arglist <- list(...)
    categories <- NULL
    if (is.null(ontology)) {
        retlist <- list()
        message("No ontology provided, performing all.")
        for (type in c("MF", "BP", "CC")) {
            retlist[[type]] <- gather_goseq_genes(goseq, ontology=type,
                                                  pval=pval, include_all=include_all, ...)
        }
        return(retlist)
    } else if (ontology == "MF") {
        categories <- goseq[["mf_subset"]]
    } else if (ontology == "BP") {
        categories <- goseq[["bp_subset"]]
    } else if (ontology == "CC") {
        categories <- goseq[["cc_subset"]]
    } else {
        retlist <- list()
        message("No ontology provided, performing all.")
        for (type in c("MF", "BP", "CC")) {
            retlist[[type]] <- gather_goseq_genes(goseq, ontology=type,
                                                  pval=pval, include_all=include_all, ...)
        }
        return(retlist)
    }
    input <- goseq[["input"]]
    ##categories <- subset(categories, over_represented_pvalue <= pval)
    categories <- categories[ categories[["over_represented_pvalue"]] <= pval, ]
    cats <- rownames(categories)
    godf <- goseq[["godf"]]
    genes_per_ont <- function(cat) {
        ## all_entries <- subset(godf, GO==cat)[["ID"]]
        colnames(godf) <- c("ID", "GO")
        ## Only keep the set of entries which are filled in.
        godf <- godf[complete.cases(godf), ]
        ## Pull all rows which are of our category.
        found_idx <- godf[["GO"]] == cat
        ## Then extract those rows from the full set of go mappings
        foundlings <- godf[found_idx, ]
        ## Finally, pull those gene IDs
        all_entries <- unique(foundlings[["ID"]])
        ## Extract the limma logFC for all genes.
        all_names <- toString(all_entries)
        ## Now find the 'significant' genes in this set.
        entries_in_sig <- input[ rownames(input) %in% all_entries, ]
        ## And get their names.
        sig_names <- toString(as.character(rownames(entries_in_sig)))
        ## Along with the logFC from limma
        sig_limma <- toString(as.character(entries_in_sig[, "limma_logfc"]))
        sig_deseq <- toString(as.character(entries_in_sig[, "deseq_logfc"]))
        sig_edger <- toString(as.character(entries_in_sig[, "edger_logfc"]))
        sig_basic <- toString(as.character(entries_in_sig[, "basic_logfc"]))
        ##names <- toString(as.character(rownames(entries_in_input)))
        retlist <- list("all" = all_names,
                        "sig" = sig_names,
                        "limma_sigfc" = sig_limma,
                        "edger_sigfc" = sig_edger,
                        "deseq_sigfc" = sig_deseq,
                        "basic_sigfc" = sig_basic)
        return(retlist)
    }
    gene_list <- lapply(cats, genes_per_ont)
    names(gene_list) <- cats
    gene_df <- data.table::rbindlist(gene_list)
    rownames(gene_df) <- cats
    return(gene_df)
}

## EOF
