#' Enhance the goseq table of gene ontology information.
#'
#' While goseq has some nice functionality, the table of outputs it provides is somewhat lacking.
#' This attempts to increase that with some extra helpful data like ontology categories,
#' definitions, etc.
#'
#' @param df Dataframe of ontology information.  This is intended to be the output from goseq
#'     including information like numbers/category, GOids, etc.  It requires a column 'category'
#'     which contains: GO:000001 and such.
#' @param file Csv file to which to write the table.
#' @return Ontology table with annotation information included.
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
    df$good <- gotest(df[["category"]])
    message("Removing undefined categories.")
    ## df = subset(df, good == 1)
    df <- df[ which(df[["good"]] == 1), ]
    message("Gathering synonyms.")
    df$synonym <- gosyn(df[["category"]])
    ##message("Gathering secondary ids.")
    ##secondary <- try(gosec(df$category), silent=TRUE)
    ##if (class(secondary) != 'try-error') {
    ##    df$secondary <- secondary
    ##}
##    print("Gather approximate go levels.")  ## This function is too slow, commented it out.
##    df$level = golevel(df$categoy)
    message("Gathering category definitions.")
    df[["definition"]] <- godef(df[["category"]])
    df <- df[, c("category","numDEInCat","numInCat","over_represented_pvalue",
                 "under_represented_pvalue","qvalue","ontology","term",
                 "synonym","definition")]
    if (!is.null(file)) {
        write.csv(df, file=file)
    }
    return(df)
}


extract_lengths <- function(db=NULL, gene_list=NULL,
                            types=c("GenomicFeatures::genes",
                                    "GenomicFeatures::cds",
                                    "GenomicFeatures::transcripts")) {
    tmpdb <- db
    metadf <- NULL
    success <- FALSE
    for (type in types) {
        test_string <- paste0("testing <- ", type, "(tmpdb)")
        eval(parse(text=test_string))
        ## as.data.frame is not only base, but also biocgenerics!!!
        test_meta <- BiocGenerics::as.data.frame(testing)
        overlap <- gene_list %in% rownames(test_meta)
        if (sum(overlap) == 0) {
            ## sometimes the rownames get set to 1..rows, in that case, try the last column which is the transcript name
            overlap <- gene_list %in% test_meta[, ncol(test_meta)]
            rownames(test_meta) <- test_meta[, ncol(test_meta)]
        }
        if (sum(overlap) == length(gene_list)) {
            success <- TRUE
            if (!is.null(test_meta[["width"]])) {
                metadf <- as.data.frame(cbind(rownames(test_meta), test_meta[["width"]]))
            } else if (!is.null(test_meta[["length"]])) {
                metadf <- as.data.frame(cbind(rownames(test_meta), test_meta[["length"]]))
            } else {
                stop("This requires either length or width columns.")
            }
            colnames(metadf) <- c("ID","length")
            rownames(metadf) <- metadf[["ID"]]
            message("Returning metadf.")
            return(metadf)
        }
    }
    if (!isTRUE(success)) {
        stop("Unable to get the set of lengths.")
    }
    return(metadf)
}

## The last step is to gather GO ids
extract_go <- function(db) {
    possible_keytypes <- keytypes(db)
    godf <- data.frame()
    success <- FALSE
    for (type in possible_keytypes) {
        test <- head(keys(x=db, keytype=type))
        found <- sum(test %in% metadf[["ID"]])
        if (found == length(test)) {
            ## Use this keytype to get geneIDs
            ids <- keys(x=db, keytype=type)
            if ("GOID" %in% possible_keytypes) {
                godf <- select(x=db, keys=ids, keytype=type, columns=c("GOID"))
                godf[["ID"]] <- godf[[1]]
                godf <- godf[, c("ID","GOID")]
                colnames(godf) <- c("ID","GO")
            } else if ("GO" %in% possible_keytypes) {
                godf <- select(x=db, keys=ids, keytype=type, columns=c("GO"))
                godf[["ID"]] <- godf[[1]]
                godf <- godf[, c("ID","GO")]
                colnames(godf) <- c("ID","GO")
            } else if ("GOALL" %in% possible_keytypes) {
                godf <- select(x=db, keys=ids, keytype=type, columns=c("GOALL"))
                godf[["ID"]] <- godf[[1]]
                godf <- godf[, c("ID","GO")]
                colnames(godf) <- c("ID","GO")
            }
            return(godf)
        }
    }
}

write_new_goseq_xlsx <- function(goseq, file="excel/goseq.xlsx", pval=0.1, add_plots=TRUE) {
    excel_dir <- dirname(file)
    if (!file.exists(excel_dir)) {
        dir.create(excel_dir, recursive=TRUE)
    }

    ## Pull out the relevant portions of the goseq data
    ## For this I am using the same (arbitrary) rules as in gather_goseq_genes()
    goseq_mf <- goseq[["mf_subset"]]
    goseq_mf <- goseq_mf[ goseq_mf[["over_represented_pvalue"]] <= pval, ]
    goseq_mf_genes <- gather_goseq_genes(goseq, ontology="MF", pval=pval)
    mf_genes <- as.data.frame(goseq_mf_genes)
    rownames(mf_genes) <- rownames(goseq_mf_genes)
    goseq_mf <- merge(goseq_mf, mf_genes, by="row.names")
    rownames(goseq_mf) <- goseq_mf[["Row.names"]]
    goseq_mf <- goseq_mf[-1]

    goseq_bp <- goseq[["bp_subset"]]
    goseq_bp <- goseq_bp[ goseq_bp[["over_represented_pvalue"]] <= pval, ]
    goseq_bp_genes <- gather_goseq_genes(goseq, ontology="BP", pval=pval)
    bp_genes <- as.data.frame(goseq_bp_genes)
    rownames(bp_genes) <- rownames(goseq_bp_genes)
    goseq_bp <- merge(goseq_bp, bp_genes, by="row.names")
    rownames(goseq_bp) <- goseq_bp[["Row.names"]]
    goseq_bp <- goseq_bp[-1]

    goseq_cc <- goseq[["cc_subset"]]
    goseq_cc <- goseq_cc[ goseq_cc[["over_represented_pvalue"]] <= pval, ]
    goseq_cc_genes <- gather_goseq_genes(goseq, ontology="CC", pval=pval)
    cc_genes <- as.data.frame(goseq_cc_genes)
    rownames(cc_genes) <- rownames(goseq_cc_genes)
    goseq_cc <- merge(goseq_cc, cc_genes, by="row.names")
    rownames(goseq_cc) <- goseq_cc[["Row.names"]]
    goseq_cc <- goseq_cc[-1]

    goseq_mf <- goseq_mf[, c(7,1,6,2,8,10,9,4,5)]
    goseq_bp <- goseq_bp[, c(7,1,6,2,8,10,9,4,5)]
    goseq_cc <- goseq_cc[, c(7,1,6,2,8,10,9,4,5)]
    colnames(goseq_mf) <- c("Ontology","Category","Term","Over p-value", "Q-value",
                            "DE genes in cat", "All genes in cat", "Num. DE", "Num. in cat.")
    colnames(goseq_bp) <- c("Ontology","Category","Term","Over p-value", "Q-value",
                            "DE genes in cat", "All genes in cat", "Num. DE", "Num. in cat.")
    colnames(goseq_cc) <- c("Ontology","Category","Term","Over p-value", "Q-value",
                            "DE genes in cat", "All genes in cat", "Num. DE", "Num. in cat.")

    wb <- openxlsx::createWorkbook(creator="atb")
    hs1 <- openxlsx::createStyle(fontColour="#000000", halign="LEFT", textDecoration="bold", border="Bottom", fontSize="30")
    ## This stanza will be repeated so I am just incrementing the new_row

    new_row <- 1
    sheet <- "BP"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "BP Results from goseq.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=goseq_bp, tableStyle="TableStyleMedium9", startRow=new_row)
    ## I want to add the pvalue plots, which are fairly deeply embedded in kept_ontology
    if (isTRUE(add_plots)) {
        a_plot <- goseq[["pvalue_plots"]][["bpp_plot_over"]]
        print(a_plot)
        openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(goseq_bp) + 2, startRow=new_row, fileType="png", units="in")
    }
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    new_row <- new_row + nrow(goseq_bp) + 2


    new_row <- 1
    sheet <- "MF"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "MF Results from goseq.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=goseq_mf, tableStyle="TableStyleMedium9", startRow=new_row)
    if (isTRUE(add_plots)) {
        a_plot <- goseq[["pvalue_plots"]][["mfp_plot_over"]]
        print(a_plot)
        openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(goseq_mf) + 2, startRow=new_row, fileType="png", units="in")
    }
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    new_row <- new_row + nrow(goseq_mf) + 2

    new_row <- 1
    sheet <- "CC"
    openxlsx::addWorksheet(wb, sheetName=sheet)
    openxlsx::writeData(wb, sheet, "CC Results from goseq.", startRow=new_row)
    openxlsx::addStyle(wb, sheet, hs1, new_row, 1)
    new_row <- new_row + 1
    openxlsx::writeDataTable(wb, sheet, x=goseq_cc, tableStyle="TableStyleMedium9", startRow=new_row)
    if (isTRUE(add_plots)) {
        a_plot <- goseq[["pvalue_plots"]][["ccp_plot_over"]]
        print(a_plot)
        openxlsx::insertPlot(wb, sheet, width=6, height=6, startCol=ncol(goseq_cc) + 2, startRow=new_row, fileType="png", units="in")
    }
    openxlsx::setColWidths(wb, sheet=sheet, cols=2:7, widths="auto")
    new_row <- new_row + nrow(goseq_cc) + 2

    res <- openxlsx::saveWorkbook(wb, file, overwrite=TRUE)
    return(res)
}

new_simple_goseq <- function(de_genes, go_db, length_db, doplot=TRUE,
                             adjust=0.1, pvalue=0.1, qvalue=0.1,
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
    ## de_genes may be a list, character list, or data frame.
    gene_list <- NULL
    if (class(de_genes) == "character") { ## Then this is a character list of gene ids
        gene_list <- de_genes
    } else if (class(de_genes) == "list") {
        gene_list <- names(de_genes)
    } else if (class(de_genes) == "data.frame") {
        if (is.null(rownames(de_genes)) & is.null(de_genes[["ID"]])) {
            stop("This requires a set of gene IDs either from the rownames or a column named 'ID'.")
        } else if (!is.null(rownames(de_genes))) {
            gene_list <- rownames(de_genes)
        } else {
            gene_list <- de_genes[["ID"]]
        }
    } else {
        stop("Not sure how to handle your set of gene ids.")
    }
    ## At this point I should have a character list of gene ids named 'gene_list'
    de_genelist <- as.data.frame(gene_list)
    de_genelist[["DE"]] <- 1
    colnames(de_genelist) <- c("ID","DE")

    ## Database of lengths may be a gff file, TxDb, or OrganismDb
    metadf <- NULL
    if (class(length_db) == "character")  {  ## Then this should be either a gff file or species name.
        if (grepl(pattern="\\.gff", x=length_db, perl=TRUE) | grepl(pattern="\\.gtf", x=length_db, perl=TRUE)) { ## gff file
            txdb <- GenomicFeatures::makeTxDbFromGFF(length_db)
            metadf <- extract_lengths(db=txdb, gene_list=gene_list)
        } else {  ## Then species name
            message("A species name.")
        }
    } else if (class(length_db)[[1]] == "AnnotationDbi") {
        stop("This currently requires an actual OrganismDb, not AnnotationDbi.")
    } else if (class(length_db)[[1]] == "OrgDb") {
        stop("OrgDb objects contain links to other databases, but sadly are missing gene lengths.")
    } else if (class(length_db)[[1]] == "OrganismDb" | class(length_db)[[1]] == "AnnotationDbi") {
        metadf <- extract_lengths(db=length_db, gene_list=gene_list)
        stop("I can't work with OrganismDb just yet.")
    } else if (class(length_db)[[1]] == "TxDb") {
        metadf <- extract_lengths(db=length_db, gene_list=gene_list)
    } else {
        stop("This requires either the name of a goseq supported species or an orgdb instance.")
    }
    ## Now I should have the gene list and gene lengths

    godf <- data.frame()
    if (class(go_db) == "character") {  ## A text table or species name
        if (grepl(pattern="\\.csv", x=go_db, perl=TRUE) | grepl(pattern="\\.tab", x=go_db, perl=TRUE)) { ## table
            godf <- read.table(go_db, ...)
            colnames(godf) <- c("ID","GO")
        } else {  ## Assume species name
            message("Using goseq supported species name.")
            supported <- TRUE
            species <- go_db
        }
    } else if (class(go_db)[[1]] == "OrganismDb") {
        message("Using Organismdb")
        godf <- extract_go(go_db)
    } else if (class(go_db)[[1]] == "OrgDb") {
        message("Using Orgdb.")
        godf <- extract_go(go_db)
    } else {
        message("Not sure what to do here.")
    }
    ## Ok, now I have a df of GOids, all gene lengths, and DE gene list. That is everything I am supposed to need for goseq.

    ## So lets merge the de genes and gene lengths to ensure that they are consistent.
    ## Then make the vectors expected by goseq
    merged_ids_lengths <- metadf
    merged_ids_lengths <- merge(merged_ids_lengths, de_genelist, by.x="ID", by.y="ID", all.x=TRUE)
    merged_ids_lengths[is.na(merged_ids_lengths)] <- 0
    ## Not casing the next lines as character/numeric causes weird errors like 'names' attribute must be the same length as the vector
    de_vector <- as.vector(as.numeric(merged_ids_lengths[["DE"]]))
    names(de_vector) <- as.character(merged_ids_lengths[["ID"]])
    length_vector <- as.vector(as.numeric(merged_ids_lengths[["length"]]))
    names(length_vector) <- as.character(merged_ids_lengths[["ID"]])

    pwf_plot <- NULL
    pwf <- goseq::nullp(DEgenes=de_vector, bias.data=length_vector, plot.fit=doplot)
    if (isTRUE(doplot)) {
        pwf_plot <- recordPlot()
    }
    godata <- goseq::goseq(pwf, gene2cat=godf, use_genes_without_cat=TRUE, method=goseq_method)

    goseq_p <- try(plot_histogram(godata$over_represented_pvalue, bins=20))
    goseq_p_second <- sort(unique(table(goseq_p[["data"]])), decreasing=TRUE)[2]
    ## Set the y scale to 2x the second highest number
    ## (assuming always that the highest is a p-value of 1)
    goseq_y_limit <- goseq_p_second * 2
    goseq_p <- goseq_p + ggplot2::scale_y_continuous(limits=c(0, goseq_y_limit))
    message("simple_goseq(): Calculating q-values")
    qdata <- godata[["over_represented_pvalue"]]
    qdata[qdata > 1] <- 1 ## For scientific numbers which are 1.0000E+00 it might evaluate to 1.0000000000000001
    qvalues <- tryCatch(
    {
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
    colnames(godata) <- c("category","over_represented_pvalue","under_represented_pvalue",
                          "numDEInCat","numInCat","term","ontology","qvalue")
    if (is.null(adjust)) {
        godata_interesting <- subset(godata, godata$over_represented_pvalue <= pvalue)
        padjust_method <- "none"
    } else {  ## There is a requested pvalue adjustment
        godata_interesting <- subset(godata, stats::p.adjust(godata[["over_represented_pvalue"]], method=padjust_method) <= adjust)
        if (dim(godata_interesting)[1] < minimum_interesting) {
            message(paste("simple_goseq(): There are no genes with an adj.p<", adjust, " using: ", padjust_method, ".", sep=""))
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
    mf_interesting <- mf_interesting[, c("ontology","numDEInCat","numInCat","over_represented_pvalue","qvalue","term")]
    ##bp_interesting <- subset(godata_interesting, ontology == "BP")
    bp_interesting <- godata_interesting[godata_interesting[["ontology"]] == "BP", ]
    rownames(bp_interesting) <- bp_interesting[["category"]]
    bp_interesting <- bp_interesting[ ,c("ontology","numDEInCat","numInCat","over_represented_pvalue","qvalue","term")]
    ##cc_interesting <- subset(godata_interesting, ontology == "CC")
    cc_interesting <- godata_interesting[godata_interesting[["ontology"]] == "CC", ]
    rownames(cc_interesting) <- cc_interesting[["category"]]
    cc_interesting <- cc_interesting[, c("ontology","numDEInCat","numInCat","over_represented_pvalue","qvalue","term")]

    pval_plots <- list(
        "bpp_plot_over" = pvalue_plots[["bpp_plot_over"]],
        "mfp_plot_over" = pvalue_plots[["mfp_plot_over"]],
        "ccp_plot_over" = pvalue_plots[["ccp_plot_over"]])

    return_list <- list("input" = de_genes,
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
#' @param goseq_data List of goseq specific results as generated by simple_goseq().
#' @param ontology Ontology to search (MF/BP/CC).
#' @param pval Maximum accepted pvalue to include in the list of categories to cross reference.
#' @param include_all Include all genes in the ontology search?
#' @param ... Extra options without a purpose just yet.
#' @return Data frame of categories/genes.
#' @seealso \link{simple_goseq} \code{\link[clusterProfiler]{buildGOmap}},
#' @examples
#' \dontrun{
#'  data = simple_goseq(de_genes=limma_output, lengths=annotation_df, goids=goids_df)
#'  genes_in_cats = gather_genes(data, ont='BP')
#' }
#' @export
gather_goseq_genes <- function(goseq_data, ontology=NULL, pval=0.1, include_all=FALSE, ...) {
    arglist <- list(...)
    categories <- NULL
    if (is.null(ontology)) {
        retlist <- list()
        message("No ontology provided, performing all.")
        for (type in c("MF","BP","CC")) {
            retlist[[type]] <- gather_goseq_genes(goseq_data, ontology=type, pval=pval, include_all=include_all, ...)
        }
        return(retlist)
    } else if (ontology == "MF") {
        categories <- goseq_data[["mf_subset"]]
    } else if (ontology == "BP") {
        categories <- goseq_data[["bp_subset"]]
    } else if (ontology == "CC") {
        categories <- goseq_data[["cc_subset"]]
    } else {
        retlist <- list()
        message("No ontology provided, performing all.")
        for (type in c("MF","BP","CC")) {
            retlist[[type]] <- gather_goseq_genes(goseq_data, ontology=type, pval=pval, include_all=include_all, ...)
        }
        return(retlist)
    }
    input <- goseq_data[["input"]]
    ##categories <- subset(categories, over_represented_pvalue <= pval)
    categories <- categories[ categories[["over_represented_pvalue"]] <= pval, ]
    cats <- rownames(categories)
    godf <- goseq_data[["godf"]]
    genes_per_ont <- function(cat) {
        all_entries <- subset(godf, GO==cat)[["ID"]]
        entries_in_input <- input[ rownames(input) %in% all_entries, ]
        names <- toString(as.character(rownames(entries_in_input)))
        all <- toString(all_entries)
        retlist <- list("all" = all,
                        "sig" = names)
        return(retlist)
    }
    gene_list <- lapply(cats, genes_per_ont)
    names(gene_list) <- cats
    gene_df <- data.table::rbindlist(gene_list)
    rownames(gene_df) <- cats
    return(gene_df)
}

#' Perform a simplified goseq analysis.
#'
#' goseq can be pretty difficult to get set up for non-supported organisms.  This attempts to make
#' that process a bit simpler as well as give some standard outputs which should be similar to those
#' returned by clusterprofiler/topgo/gostats/gprofiler.
#'
#' @param de_genes Data frame of differentially expressed genes, containing IDs etc.
#' @param all_genes Universe of possible genes.
#' @param lengths Length of each gene with an ID in de_genes.
#' @param goids_df List of ontology accessions to gene accessions.
#' @param doplot Include pwf plots?
#' @param adjust Minimum adjusted pvalue for 'significant.'
#' @param pvalue Minimum pvalue for 'significant.'
#' @param qvalue Minimum qvalue for 'significant.'
#' @param goseq_method Statistical test for goseq to use.
#' @param padjust_method Which method to use to adjust the pvalues.
#' @param species Optionally choose a species from supportedOrganisms().
#' @param length_db Source of gene lengths?
#' @param gff gff file source of gene lengths.
#' @param ... Extra parameters which I do not recall
#' @return Big list including:
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
simple_goseq <- function(de_genes, all_genes=NULL, lengths=NULL, goids_df=NULL, doplot=TRUE,
                         adjust=0.1, pvalue=0.1, qvalue=0.1, goseq_method="Wallenius",
                         padjust_method="BH", species=NULL, length_db="ensGene", gff=NULL, ...) {
    message("simple_goseq() makes some pretty hard assumptions about the data it is fed:")
    message("It requires 2 tables, one of GOids which must have columns (gene)ID and GO(category)")
    message("The other table is of gene lengths with columns (gene)ID and (gene)width.")
    message("Other columns are fine, but ignored.")
    arglist <- list(...)
    if (is.null(arglist[["minimum_interesting"]])) {
        arglist[["minimum_interesting"]] = 10
    }
    if (is.null(lengths) & is.null(gff) & is.null(species)) {
        stop("simple_goseq() requires either a length dataframe or gff file for gene lengths; or a supported species.")
    } else if (!is.null(lengths)) {
        message("simple_goseq(): Using the explicit lengths df for gene lengths.")
    } else if (!is.null(species)) {
        message(paste0("simple_goseq(): Hopefully ", species, " is supported by goseq."))
    } else {
        ## This is probably hopelessly fragile and requires further thought
        length_df <- gff2df(gff)
        lengths <- length_df[length_df[["type"]] == 'CDS', ]
        lengths <- length_df[, c("gene_id", "width")]
        colnames(lengths) <- c("ID","width")
    }
    if (is.null(de_genes[["ID"]])) {
        de_genes[["ID"]] <- make.names(rownames(de_genes), unique=TRUE)
    }
    if (is.null(de_genes[["DE"]])) {
        de_genes[["DE"]] <- 1
    }
    de_vector <- NULL
    de_table <- de_genes[, c("ID","DE")]
    if (is.null(lengths) & is.null(all_genes) & is.null(species)) {
        stop("simple_goseq(): Need either a set of all genes or gene lengths")
    } else if (!is.null(lengths)) {
        message("simple_goseq(): Using the length data to fill in the de vector.")
        de_table <- merge(de_table, lengths, by.x="ID", by.y="ID", all.y=TRUE)
        de_table[is.na(de_table)] <- 0  ## Set the new entries DE status to 0
        rownames(de_table) <- make.names(de_table[["ID"]], unique=TRUE)
        de_vector <- as.vector(de_table[["DE"]])
        names(de_vector) <- rownames(de_table)
    } else if (!is.null(species)) {
        message("simple_goseq(): Using species and length_db to get metadata.")
        if (species == 'hsapiens') {
            message("Replacing hsapiens with hg19.")
            species <- "hg19"
        }
        db_name <- paste0(species, ".", length_db, ".LENGTH")
        db_invocation <- paste0("data(", db_name, ", package='geneLenDataBase')")
        eval(parse(text=db_invocation))

        gene_names <- as.data.frame(get(paste(species, length_db, "LENGTH", sep = "."))$Gene)
        colnames(gene_names) <- c("ID")
        de_table <- merge(de_table, gene_names, by.x="ID", by.y="ID", all.y=TRUE)
        de_table[is.na(de_table)] <- 0
        rownames(de_table) <- make.names(de_table[["ID"]], unique=TRUE)
        de_vector <- as.vector(de_table[["DE"]])
        names(de_vector) <- rownames(de_table)
    } else { ## If both lengths and all_genes are defined, use all_genes.
        message("simple_goseq(): Using all genes to fill in the de vector.")
        de_table <- merge(de_table, all_genes, by.x="ID", by.y="row.names", all.y=TRUE)
        ##de_table[is.na(de_table)] <- 0  ## Set the new entries DE status to 0
        de_table <- merge(de_table, lengths, by.x="ID", by.y="ID", all.x=TRUE)
        de_table[["DE"]][is.na(de_table[["DE"]])] <- 0  ## Set the new entries DE status to 0
        rownames(de_table) <- make.names(de_table[["ID"]], unique=TRUE)
        de_vector <- as.vector(de_table[["DE"]])
        names(de_vector) <- rownames(de_table)
    }
    pwf <- NULL
    if (is.null(species)) {
        ## length_table = lengths[,c("ID","width")]
        if (!is.null(de_table[["length"]])) {
            de_table[["width"]] <- de_table[["length"]]
        }
        width_vector <- as.vector(de_table[["width"]])
        names(width_vector) <- de_table[["ID"]]
        if (is.null(goids_df)) {
            stop("simple_goseq(): The goids are not defined.")
        }
        goids_df <- goids_df[, c("ID","GO")]
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
        godata <- goseq::goseq(pwf, gene2cat=goids_df, use_genes_without_cat=TRUE, method=goseq_method)
    } else {
        godata <- goseq::goseq(pwf, species, length_db, use_genes_without_cat=TRUE, method=goseq_method)
    }
    goseq_p <- try(plot_histogram(godata$over_represented_pvalue, bins=20))
    goseq_p_second <- sort(unique(table(goseq_p[["data"]])), decreasing=TRUE)[2]
    ## Set the y scale to 2x the second highest number
    ## (assuming always that the highest is a p-value of 1)
    goseq_y_limit <- goseq_p_second * 2
    goseq_p <- goseq_p + ggplot2::scale_y_continuous(limits=c(0, goseq_y_limit))
    message("simple_goseq(): Calculating q-values")
    qdata <- godata[["over_represented_pvalue"]]
    qdata[qdata > 1] <- 1 ## For scientific numbers which are 1.0000E+00 it might evaluate to 1.0000000000000001
    qvalues <- tryCatch(
    {
        ttmp <- as.numeric(qdata)
        ttmp <- qvalue::qvalue(ttmp)[["qvalues"]]
    },
    error=function(cond) {
        message(paste0("The qvalue estimate failed."))
        return(1)
    },
    finally={
    })
    godata[["term"]] <- goterm(godata[["category"]])
    godata[["ontology"]] <- goont(godata[["category"]])
    godata <- cbind(godata, qvalues)
    colnames(godata) <- c("category","over_represented_pvalue","under_represented_pvalue",
                          "numDEInCat","numInCat","term","ontology","qvalue")
    if (is.null(adjust)) {
        godata_interesting <- subset(godata, godata$over_represented_pvalue < pvalue)
        padjust_method <- "none"
    } else {  ## There is a requested pvalue adjustment
        godata_interesting <- subset(godata, p.adjust(godata[["over_represented_pvalue"]], method=padjust_method) <= adjust)
        if (dim(godata_interesting)[1] < arglist[["minimum_interesting"]]) {
            message(paste("simple_goseq(): There are no genes with an adj.p<", adjust, " using: ", padjust_method, ".", sep=""))
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
    mf_interesting <- mf_interesting[, c("ontology","numDEInCat","numInCat","over_represented_pvalue","qvalue","term")]
    ##bp_interesting <- subset(godata_interesting, ontology == "BP")
    bp_interesting <- godata_interesting[godata_interesting[["ontology"]] == "BP", ]
    rownames(bp_interesting) <- bp_interesting[["category"]]
    bp_interesting <- bp_interesting[ ,c("ontology","numDEInCat","numInCat","over_represented_pvalue","qvalue","term")]
    ##cc_interesting <- subset(godata_interesting, ontology == "CC")
    cc_interesting <- godata_interesting[godata_interesting[["ontology"]] == "CC", ]
    rownames(cc_interesting) <- cc_interesting[["category"]]
    cc_interesting <- cc_interesting[, c("ontology","numDEInCat","numInCat","over_represented_pvalue","qvalue","term")]

    pval_plots <- list(
        "bpp_plot_over" = pvalue_plots[["bpp_plot_over"]],
        "mfp_plot_over" = pvalue_plots[["mfp_plot_over"]],
        "ccp_plot_over" = pvalue_plots[["ccp_plot_over"]])

    return_list <- list("input" = de_genes,
                        "pwf" = pwf,
                        "pwf_plot" = pwf_plot,
                        "alldata" = godata,
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

#' Make a pvalue plot from goseq data.
#'
#' With minor changes, it is possible to push the goseq results into a clusterProfiler-ish pvalue
#' plot.  This handles those changes and returns the ggplot results.
#'
#' @param goterms Some data from goseq!
#' @param wrapped_width Number of characters before wrapping to help legibility.
#' @param cutoff Pvalue cutoff for the plot.
#' @param n How many groups to include?
#' @param mincat Minimum size of the category for inclusion.
#' @param level Levels of the ontology tree to use.
#' @return Plots!
#' @seealso \link[goseq]{goseq} \pkg{clusterProfiler} \code{\link{plot_ontpval}}
#' @export
plot_goseq_pval <- function(goterms, wrapped_width=20, cutoff=0.1, n=10, mincat=10, level=NULL) {
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
    plotting_mf[["score"]] <- plotting_mf[["numDEInCat"]] / plotting_mf[["numInCat"]]
    plotting_mf <- plotting_mf[ plotting_mf[["ontology"]] == "MF", ]
    plotting_mf <- plotting_mf[ plotting_mf[["term"]] != "NULL", ]
    plotting_mf <- plotting_mf[ plotting_mf[["over_represented_pvalue"]] <= cutoff, ]
    plotting_mf <- plotting_mf[ plotting_mf[["numInCat"]] > mincat, ]
    plotting_mf <- plotting_mf[order(plotting_mf[["over_represented_pvalue"]]),]
    plotting_mf <- head(plotting_mf, n=n)
    plotting_mf <- plotting_mf[, c("term","over_represented_pvalue","score")]
    plotting_mf[["term"]] <- as.character(lapply(strwrap(plotting_mf[["term"]], wrapped_width, simplify=FALSE), paste, collapse="\n"))
    colnames(plotting_mf) <- c("term","pvalue","score")
    mf_pval_plot <- plot_ontpval(plotting_mf, ontology="MF")

    plotting_bp <- subset(goterms, complete.cases(goterms))
    plotting_bp[["score"]] <- plotting_bp[["numDEInCat"]] / plotting_bp[["numInCat"]]
    plotting_bp <- plotting_bp[ plotting_bp[["ontology"]] == "BP", ]
    plotting_bp <- plotting_bp[ plotting_bp[["term"]] != "NULL", ]
    plotting_bp <- plotting_bp[ plotting_bp[["over_represented_pvalue"]] <= cutoff, ]
    plotting_bp <- plotting_bp[ plotting_bp[["numInCat"]] > mincat, ]
    plotting_bp <- plotting_bp[order(plotting_bp[["over_represented_pvalue"]]),]
    plotting_bp <- head(plotting_bp, n=n)
    plotting_bp <- plotting_bp[, c("term","over_represented_pvalue","score")]
    colnames(plotting_bp) <- c("term","pvalue","score")
    plotting_bp[["term"]] <- as.character(lapply(strwrap(plotting_bp[["term"]], wrapped_width, simplify=FALSE), paste, collapse="\n"))
    bp_pval_plot <- plot_ontpval(plotting_bp, ontology="BP")

    plotting_cc <- subset(goterms, complete.cases(goterms))
    plotting_cc[["score"]] <- plotting_cc[["numDEInCat"]] / plotting_cc[["numInCat"]]
    plotting_cc <- plotting_cc[ plotting_cc[["ontology"]] == "CC", ]
    plotting_cc <- plotting_cc[ plotting_cc[["term"]] != "NULL", ]
    plotting_cc <- plotting_cc[ plotting_cc[["over_represented_pvalue"]] <= cutoff, ]
    plotting_cc <- plotting_cc[ plotting_cc[["numInCat"]] > mincat, ]
    plotting_cc <- plotting_cc[order(plotting_cc[["over_represented_pvalue"]]),]
    plotting_cc <- head(plotting_cc, n=n)
    plotting_cc <- plotting_cc[, c("term","over_represented_pvalue","score")]
    colnames(plotting_cc) <- c("term","pvalue","score")
    plotting_cc[["term"]] <- as.character(lapply(strwrap(plotting_cc[["term"]], wrapped_width, simplify=FALSE), paste, collapse="\n"))
    cc_pval_plot <- plot_ontpval(plotting_cc, ontology="CC")

    pval_plots <- list(
        "mfp_plot_over" = mf_pval_plot,
        "bpp_plot_over" = bp_pval_plot,
        "ccp_plot_over" = cc_pval_plot,
        "mf_subset_over" = plotting_mf,
        "bp_subset_over" = plotting_bp,
        "cc_subset_over" = plotting_cc)
    return(pval_plots)
}

#' Make fun trees a la topgo from goseq data.
#'
#' This seeks to force goseq data into a format suitable for topGO and then use its tree plotting
#' function to make it possible to see significantly increased ontology trees.
#'
#' @param de_genes Some differentially expressed genes.
#' @param godata Data from goseq.
#' @param goid_map File to save go id mapping.
#' @param score_limit Score limit for the coloring.
#' @param goids_df Mapping of IDs to GO in the Ramigo expected format.
#' @param overwrite Overwrite the trees?
#' @param selector Function for choosing genes.
#' @param pval_column Column to acquire pvalues.
#' @return A plot!
#' @seealso \pkg{Ramigo}
#' @export
goseq_trees <- function(de_genes, godata, goid_map="reference/go/id2go.map",
                        score_limit=0.01, goids_df=NULL, overwrite=FALSE,
                        selector="topDiffGenes", pval_column="adj.P.Val") {
    mapping <- make_id2gomap(goid_map=goid_map, goids_df=goids_df, overwrite=overwrite)
    geneID2GO <- topGO::readMappings(file=goid_map)
    annotated_genes <- names(geneID2GO)
    if (is.null(de_genes[["ID"]])) {
        de_genes[["ID"]] <- make.names(rownames(de_genes), unique=TRUE)
    }
    interesting_genes <- factor(annotated_genes %in% de_genes[["ID"]])
    names(interesting_genes) <- annotated_genes

    if (is.null(de_genes[[pval_column]])) {
        mf_GOdata <- new("topGOdata", ontology="MF", allGenes=interesting_genes, annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata <- new("topGOdata", ontology="BP", allGenes=interesting_genes, annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata <- new("topGOdata", ontology="CC", allGenes=interesting_genes, annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
    } else {
        pvals <- as.vector(as.numeric(de_genes[[pval_column]]))
        names(pvals) <- rownames(de_genes)
        requireNamespace("topGO")
        attachNamespace("topGO")
        mf_GOdata <- new("topGOdata", description="MF", ontology="MF", allGenes=pvals,
                         geneSel=get(selector), annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        bp_GOdata <- new("topGOdata", description="BP", ontology="BP", allGenes=pvals,
                         geneSel=get(selector), annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
        cc_GOdata <- new("topGOdata", description="CC", ontology="CC", allGenes=pvals,
                         geneSel=get(selector), annot=topGO::annFUN.gene2GO, gene2GO=geneID2GO)
    }

    enriched_ids <- godata$alldata[["category"]]
    enriched_scores <- godata$alldata[["over_represented_pvalue"]]
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
    trees <- list(
        "MF_over" = mf_tree,
        "BP_over" = bp_tree,
        "CC_over" = cc_tree,
        "MF_overdata" = mf_tree_data,
        "BP_overdata" = bp_tree_data,
        "CC_overdata" = cc_tree_data)
    return(trees)
}

## EOF
