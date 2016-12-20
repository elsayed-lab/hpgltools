#' Simplification function for gostats, in the same vein as those written for clusterProfiler,
#' goseq, and topGO.
#'
#' GOstats has a couple interesting peculiarities:  Chief among them: the gene IDs must be
#' integers. As a result, I am going to have this function take a gff file in order to get the go
#' ids and gene ids on the same page.
#'
#' @param sig_genes Input list of differentially expressed genes.
#' @param gff Annotation information for this genome.
#' @param goids_df Set of GOids, as before in the format ID/GO.
#' @param universe_merge Column from which to create the universe of genes.
#' @param second_merge_try If the first universe merge fails, try this.
#' @param species Genbank organism to use.
#' @param pcutoff Pvalue cutoff for deciding significant.
#' @param direction Under or over represented categories.
#' @param conditional Perform a conditional search?
#' @param categorysize Category size below which to not include groups.
#' @param gff_type Gff column to use for creating the universe.
#' @param ... More parameters!
#' @return List of returns from GSEABase, Category, etc.
#' @seealso \pkg{GSEABase} \pkg{Category}
#' @export
simple_gostats <- function(sig_genes, gff, goids_df, universe_merge="id", second_merge_try="locus_tag",
                           species="fun", pcutoff=0.10, direction="over", conditional=FALSE,
                           categorysize=NULL, gff_type="cds", ...) {
    ## The import(gff) is being used for this primarily because it uses integers for the rownames and because it (should)
    ## contain every gene in the 'universe' used by GOstats, as much it ought to be pretty much
    ## perfect.
    arglist = list(...)
    if (!is.null(arglist[["gff_df"]])) {
        annotation <- arglist[["gff_df"]]
    } else {
        annotation <- gff2df(gff, type=gff_type)
    }
    colnames(annotation) <- tolower(colnames(annotation))
    colnames(annotation) <- gsub(pattern="length", replacement="width", x=colnames(annotation))

    ## This is similar to logic in ontology_goseq and is similarly problematic.
    ## Some gff files I use have all the annotation data in the type called 'gene', others use 'CDS', others use 'exon'
    ## I need a robust method of finding the correct feature type to call upon.

    ## I think there might be a weird environment collision occuring which is causing
    ## some gostats functionality to fail when functions are called using Category:: explicitly.
    ## Therefore I am loading these environments here and calling the functions without ::
    ## For a further discussion of what is happening:
    ## https://stat.ethz.ch/pipermail/bioconductor/2009-November/030348.html
    try(detach("package:GOstats", unload=TRUE), silent=TRUE)
    try(detach("package:Category", unload=TRUE), silent=TRUE)
    ## In theory, requireNamespace is sufficient, but that is not true.
    requireNamespace("GOstats")
    requireNamespace("GSEABase")
    requireNamespace("AnnotationDbi")
    requireNamespace("Category")
    message("The namespaces/environments uses by GOstats are entirely too complex.")
    message("If I try to call functions with Category:: or GOstats:: then they collide")
    message("And things fail without error, but if I try library() then R CMD check")
    message("gets pissed, well I tried both ways and I am calling library().")
    message("R CMD check can bit my shiny metal ass.")
    library("GOstats")
    ## library("Category")
    ## library("GSEABase")
    ## library("GOstats")
    ## library("AnnotationDbi")
    message(paste0("simple_gostats(): gff_type is: ", gff_type, ". Change that if there are bad merges."))
    types <- c("cds","gene","exon","protein_coding")
    for (type in types) {
        message(paste0("simple_gostats(): type ", type, " has ", sum(annotation$type == type), " annotations."))
    }

    annotation <- annotation[annotation$type == gff_type, ]
    message(paste0("simple_gostats(): the current annotations has: ", nrow(annotation),
                   " rows and ", ncol(annotation), " columns."))
    annotation[, universe_merge] <- make.names(annotation[, universe_merge], unique=TRUE)
    if (universe_merge %in% names(annotation)) {
        universe <- annotation[, c(universe_merge, "width")]
    } else if (second_merge_try %in% names(annotation)) {
        universe <- annotation[, c(second_merge_try, "width")]
    } else if ("transcript_name" %in% names(annotation)) {
        universe <- annotation[, c("transcript_name", "width")]
    } else {
        stop("simple_gostats(): Unable to cross reference annotations into universe, perhaps change gff_type to make the merge work.")
    }
    ## This section is a little odd
    ## The goal is to collect a consistent set of numeric gene IDs
    ## In addition, one must cross reference those IDs consistently with the universe of all genes.
    ## Thus in a few linues I will be doing a merge of all genes against the sig_genes and another merge
    ## of the gene<->go mappings, finally extracting the portions of the resulting dataframe into a format suitable for
    ## casting as a GOFrame/GOAllFrame
    colnames(universe) <- c("geneid","width")
    universe$id <- rownames(universe)
    universe <- universe[complete.cases(universe), ]

    if (is.null(sig_genes$ID)) {
        sig_genes$ID <- rownames(sig_genes)
    }
    universe_cross_de <- merge(universe, sig_genes, by.x="geneid", by.y="ID")
    degenes_ids <- universe_cross_de$id
    universe_ids <- universe$id
    ## Sometimes I have the columns set to 'ID','GO' -- others I have 'ORF','GO'
    ## FIXME!  This should be standardized.
    colnames(goids_df) <- c("ID","GO")
    gostats_go <- merge(universe, goids_df, by.x="geneid", by.y="ID")
    if (nrow(gostats_go) == 0) {
        stop("simple_gostats(): The merging of the universe vs. goids failed.")
    }
    if (ncol(gostats_go) == 5) {
        colnames(gostats_go) <- c("sysName","width", "frame.gene_id", "frame.go_id", "ID")
    } else if (ncol(gostats_go) == 4) {
        colnames(gostats_go) <- c("sysName","width", "frame.gene_id", "frame.go_id")
    } else {
        stop("Cannot set the columns for the gostats df.")
    }
    gostats_go[["frame.Evidence"]] <- "TAS"
    gostats_go <- gostats_go[, c("frame.go_id", "frame.Evidence", "frame.gene_id")]
    gostats_frame <- AnnotationDbi::GOFrame(gostats_go, organism=species)
    gostats_all <- AnnotationDbi::GOAllFrame(gostats_frame)
    message("simple_gostats(): Creating the gene set collection.  This is slow.")
    gsc <- GSEABase::GeneSetCollection(gostats_all,
                                       setType=GSEABase::GOCollection())

    mf_over <- bp_over <- cc_over <- NULL
    mf_under <- bp_under <- cc_under <- NULL
    message("simple_gostats(): Performing MF GSEA.")
    mf_params <- Category::GSEAGOHyperGParams(
        name=paste("GSEA of ", species, sep=""), geneSetCollection=gsc,
        geneIds=degenes_ids, universeGeneIds=universe_ids,
        ontology="MF", pvalueCutoff=pcutoff,
        conditional=conditional, testDirection="over")
    ## This is where it fell over
    mf_over <- Category::hyperGTest(mf_params)
    message(paste0("Found ", nrow(GOstats::summary(mf_over)), " over MF categories."))
    message("simple_gostats(): Performing BP GSEA.")
    bp_params <- Category::GSEAGOHyperGParams(
        name=paste("GSEA of ", species, sep=""),
        geneSetCollection=gsc, geneIds=degenes_ids, universeGeneIds=universe_ids,
        ontology="BP", pvalueCutoff=pcutoff,
        conditional=FALSE, testDirection="over")
    bp_over <- Category::hyperGTest(bp_params)
    message(paste0("Found ", nrow(GOstats::summary(bp_over)), " over BP categories."))
    message("simple_gostats(): Performing CC GSEA.")
    cc_params <- Category::GSEAGOHyperGParams(
        name=paste("GSEA of ", species, sep=""), geneSetCollection=gsc,
        geneIds=degenes_ids, universeGeneIds=universe_ids,
        ontology="CC", pvalueCutoff=pcutoff,
        conditional=FALSE, testDirection="over")
    cc_over <- Category::hyperGTest(cc_params)
    message(paste0("Found ", nrow(GOstats::summary(cc_over)), " over CC categories."))
    message("simple_gostats(): Performing under MF GSEA.")
    mf_params <- Category::GSEAGOHyperGParams(
        name=paste("GSEA of ", species, sep=""), geneSetCollection=gsc,
        geneIds=degenes_ids, universeGeneIds=universe_ids,
        ontology="MF", pvalueCutoff=pcutoff,
        conditional=conditional, testDirection="under")
    mf_under <- Category::hyperGTest(mf_params)
    message(paste0("Found ", nrow(GOstats::summary(mf_under)), " under MF categories."))
    message("simple_gostats(): Performing under BP GSEA.")
    bp_params <- Category::GSEAGOHyperGParams(
        name=paste("GSEA of ", species, sep=""), geneSetCollection=gsc,
        geneIds=degenes_ids, universeGeneIds=universe_ids,
        ontology="BP", pvalueCutoff=pcutoff,
        conditional=FALSE, testDirection="under")
    bp_under <- Category::hyperGTest(bp_params)
    message(paste0("Found ", nrow(GOstats::summary(bp_under)), " under BP categories."))
    message("simple_gostats(): Performing under CC GSEA.")
    cc_params <- Category::GSEAGOHyperGParams(
        name=paste("GSEA of ", species, sep=""), geneSetCollection=gsc,
        geneIds=degenes_ids, universeGeneIds=universe_ids,
        ontology="CC", pvalueCutoff=pcutoff,
        conditional=FALSE, testDirection="under")
    cc_under <- Category::hyperGTest(cc_params)
    message(paste0("Found ", nrow(GOstats::summary(cc_under)), " under CC categories."))
    mf_over_table <- bp_over_table <- cc_over_table <- NULL
    mf_under_table <- bp_under_table <- cc_under_table <- NULL
    mf_over_table <- GOstats::summary(mf_over, pvalue=1.0, htmlLinks=TRUE)
    bp_over_table <- GOstats::summary(bp_over, pvalue=1.0, htmlLinks=TRUE)
    cc_over_table <- GOstats::summary(cc_over, pvalue=1.0, htmlLinks=TRUE)
    mf_under_table <- GOstats::summary(mf_under, pvalue=1.0, htmlLinks=TRUE)
    bp_under_table <- GOstats::summary(bp_under, pvalue=1.0, htmlLinks=TRUE)
    cc_under_table <- GOstats::summary(cc_under, pvalue=1.0, htmlLinks=TRUE)
    if (!is.null(dim(mf_over_table))) {
        mf_over_table$qvalue <- tryCatch(
        {
            ttmp <- as.numeric(mf_over_table$Pvalue)
            ttmp <- qvalue::qvalue(ttmp, robust=TRUE)$qvalues
            signif(x=ttmp, digits=4)
        },
        error=function(cond) {
            return(1)
        },
        finally={
        })
    }
    if (!is.null(dim(bp_over_table))) {
        bp_over_table$qvalue <- tryCatch(
        {
            ttmp <- as.numeric(bp_over_table$Pvalue)
            ttmp <- qvalue::qvalue(ttmp, robust=TRUE)$qvalues
            signif(x=ttmp, digits=4)
        },
        error=function(cond) {
            return(1)
        },
        finally={
        })
    }
    if (!is.null(dim(cc_over_table))) {
        cc_over_table$qvalue <- tryCatch(
        {
            ttmp <- as.numeric(cc_over_table$Pvalue)
            ttmp <- qvalue::qvalue(ttmp, robust=TRUE)$qvalues
            signif(x=ttmp, digits=4)
        },
        error=function(cond) {
            return(1)
        },
        finally={
        })
    }
    if (!is.null(dim(mf_under_table))) {
        mf_under_table$qvalue <- tryCatch(
        {
            ttmp <- as.numeric(mf_under_table$Pvalue)
            ttmp <- qvalue::qvalue(ttmp, robust=TRUE)$qvalues
            signif(x=ttmp, digits=4)
        },
        error=function(cond) {
            return(1)
        },
        finally={
        })
    }
    if (!is.null(dim(bp_under_table))) {
        bp_under_table$qvalue <- tryCatch(
        {
            ttmp <- as.numeric(bp_under_table$Pvalue)
            ttmp <- qvalue::qvalue(ttmp, robust=TRUE)$qvalues
            signif(x=ttmp, digits=4)
        },
        error=function(cond) {
            return(1)
        },
        finally={
        })
    }
    if (!is.null(dim(cc_under_table))) {
        cc_under_table$qvalue <- tryCatch(
        {
            ttmp <- as.numeric(cc_under_table$Pvalue)
            ttmp <- qvalue::qvalue(ttmp, robust=TRUE)$qvalues
            signif(x=ttmp, digits=4)
        },
        error=function(cond) {
            return(1)
        },
        finally={
        })
    }

    if (is.null(categorysize)) {
        mf_over_sig <- GOstats::summary(mf_over)
        bp_over_sig <- GOstats::summary(bp_over)
        cc_over_sig <- GOstats::summary(cc_over)
        mf_under_sig <- GOstats::summary(mf_under)
        bp_under_sig <- GOstats::summary(bp_under)
        cc_under_sig <- GOstats::summary(cc_under)
    } else {
        mf_over_sig <- GOstats::summary(mf_over, categorySize=categorysize)
        bp_over_sig <- GOstats::summary(bp_over, categorySize=categorysize)
        cc_over_sig <- GOstats::summary(cc_over, categorySize=categorysize)
        mf_under_sig <- GOstats::summary(mf_under, categorySize=categorysize)
        bp_under_sig <- GOstats::summary(bp_under, categorySize=categorysize)
        cc_under_sig <- GOstats::summary(cc_under, categorySize=categorysize)
    }
    if (!is.null(dim(mf_over_sig))) {
        mf_over_sig$definition <- try(godef(mf_over_sig$GOMFID), silent=TRUE)
    } else {
        mf_over_sig <- NULL
    }
    if (!is.null(dim(bp_over_sig))) {
        bp_over_sig$definition <- try(godef(bp_over_sig$GOBPID), silent=TRUE)
    } else {
        bp_over_sig <- NULL
    }
    if (!is.null(dim(cc_over_sig))) {
        cc_over_sig$definition <- try(godef(cc_over_sig$GOCCID), silent=TRUE)
    } else {
        bp_over_sig <- NULL
    }
    if (!is.null(dim(mf_under_sig))) {
        mf_under_sig$definition <- try(godef(mf_under_sig$GOMFID), silent=TRUE)
    } else {
        mf_under_sig <- NULL
    }
    if (!is.null(dim(bp_under_sig))) {
        bp_under_sig$definition <- try(godef(bp_under_sig$GOBPID), silent=TRUE)
    } else {
        bp_under_sig <- NULL
    }
    if (!is.null(dim(cc_under_sig))) {
        cc_under_sig$definition <- try(godef(cc_under_sig$GOCCID), silent=TRUE)
    } else {
        bp_under_sig <- NULL
    }

    gostats_p_mf_over <- try(plot_histogram(mf_over_table$Pvalue, bins=20), silent=TRUE)
    gostats_p_mf_under <- try(plot_histogram(mf_under_table$Pvalue, bins=20), silent=TRUE)
    gostats_p_bp_over <- try(plot_histogram(bp_over_table$Pvalue, bins=20), silent=TRUE)
    gostats_p_bp_under <- try(plot_histogram(bp_under_table$Pvalue, bins=20), silent=TRUE)
    gostats_p_cc_over <- try(plot_histogram(cc_over_table$Pvalue, bins=20), silent=TRUE)
    gostats_p_cc_under <- try(plot_histogram(cc_under_table$Pvalue, bins=20), silent=TRUE)

    ret_list <- list(
        "mf_over_all" = mf_over_table,
        "bp_over_all" = bp_over_table,
        "cc_over_all" = cc_over_table,
        "mf_under_all" = mf_under_table,
        "bp_under_all" = bp_under_table,
        "cc_under_all" = cc_under_table,
        "mf_over_enriched" = mf_over_sig,
        "bp_over_enriched" = bp_over_sig,
        "cc_over_enriched" = cc_over_sig,
        "mf_under_enriched" = mf_under_sig,
        "bp_under_enriched" = bp_under_sig,
        "cc_under_enriched" = cc_under_sig,
        "gostats_mfp_over" = gostats_p_mf_over,
        "gostats_bpp_over" = gostats_p_bp_over,
        "gostats_ccp_over" = gostats_p_cc_over,
        "gostats_mfp_under" = gostats_p_mf_under,
        "gostats_bpp_under" = gostats_p_bp_under,
        "gostats_ccp_under" = gostats_p_cc_under)

    pvalue_plots <- try(plot_gostats_pval(ret_list))
    ret_list[["pvalue_plots"]] <- pvalue_plots
    return(ret_list)
}

## EOF
