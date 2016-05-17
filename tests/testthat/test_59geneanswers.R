library(testthat)
library(hpgltools)

context("Does seqAnswers work?")
## seqAnswers is new to me.  I have decided that henceforth all new functions will be implemented first as a test, then moved into the actual codebase.

## I want to do some much longer tests using goseq/clusterprofiler/topgo/gostats/gprofiler
## These will likely not work with Travis as they take forever.
## I am not however, certain how to use skip_on_travis(), so I printed it and am copying the
## useful bits here.
if (!identical(Sys.getenv("TRAVIS"), "true")) {

    limma <- new.env()
    load("de_limma.rda", envir=limma)
    table <- limma$hpgl_table
    sig_genes <- sp(get_sig_genes(table, column="untreated")$up_genes)$result

    tt <- sp(library("org.Dm.eg.db"))
    ## goseq_result <- simple_goseq(sig_genes, species="hsapiens")
    ## gprofiler_result <- suppressMessages(simple_gprofiler(sig_genes, species="dmelanogaster"))

    tt <- sp(require.auto("GeneAnswers"))
    tt <- sp(library("GeneAnswers"))
    tt <- sp(library("org.Dm.eg.db"))
    tt <- sp(library("GO.db"))

    ## get named vector of entrez ids
    fb.entrez <- unlist(as.list(org.Dm.egFLYBASE2EG))

    ## for a data frame x with flybase ids (column 1) and data values (column 2)
    ## match the flybase names against the vector of entrez ids
    sig_genes$ids = rownames(sig_genes)
    iv <- match(sig_genes$ids, names(fb.entrez))

    ## add a column for entrez ids
    sig_genes <- cbind(sig_genes, rep(NA, nrow(sig_genes)))

    ## fill it in by mapping the entez ids onto the matching flybase ids
    sig_genes$entrezid <- fb.entrez[iv]

    ## now you can do some GO analysis
    ## for an index vector "myTopHits" of your top data
    ##topset <- sig_genes[myTopHits,3]
    ## remove entries that had no matching entrez id(NA)
    ##topset <- topset[!is.na(topset)]

    ## Get BP enrichment
    ##foo <- geneAnswersBuilder(topset, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG')
    ##go.bp <- foo@enrichmentInfo

}
