start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("99cleanup.R: Cleaning up this mess.\n")

de_removed <- file.remove("de_deseq.rda")
ed_removed <- file.remove("de_edger.rda")
li_removed <- file.remove("de_limma.rda")
ba_removed <- file.remove("de_basic.rda")

limma_combat_removed <- file.remove("de_limma_combat.rda")
pasilla_removed <- file.remove("pasilla.Rdata")
plots_removed <- file.remove("Rplots.pdf")
map_removed <- file.remove("id2go.map")
gff_removed <- file.remove("dmel.gff")


remove_directories <- c("organdb", "pathview", "pathview_in", "circos", "test_excel", "test_excel_sig")
for (dir in remove_directories) {
    if (file.exists(dir)) {
        directory_removed <- unlink(dir, recursive=TRUE)
    }
}

## It annoys me greatly when something completes and doesn't add a newline, perhaps this will fix
## that.  Also, I would like to have something which counts the tests run.

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 99cleanup.R in ", elapsed,  " seconds."))
