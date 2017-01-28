start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("99cleanup.R: Cleaning up this mess.\n")

de_removed <- try(file.remove("de_deseq.rda"))
ed_removed <- try(file.remove("de_edger.rda"))
li_removed <- try(file.remove("de_limma.rda"))
ba_removed <- try(file.remove("de_basic.rda"))

limma_combat_removed <- try(file.remove("de_limma_combat.rda"))
pasilla_removed <- try(file.remove("pasilla.Rdata"))
plots_removed <- try(file.remove("Rplots.pdf"))
map_removed <- try(file.remove("id2go.map"))
gff_removed <- try(file.remove("dmel.gff"))
gb_removed <- try(file.remove("AE009949.gb"))

basic_removed <- try(file.remove("test_basic.xlsx"))
deseq_removed <- try(file.remove("test_deseq.xlsx"))
edger_removed <- try(file.remove("test_edger.xlsx"))
limma_removed <- try(file.remove("test_limma.xlsx"))

xlsx_table_removed <- try(file.remove("test_excel.xlsx"))
xlsx_sig_semoved <- try(file.remove("test_excel_sig.xlsx"))

remove_directories <- c("organdb", "pathview", "pathview_in", "circos", "test_excel", "test_excel_sig")
for (dir in remove_directories) {
    if (file.exists(dir)) {
        directory_removed <- try(unlink(dir, recursive=TRUE))
    }
}

## It annoys me greatly when something completes and doesn't add a newline, perhaps this will fix
## that.  Also, I would like to have something which counts the tests run.

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 99cleanup.R in ", elapsed,  " seconds."))
