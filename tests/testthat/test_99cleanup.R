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


directory_removed <- file.remove("organdb", recursive=TRUE)
directory_removed <- file.remove("pathview", recursive=TRUE)
directory_removed <- file.remove("pathview_in", recursive=TRUE)


## It annoys me greatly when something completes and doesn't add a newline, perhaps this will fix
## that.  Also, I would like to have something which counts the tests run.
test_that("Did this finish?", {
    expect_equal(1, 1)
})

message("\nFinished 99cleanup.R")
