de_removed <- file.remove("de_deseq.rda")
ed_removed <- file.remove("de_edger.rda")
## li_removed <- file.remove("de_limma.rda")
ba_removed <- file.remove("de_basic.rda")

limma_combat_removed <- file.remove("de_limma_combat.rda")
## pasilla_removed <- file.remove("pasilla.Rdata")
plots_removed <- file.remove("Rplots.pdf")

## It annoys me greatly when something completes and doesn't add a newline, perhaps this will fix
## that.  Also, I would like to have something which counts the tests run.
message("Finished tests.")
