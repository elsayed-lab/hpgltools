start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("50gsea_ontshared.R: Do the various ontology helper functions still work?\n")

expected_term_result <- "adenyl ribonucleotide binding"
actual_term_result <- sm(as.character(goterm()))
expected_syn_result <- "mitochondrial inheritance"
actual_syn_result <- sm(as.character(gosyn()))
expected_sec_result <- c("GO:0000141", "GO:0030482")
actual_sec_result <- sm(as.character(gosec()))
expected_def_result <- "An assembly of actin filaments that are on the same axis but may be oriented with the same or opposite polarities and may be packed with different levels of tightness."
actual_def_result <- sm(as.character(godef()))
expected_ont_result <- c("CC","CC")
actual_ont_result <- sm(as.character(goont()))
expected_level_result <- c("3", "3")
actual_level_result <- sm(as.character(golevel()))

test_that("Are GO.db functions working (term)?", {
    expect_equal(expected_term_result, actual_term_result)
})
test_that("Are GO.db functions working (synonym)?", {
    expect_equal(expected_syn_result, actual_syn_result)
})
test_that("Are GO.db functions working (secondary)?", {
    expect_equal(expected_sec_result, actual_sec_result)
})
test_that("Are GO.db functions working (definition)?", {
    expect_equal(expected_def_result, actual_def_result)
})
test_that("Are GO.db functions working (ontology)?", {
    expect_equal(expected_ont_result, actual_ont_result)
})
test_that("Are GO.db functions working (level)?", {
    expect_equal(expected_level_result, actual_level_result)
})

## I want to do some much longer tests using goseq/clusterprofiler/topgo/gostats/gprofiler
## These will likely not work with Travis as they take forever.
## I am not however, certain how to use skip_on_travis(), so I printed it and am copying the
## useful bits here.
## if (!identical(Sys.getenv("TRAVIS"), "true")) {
## }

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 50gsea_ontshared.R in ", elapsed,  " seconds."))
tt <- clear_session()
