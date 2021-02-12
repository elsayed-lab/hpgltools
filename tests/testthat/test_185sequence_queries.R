start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("185sequence_queries.R:
  123456789012\n")

a909_fasta <- system.file("gbs_tnseq/sagalactiae_a909.fasta", package = "hpgltools")
number_atg <- count_nmer(a909_fasta)
expected <- 37134
actual <- as.numeric(number_atg)

test_that("Can we count patterns in sequence data?", {
  expect_equal(expected, actual)
})

dm_annot <- new.env()
load("pasilla.rda", envir = dm_annot)
annot_df <- fData(dm_annot$expt)
bsg <- "BSgenome.Dmelanogaster.UCSC.dm6"
bsg_installed <- please_install(bsg)
lib_result <- sm(requireNamespace(bsg))
att_result <- sm(try(attachNamespace(bsg), silent = TRUE))
bsg <- get0(bsg)
annot_df[["chromosome_name"]] <- paste0("chr", annot_df[["chromosome_name"]])

dm_utrs <- gather_utrs_padding(bsg, annot_df = annot_df, name_column = "ensembl_gene_id",
                               chr_column = "chromosome_name", start_column = "start_position",
                               end_column = "end_position", strand_column = "strand",
                               type_column = "gene_biotype", gene_type = "protein_coding", padding = 100)
testthat("Do we get some utrs?", {
  expect_equal(ncol(dm_utrs), 15)
  expect_equal(nrow(dm_utrs), 6973)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 185sequence_queries.R in ", elapsed,  " seconds."))
