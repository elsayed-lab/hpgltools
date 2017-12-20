start <- as.POSIXlt(Sys.time())
context("060expt.R: Test functions in expt.r")
## Functions exported in expt.r
## create_expt(), concatenate_runs(), exclude_genes_expt(),
## features_greater_than(), make_exampledata(), median_by_factor(),
## read_counts_expt(), set_expt_batch(), set_expt_colors(), set_expt_condition(),
## set_expt_factors(), set_expt_samplenames(), subset_expt(), what_happened(),
## write_expt()
## S4 methods: exprs(), fData(), pData(), notes()

tt <- sm(library(fission))
tt <- sm(data(fission))
meta <- as.data.frame(fission@colData)
meta$condition <- paste(meta$strain, meta$minute, sep=".")

meta$batch <- meta$replicate
meta$sample.id <- rownames(meta)
fission_data <- fission@assays$data$counts
fission_expt <- create_expt(metadata=meta, count_dataframe=fission_data)

## Neat, it works, and even figures out that the default mart is incorrect by itself.
pombe_annotations <- load_biomart_annotations(
  host="fungi.ensembl.org",
  trymart="fungal_mart",
  gene_requests=c("pombase_transcript", "ensembl_gene_id", "ensembl_transcript_id",
                  "hgnc_symbol", "description", "gene_biotype"),
  species="spombe", overwrite=TRUE)
pombe_mart <- pombe_annotations[["mart"]]
annotations <- pombe_annotations[["annotation"]]
rownames(annotations) <- make.names(gsub(pattern="\\.\\d+$", replacement="", x=rownames(annotations)), unique=TRUE)

## create_expt() with a dataframe.
pombe_expt <- create_expt(metadata=meta, count_dataframe=fission_data, gene_info=annotations)

## fData()
testing <- fData(pombe_expt)
actual <- dim(testing)
expected <- c(7039, 10)
test_that("Do we get annotation data from our expt?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})
actual <- head(testing$start_position)
expected <- c("1", "7619", "11027", "15855", "21381", "23589")
test_that("Do we get annotation data from our expt?", {
  expect_equal(actual, expected)
})

## pData()
testing <- pData(pombe_expt)
actual <- dim(testing)
expected <- c(36, 8)
test_that("Do we get experimental metadata from our expt?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})

## exprs()
testing <- exprs(pombe_expt)
actual <- dim(testing)
expected <- c(7039, 36)
test_that("Do we get expression from our expt?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})
## Check that the number distribution is what we expect
actual <- as.numeric(summary(testing[, 1]))
expected <- c(0.000, 46.000, 306.000, 2225.509, 869.000, 4542884.000)
test_that("Do we get expression from our expt?", {
  expect_equal(actual[1], expected[1], tolerance=0.001)
  expect_equal(actual[2], expected[2], tolerance=0.001)
  expect_equal(actual[3], expected[3], tolerance=0.001)
  expect_equal(actual[4], expected[4], tolerance=0.001)
  expect_equal(actual[5], expected[5], tolerance=0.001)
  expect_equal(actual[6], expected[6], tolerance=0.001)
})

## concatenate_runs()
test_expt <- concatenate_runs(expt=pombe_expt, column="minute")
actual <- dim(pData(test_expt))
expected <- c(6,8)
test_that("Do we get a reasonable number of resulting samples if we collapse by time?", {
  expect_equal(actual[1], expected[1], tolerance=0.001)
  expect_equal(actual[2], expected[2], tolerance=0.001)
})


end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 060expt.R in ", elapsed,  " seconds."))
