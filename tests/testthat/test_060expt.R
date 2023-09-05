start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("060expt.R")

## make_pombe_expt() invokes create_expt()
pombe_expt <- make_pombe_expt()
chosen_genes <- c("SPAC212.11", "SPAC212.09c", "SPNCRNA.70",
                  "SPAC212.12", "SPAC212.04c", "SPAC212.03",
                  "SPAC212.02", "SPAC212.01c", "SPAC977.03")

testing <- fData(pombe_expt)
actual <- dim(testing)
expected <- c(7039, 10)
test_that("Do we get annotation data from our expt?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})
actual <- testing[chosen_genes, "start_position"]
expected <- c("1", "7619", "11027", "15855", "21381",
              "23589", "27353", "28738", "33835")
test_that("Do we get annotation data from our expt?", {
  expect_equal(actual, expected)
})

## Take a moment to test my new subset and print
high_genes <- rowSums(exprs(pombe_expt)) > 100000
high_expt <- pombe_expt[high_genes, ]
high_expt

testing <- pData(pombe_expt)
actual <- dim(testing)
expected <- c(36, 8)
test_that("Do we get experimental metadata from our expt?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})

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
  expect_equal(actual[1], expected[1], tolerance = 0.001)
  expect_equal(actual[2], expected[2], tolerance = 0.001)
  expect_equal(actual[3], expected[3], tolerance = 0.001)
  expect_equal(actual[4], expected[4], tolerance = 0.001)
  expect_equal(actual[5], expected[5], tolerance = 0.001)
  expect_equal(actual[6], expected[6], tolerance = 0.001)
})

test_expt <- concatenate_runs(expt = pombe_expt, column = "minute")
actual <- dim(pData(test_expt))
expected <- c(6, 8)
test_that("Do we get a reasonable number of resulting samples if we collapse by time?", {
  expect_equal(actual[1], expected[1], tolerance = 0.001)
  expect_equal(actual[2], expected[2], tolerance = 0.001)
})

no_rrna <- subset_genes(pombe_expt,
                        column = "gene_biotype",
                        patterns = c("ncRNA", "pseudogene"))
expected <- 5778
actual <- nrow(exprs(no_rrna))
test_that("Does exclude_genes_expt remove stuff?", {
  expect_equal(actual, expected)
})

testing_1 <- features_greater_than(pombe_expt, cutoff = 1, hard = FALSE)
testing_10 <- features_greater_than(pombe_expt, cutoff = 10, hard = FALSE)
testing_100 <- features_greater_than(pombe_expt, cutoff = 100, hard = FALSE)
testing_1000 <- features_greater_than(pombe_expt, cutoff = 1000, hard = FALSE)
testing_10000 <- features_greater_than(pombe_expt, cutoff = 10000, hard = FALSE)
testing_100000 <- features_greater_than(pombe_expt, cutoff = 100000, hard = FALSE)
test_that("Do we get expected features greater than some cutoffs?", {
  expect_equal(as.numeric(testing_1[["number"]][1]), 6457)
  expect_equal(as.numeric(testing_10[["number"]][1]), 5974)
  expect_equal(as.numeric(testing_100[["number"]][1]), 4754)
  expect_equal(as.numeric(testing_1000[["number"]][1]), 1566)
  expect_equal(as.numeric(testing_10000[["number"]][1]), 212)
  expect_equal(as.numeric(testing_100000[["number"]][1]), 5)
})

## median_by_factor()
## I changed the output of this function!
medians <- median_by_factor(pombe_expt)[["medians"]]
actual <- medians[chosen_genes, 1]
expected <- c(8, 31, 0, 0, 21, 0, 0, 2, 0)
test_that("Do we get expected medians?", {
  expect_equal(expected, actual)
})

new_batches <- c(rep(x = "a", times = 12), rep(x = "b", times = 12), rep(x = "c", times = 12))
testing <- set_expt_batches(pombe_expt, fact = new_batches)
actual <- pData(testing)[["batch"]]
test_that("Did we change the batches?", {
  expect_equal(new_batches, actual)
})

testing <- set_expt_colors(pombe_expt, chosen_palette = "Accent")
old_actual <- levels(as.factor(pombe_expt$colors))
old_expected <- c("#1B9E77", "#666666", "#7D8F31", "#7E6EA2", "#8E7037", "#93752C",
                  "#A0A811", "#B3499C", "#B78415", "#BD6332", "#CF3F76", "#E0A604")
new_actual <- levels(as.factor(testing$colors))
new_expected <- c("#5C86AB", "#666666", "#7FC97F", "#8B3B99", "#9E5F33", "#A7B7B5",
                  "#CC4233", "#CFB2BE", "#DAE49D", "#EB0A75", "#F7BE8D", "#FEE290")
test_that("Did we get some old/new colors?", {
  expect_equal(old_actual, old_expected)
  expect_equal(new_actual, new_expected)
})

testing <- set_expt_conditions(pombe_expt, fact = "minute")
expected <- levels(pombe_expt[["design"]][["minute"]])
actual <- levels(testing[["design"]][["condition"]])
test_that("Did we get some new conditions?", {
  expect_equal(actual, expected)
})

newnames <- paste0(sampleNames(pombe_expt), "test")
testing <- set_expt_samplenames(pombe_expt, newnames = newnames)
test_that("Can we change the same names?", {
  expect_equal(newnames, sampleNames(testing))
})

testing <- subset_expt(pombe_expt, subset = "replicate=='r1'")
actual <- sampleNames(testing)
expected <- c("GSM1368273", "GSM1368276", "GSM1368279", "GSM1368282",
              "GSM1368285", "GSM1368288", "GSM1368291", "GSM1368294",
              "GSM1368297", "GSM1368300", "GSM1368303", "GSM1368306")
test_that("Can we subset an expt?", {
  expect_equal(expected, actual)
})

testing <- normalize_expt(pombe_expt, transform = "log2", norm = "tmm",
                          convert = "cpm", filter = TRUE, batch = "sva")
actual <- what_happened(testing)
expected <- glue::glue("log2(sva(cpm(tmm(cbcb(data)))))")
test_that("Will an expt tell us what happened to it?", {
  expect_equal(expected, actual)
})

## Using aes(size) to define line width is deprecated now; I will need to
## hunt for the function(s) doing this.
testing <- write_expt(pombe_expt, excel = "testing_write_expt.xlsx")
test_that("Did write_expt() work?", {
  expect_true(file.exists("testing_write_expt.xlsx"))
})
tt <- file.remove("testing_write_expt.xlsx")

## I want to test and ensure that the various plots are getting created, I have
## been bitten in the butt recently due to minor formatting changes.
test_that("Did write_expt() create pictures?", {
  expect_equal("data.frame", class(testing[["legend"]])[[1]])
  expect_equal("data.frame", class(testing[["annotations"]])[[1]])
  expect_equal("data.frame", class(testing[["design"]])[[1]])
  expect_equal("gg", class(testing[["raw_libsize"]])[[1]])
  expect_equal("gg", class(testing[["raw_nonzero"]])[[1]])
  expect_equal("gg", class(testing[["raw_density"]])[[1]])
  expect_equal("gg", class(testing[["raw_cv"]])[[1]])
  expect_equal("gg", class(testing[["raw_boxplot"]])[[1]])
  expect_equal("recordedplot", class(testing[["raw_corheat"]])[[1]])
  expect_equal("recordedplot", class(testing[["raw_disheat"]])[[1]])
  expect_equal("gg", class(testing[["raw_smc"]])[[1]])
  expect_equal("gg", class(testing[["raw_smd"]])[[1]])
  expect_equal("gg", class(testing[["raw_pca"]])[[1]])
  expect_equal("data.frame", class(testing[["raw_pca_table"]])[[1]])
  expect_equal("gg", class(testing[["raw_tsne"]])[[1]])
  expect_equal("data.frame", class(testing[["raw_tsne_table"]])[[1]])
  expect_equal("gg", class(testing[["raw_scaled_pca"]])[[1]])
  expect_equal("data.frame", class(testing[["raw_scaled_pca_table"]])[[1]])
  expect_equal("gg", class(testing[["raw_scaled_tsne"]])[[1]])
  expect_equal("data.frame", class(testing[["raw_scaled_tsne_table"]])[[1]])
  ## I turned off qq plots by default.
  ## expect_equal("recordedplot", class(testing[["raw_qq"]])[[1]])
  expect_equal("gg", class(testing[["norm_libsize"]])[[1]])
  expect_equal("gg", class(testing[["norm_nonzero"]])[[1]])
  expect_equal("gg", class(testing[["norm_density"]])[[1]])
  expect_equal("gg", class(testing[["norm_boxplot"]])[[1]])
  expect_equal("recordedplot", class(testing[["norm_corheat"]])[[1]])
  expect_equal("recordedplot", class(testing[["norm_disheat"]])[[1]])
  expect_equal("gg", class(testing[["norm_smc"]])[[1]])
  expect_equal("gg", class(testing[["norm_smd"]])[[1]])
  expect_equal("gg", class(testing[["norm_pca"]])[[1]])
  expect_equal("data.frame", class(testing[["norm_pca_table"]])[[1]])
  expect_equal("gg", class(testing[["norm_tsne"]])[[1]])
  expect_equal("data.frame", class(testing[["norm_tsne_table"]])[[1]])
  expect_equal("data.frame", class(testing[["medians"]])[[1]])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 060expt.R in ", elapsed,  " seconds.")
