start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("060expt.R:
  1234567890123456789012345678901234567890123456789012345678901234567\n")
## Functions exported in expt.r
## make_pombe_expt(), create_expt(), concatenate_runs(), exclude_genes_expt(),
## features_greater_than(), make_exampledata(), median_by_factor(),
## set_expt_batches(), set_expt_colors(), set_expt_conditions(),
## set_expt_factors(), set_expt_samplenames(), subset_expt(), what_happened(),
## write_expt()
## S4 methods: exprs(), fData(), pData(), notes()

## make_pombe_expt() invokes create_expt()
pombe_expt <- make_pombe_expt()
chosen_genes <- c("SPAC212.11", "SPAC212.09c", "SPNCRNA.70",
                  "SPAC212.12", "SPAC212.04c", "SPAC212.03",
                  "SPAC212.02", "SPAC212.01c", "SPAC977.03")

## fData()
testing <- fData(pombe_expt)
actual <- dim(testing)
expected <- c(7039, 10)
## 0102
test_that("Do we get annotation data from our expt?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})
actual <- testing[chosen_genes, "start_position"]
expected <- c("1", "7619", "11027", "15855", "21381",
              "23589", "27353", "28738", "33835")
## 03
test_that("Do we get annotation data from our expt?", {
  expect_equal(actual, expected)
})

## pData()
testing <- pData(pombe_expt)
actual <- dim(testing)
expected <- c(36, 8)
## 0405
test_that("Do we get experimental metadata from our expt?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})

## exprs()
testing <- exprs(pombe_expt)
actual <- dim(testing)
expected <- c(7039, 36)
## 0607
test_that("Do we get expression from our expt?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})
## Check that the number distribution is what we expect
actual <- as.numeric(summary(testing[, 1]))
expected <- c(0.000, 46.000, 306.000, 2225.509, 869.000, 4542884.000)
## 13
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
## 1415
test_that("Do we get a reasonable number of resulting samples if we collapse by time?", {
  expect_equal(actual[1], expected[1], tolerance=0.001)
  expect_equal(actual[2], expected[2], tolerance=0.001)
})

## exclude_genes_expt()
no_rrna <- exclude_genes_expt(pombe_expt,
                              column="gene_biotype",
                              patterns=c("ncRNA", "pseudogene"))
expected <- 5778
actual <- nrow(exprs(no_rrna))
## 16
test_that("Does exclude_genes_expt remove stuff?", {
  expect_equal(actual, expected)
})

## features_greater_than()
testing_1 <- features_greater_than(pombe_expt, cutoff=1, hard=FALSE)
testing_10 <- features_greater_than(pombe_expt, cutoff=10, hard=FALSE)
testing_100 <- features_greater_than(pombe_expt, cutoff=100, hard=FALSE)
testing_1000 <- features_greater_than(pombe_expt, cutoff=1000, hard=FALSE)
testing_10000 <- features_greater_than(pombe_expt, cutoff=10000, hard=FALSE)
testing_100000 <- features_greater_than(pombe_expt, cutoff=100000, hard=FALSE)
## 22
test_that("Do we get expected features greater than some cutoffs?", {
  expect_equal(as.numeric(testing_1[["number"]][1]), 6457)
  expect_equal(as.numeric(testing_10[["number"]][1]), 5974)
  expect_equal(as.numeric(testing_100[["number"]][1]), 4754)
  expect_equal(as.numeric(testing_1000[["number"]][1]), 1566)
  expect_equal(as.numeric(testing_10000[["number"]][1]), 212)
  expect_equal(as.numeric(testing_100000[["number"]][1]), 5)
})

## make_exampledata()
testing <- DESeq2::counts(make_exampledata())
actual <- dim(testing)
expected <- c(1000, 5)
## 2324
test_that("Do we get some example data from make_exampledata()?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})

## median_by_factor()
medians <- median_by_factor(pombe_expt)
actual <- medians[chosen_genes, 1]
expected <- c(8, 31, 0, 0, 21, 0, 0, 2, 0)
## 25
test_that("Do we get expected medians?", {
  expect_equal(expected, actual)
})

## set_expt_batches()
new_batches <- c(rep(x="a", times=12), rep(x="b", times=12), rep(x="c", times=12))
testing <- set_expt_batches(pombe_expt, fact=new_batches)
actual <- pData(testing)[["batch"]]
## 26
test_that("Did we change the batches?", {
  expect_equal(new_batches, actual)
})

## set_expt_colors()
testing <- set_expt_colors(pombe_expt, chosen_palette="Accent")
old_actual <- levels(as.factor(pombe_expt$colors))
old_expected <- c("#1B9E77", "#666666", "#7D8F31", "#7E6EA2", "#8E7037", "#93752C",
                  "#A0A811", "#B3499C", "#B78415", "#BD6332", "#CF3F76", "#E0A604")
new_actual <- levels(as.factor(testing$colors))
new_expected <- c("#5C86AB", "#666666", "#7FC97F", "#8B3B99", "#9E5F33", "#A7B7B5",
                  "#CC4233", "#CFB2BE", "#DAE49D", "#EB0A75", "#F7BE8D", "#FEE290")
## 2728
test_that("Did we get some old/new colors?", {
  expect_equal(old_actual, old_expected)
  expect_equal(new_actual, new_expected)
})

## set_expt_conditions()
testing <- set_expt_conditions(pombe_expt, fact="minute")
expected <- levels(pombe_expt[["design"]][["minute"]])
actual <- levels(testing[["design"]][["condition"]])
## 29
test_that("Did we get some new conditions?", {
  expect_equal(actual, expected)
})

## set_expt_factors()
## This function just calls set_expt_batches and set_expt_conditions.  I am not testing it.

## set_expt_samplenames()
newnames <- paste0(sampleNames(pombe_expt), "test")
testing <- set_expt_samplenames(pombe_expt, newnames=newnames)
## 30
test_that("Can we change the same names?", {
  expect_equal(newnames, sampleNames(testing))
})

## subset_expt()
testing <- subset_expt(pombe_expt, subset="replicate=='r1'")
actual <- sampleNames(testing)
expected <- c("GSM1368273", "GSM1368276", "GSM1368279", "GSM1368282",
              "GSM1368285", "GSM1368288", "GSM1368291", "GSM1368294",
              "GSM1368297", "GSM1368300", "GSM1368303", "GSM1368306")
## 31
test_that("Can we subset an expt?", {
  expect_equal(expected, actual)
})

## what_happened()
testing <- normalize_expt(pombe_expt, transform="log2", norm="tmm",
                          convert="cpm", filter=TRUE, batch="sva")
actual <- what_happened(testing)
expected <- glue::glue("log2(sva(cpm(tmm(hpgl(data)))))")
## 32
test_that("Will an expt tell us what happened to it?", {
  expect_equal(expected, actual)
})

## write_expt()
testing <- write_expt(pombe_expt, excel="testing_write_expt.xlsx")
## 33
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
  ## 50
  expect_equal("gg", class(testing[["raw_scaled_pca"]])[[1]])
  expect_equal("data.frame", class(testing[["raw_scaled_pca_table"]])[[1]])
  expect_equal("gg", class(testing[["raw_scaled_tsne"]])[[1]])
  expect_equal("data.frame", class(testing[["raw_scaled_tsne_table"]])[[1]])
  expect_equal("recordedplot", class(testing[["raw_qq"]])[[1]])
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
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 060expt.R in ", elapsed,  " seconds."))
