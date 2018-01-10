start <- as.POSIXlt(Sys.time())
context("060expt.R:\n")
## Functions exported in expt.r
## make_pombe_expt(), create_expt(), concatenate_runs(), exclude_genes_expt(),
## features_greater_than(), make_exampledata(), median_by_factor(),
## set_expt_batches(), set_expt_colors(), set_expt_conditions(),
## set_expt_factors(), set_expt_samplenames(), subset_expt(), what_happened(),
## write_expt()
## S4 methods: exprs(), fData(), pData(), notes()

## make_pombe_expt() invokes create_expt()
pombe_expt <- make_pombe_expt()

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

## exclude_genes_expt()
no_rrna <- sm(exclude_genes_expt(pombe_expt,
                                 column="gene_biotype",
                                 patterns=c("ncRNA", "pseudogene")))
expected <- 5778
actual <- nrow(exprs(no_rrna))
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
test_that("Do we get expected features greater than some cutoffs?", {
  expect_equal(as.numeric(testing_1[1]), 6457)
  expect_equal(as.numeric(testing_10[1]), 5974)
  expect_equal(as.numeric(testing_100[1]), 4754)
  expect_equal(as.numeric(testing_1000[1]), 1566)
  expect_equal(as.numeric(testing_10000[1]), 212)
  expect_equal(as.numeric(testing_100000[1]), 5)
})

## make_exampledata()
testing <- counts(make_exampledata())
actual <- dim(testing)
expected <- c(1000, 5)
test_that("Do we get some example data from make_exampledata()?", {
  expect_equal(actual[1], expected[1])
  expect_equal(actual[2], expected[2])
})

## median_by_factor()
medians <- sm(median_by_factor(pombe_expt))
actual <- head(medians)
row1 <- c(7, 2, 12, 6, 6, 1, 8, 3, 19, 10, 9, 8)
row2 <- c(41, 51, 34, 48, 61, 86, 31, 55, 33, 43, 75, 73)
row3 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
test_that("Do we get expected medians?", {
  expect_equal(row1, as.numeric(actual[1, ]))
  expect_equal(row2, as.numeric(actual[2, ]))
  expect_equal(row3, as.numeric(actual[3, ]))
})

## set_expt_batches()
new_batches <- c(rep(x="a", times=12), rep(x="b", times=12), rep(x="c", times=12))
testing <- set_expt_batches(pombe_expt, fact=new_batches)
actual <- pData(testing)[["batch"]]
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
test_that("Did we get some old/new colors?", {
  expect_equal(old_actual, old_expected)
  expect_equal(new_actual, new_expected)
})

## set_expt_conditions()
testing <- set_expt_conditions(pombe_expt, fact="minute")
expected <- levels(pombe_expt[["design"]][["minute"]])
actual <- levels(testing[["design"]][["condition"]])
test_that("Did we get some new conditions?", {
  expect_equal(actual, expected)
})

## set_expt_factors()
## This function just calls set_expt_batches and set_expt_conditions.  I am not testing it.

## set_expt_samplenames()
newnames <- paste0(sampleNames(pombe_expt), "test")
testing <- set_expt_samplenames(pombe_expt, newnames=newnames)
test_that("Can we change the same names?", {
  expect_equal(newnames, sampleNames(testing))
})

## subset_expt()
testing <- subset_expt(pombe_expt, subset="replicate=='r1'")
actual <- sampleNames(testing)
expected <- c("GSM1368273", "GSM1368276", "GSM1368279", "GSM1368282",
              "GSM1368285", "GSM1368288", "GSM1368291", "GSM1368294",
              "GSM1368297", "GSM1368300", "GSM1368303", "GSM1368306")
test_that("Can we subset an expt?", {
  expect_equal(expected, actual)
})

## what_happened()
testing <- sm(normalize_expt(pombe_expt, transform="log2", norm="quant",
                             convert="cpm", filter=TRUE, batch="sva"))
actual <- what_happened(testing)
expected <- "log2(sva(cpm(quant(hpgl(data)))))"
test_that("Will an expt tell us what happened to it?", {
  expect_equal(expected, actual)
})

## write_expt()
testing <- sm(write_expt(pombe_expt, excel="testing_write_expt.xlsx"))
test_that("Did write_expt() work?", {
  expect_true(file.exists("testing_write_expt.xlsx"))
})
file.remove("testing_write_expt.xlsx")

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 060expt.R in ", elapsed,  " seconds."))
