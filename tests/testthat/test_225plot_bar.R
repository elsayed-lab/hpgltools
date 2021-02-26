start <- as.POSIXlt(Sys.time())
context("225plot_bar.R:
  1234\n")
## 2017-12, exported functions in plot_bar:
## plot_libsize, plot_libsize_prepost, plot_pct_kept, plot_significant_bar

pombe_expt <- make_pombe_expt(annotation = FALSE)

testing <- plot_libsize(pombe_expt)

actual_class <- class(testing[["plot"]])[[1]]
expected_class <- "gg"
actual_numbers <- testing[["table"]][["sum"]]
expected_numbers <- c(
  15665355, 9466415, 13701028, 18477590, 10352387, 9574669, 16337161, 10759148,
  13389388, 7013986, 17257211, 11113047, 12562371,  9761223, 10190649, 13047056,
  16554775, 10594463, 9331989, 12819822, 16405435, 10372996, 12060808, 10965123,
  18161629, 12115277, 12538806, 14260783, 12987616, 14589780, 13266122, 11760941,
  12202025, 12661612, 9968369,  9753854)
test_that("Do we get expected plot_libsize information?", {
  expect_equal(actual_class, expected_class)
  expect_equal(actual_numbers, expected_numbers)
})

## I think I might have something wrong here.
testing <- plot_libsize_prepost(pombe_expt)
actual_class <- class(testing[["count_plot"]])[[1]]
expected_class <- "gg"
actual_numbers <- testing[["table"]][["low"]]
expected_numbers <- c(
  782, 963, 839, 1060, 904, 1312, 845, 929, 893, 1006, 708, 789, 786, 906,
  928, 800, 822, 908, 908, 805, 771, 922, 834, 989, 790, 838, 771, 774,
  824, 782, 762, 795, 795, 858, 870, 949, 17, 38, 18, 113, 32, 259, 23,
  33, 33, 53, 4, 6, 11, 31, 29, 12, 19, 29, 35, 18, 12, 37, 13, 63, 13,
  17, 6, 10, 14, 7, 8, 9, 8, 18, 17, 38)
test_that("Do we get expected plot_libsize_prepost information?", {
  expect_equal(actual_class, expected_class)
  expect_equal(actual_numbers, expected_numbers)
})

## plot_significant_bar is run in extract_significant_genes.  Setting up to run it
## is a pita, I am skipping that here therefore.

## exported functions which will probably not get tested: plot_rpm

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 225plot_bar.R in ", elapsed,  " seconds.")
