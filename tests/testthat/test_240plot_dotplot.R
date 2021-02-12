start <- as.POSIXlt(Sys.time())
context("240plot_dotplot.R:\n")

pombe_expt <- make_pombe_expt(annotation = FALSE)
pombe_subset <- subset_expt(
  pombe_expt,
  subset = "minute == 0 | minute == 15 | minute == 30")
pombe_norm <- normalize_expt(pombe_expt, transform = "log2", convert = "cpm", filter = TRUE)

pombe_estimates <- all_adjusters(pombe_norm, estimate_type = "sva")
adjustments <- pombe_estimates[["model_adjust"]]

test <- plot_svfactor(pombe_norm, adjustments)
expected <- "gg"
actual <- class(test)[1]
test_that("Can we plot surrogate variables by factor?", {
  expect_equal(expected, actual)
})

test <- plot_batchsv(pombe_norm, adjustments)
actual1 <- class(test[[1]])[1]
actual2 <- class(test[[2]])[1]
actual3 <- class(test[[3]])[1]
test_that("Can we plot surrogates by batch?", {
  expect_equal(expected, actual1)
  expect_equal(expected, actual2)
  expect_equal(expected, actual3)
})

pombe_pca <- plot_pca(pombe_norm)
pombe_pcs <- pombe_pca[["table"]]
test <- plot_pcfactor(pombe_pcs, pombe_norm)
actual <- class(test)[1]
test_that("Can we plot principle components by factor?", {
  expect_equal(expected, actual)
})

test <- plot_sm(pombe_norm)
actual <- class(test[["plot"]])[1]
test_that("Can we plot standard medians?", {
  expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 240plot_dotplot.R in ", elapsed,  " seconds."))
