start <- as.POSIXlt(Sys.time())
context("115model_pca.R:\n")
## 2017-12, exported functions in model_pca:
## pca_information() pca_highscores() pcRes() plot_pca()
## plot_pcs() test_pca_methods()

pombe_expt <- make_pombe_expt(annotation=FALSE)

## 01 pca_information()
testing <- sm(pca_information(pombe_expt, plot_pcas=TRUE,
                              expt_factors=c("strain", "minute", "replicate")))
expected <- c(94.150, 4.612, 0.791, 0.120, 0.108, 0.081)
actual <- head(testing[["rsquared_table"]][["prop_var"]])
test_that("pca_information() provides a rsquared table?", {
  expect_equal(expected, actual)
})

## 02 pca_highscores()
testing <- pca_highscores(pombe_expt)
expected <- c(4.0946442, 2.3214251, 1.0686429, 1.1414550, 2.7090227, 0.5437997)
actual <- head(as.numeric(testing[["scores"]][, "Comp.1"]))
test_that("pca_highscores() provides some scores by PC?", {
  expect_equal(expected, actual)
})

## 03 pcRes()  This is called from plot_pca() and friends, test it there.
## 04 plot_pca()
## 05 plot_pcs()  This is also called from plot_pca()
testing <- plot_pca(pombe_expt)
expected <- c(6008536.5, 1329859.1, 550763.0, 214403.7, 203382.4, 176304.0)
actual <- head(testing$pca$d)
test_that("plot_pca() provides expected SVD data?", {
  expect_equal(expected, actual, tolerance=1.0)
})

test_that("plot_pca() provides a plot!?", {
  expect_equal(class(testing$plot)[1], "gg")
})

expected <- c(-0.007178318, 0.070735209, -0.136645507, -0.374502211, 0.072181094, -0.028458656)
actual <- head(testing$table)$PC1
test_that("We acquire expected values for PC1?", {
  expect_equal(expected, actual)
})

expected <- c(0.20199545, -0.19923795, -0.10775070, 0.28537150, -0.06148581, -0.25190759)
actual <- head(testing$table)$PC2
test_that("We acquire expected values for PC2?", {
  expect_equal(expected, actual)
})

expected <- c(94.15, 4.61, 0.79, 0.12, 0.11, 0.08)
actual <- head(testing$variance)
test_that("We get variances from pcRes?", {
  expect_equal(expected, actual)
})

## 06 test_pca_methods()
testing <- sm(test_pca_methods(pombe_expt))
expected <- c("svd", "ppca", "bpca", "svdi", "nipals", "tables", "plots")
actual <- names(testing)
test_that("We get some information from test_pca_methods?", {
  expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 115model_pca.R in ", elapsed,  " seconds."))
