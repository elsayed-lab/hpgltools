start <- as.POSIXlt(Sys.time())
context("115model_pca.R:\n")
## 2017-12, exported functions in model_pca:
## pca_information() pca_highscores() pcRes() plot_pca()
## plot_pcs() test_pca_methods()

pombe_expt <- make_pombe_expt(annotation=FALSE)

testing <- sm(pca_information(pombe_expt, plot_pcas=TRUE,
                              expt_factors=c("strain", "minute", "replicate")))
expected <- c(94.150, 4.612, 0.791, 0.120, 0.108, 0.081)
actual <- head(testing[["rsquared_table"]][["prop_var"]])
test_that("pca_information() provides a rsquared table?", {
  expect_equal(expected, actual)
})
               


end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 115model_pca.R in ", elapsed,  " seconds."))
