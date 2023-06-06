start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("145nmer.R")
## 2017-12, exported functions in annotation_gff:

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 145nmer.R in ", elapsed,  " seconds.")
