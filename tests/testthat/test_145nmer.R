start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("145nmer.R:\n")
## 2017-12, exported functions in annotation_gff:

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 145nmer.R in ", elapsed,  " seconds."))
