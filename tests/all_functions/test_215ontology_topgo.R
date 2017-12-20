start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("215ontology_topgo.R:\n")
## 2017-12, exported functions in annotation_gff:

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 215ontology_topgo.R in ", elapsed,  " seconds."))
