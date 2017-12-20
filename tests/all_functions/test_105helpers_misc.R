start <- as.POSIXlt(Sys.time())
context("105helpers_misc.R:\n")
## 2017-12, exported functions in annotation_gff:

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 105helpers_misc.R in ", elapsed,  " seconds."))
