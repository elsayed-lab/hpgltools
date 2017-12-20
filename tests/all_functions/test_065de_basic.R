start <- as.POSIXlt(Sys.time())
context("065de_basic.R:\n")
## 2017-12, exported functions in de_basic:

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 065de_basic.R in ", elapsed,  " seconds."))
