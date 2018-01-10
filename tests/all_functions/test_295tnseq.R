start <- as.POSIXlt(Sys.time())
context("295tnseq.R:\n")
## 2017-12, exported functions in annotation_gff:

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 295tnseq.R in ", elapsed,  " seconds."))
