start <- as.POSIXlt(Sys.time())
context("095de_xlsx.R: Test functions in annotation_genbank.r")

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 095de_xlsx.R in ", elapsed,  " seconds."))
