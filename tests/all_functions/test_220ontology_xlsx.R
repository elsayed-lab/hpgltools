start <- as.POSIXlt(Sys.time())
context("220ontology_xlsx.R:\n")
## 2017-12, exported functions in annotation_gff:

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 220ontology_xlsx.R in ", elapsed,  " seconds."))
