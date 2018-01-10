start <- as.POSIXlt(Sys.time())
context("190ontology_gostats.R:\n")
## 2017-12, exported functions in annotation_gff:

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 190ontology_gostats.R in ", elapsed,  " seconds."))
