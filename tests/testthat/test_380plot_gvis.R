start <- as.POSIXlt(Sys.time())
context("380plot_gvis.R:\n")

deseq <- new.env()
load("324_de_deseq.rda", envir=deseq)
stop("I need this to stop for now.")

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 380plot_gvis.R in ", elapsed,  " seconds."))
