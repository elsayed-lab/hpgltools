#' Automatic loading and/or installing of packages.
#'
#' \code{require.auto} loads a library, and installs it first if necessary.
#'
#' This was taken from:
#' http://sbamin.com/2012/11/05/tips-for-working-in-r-automatically-install-missing-package/
#' 
#' @param lib string name of a library
#' @return NULL currently
#' @seealso \code{\link{biocLite}} and \code{\link{install.packages}}
#' @export
#' @examples
#' ## require.auto("ggplot2")
require.auto = function(lib, github_path=NULL, verbose=FALSE, update=FALSE) {
    local({r <- getOption("repos")
           r["CRAN"] <- "http://cran.r-project.org" 
           options(repos=r)
       })
    if (isTRUE(update)) {
        update.packages(ask=FALSE)
    }
    if (isTRUE(lib %in% .packages(all.available=TRUE))) {
        if (verbose) {
            message(sprintf("Loading %s", lib))
        }
        eval(parse(text=paste("suppressPackageStartupMessages(require(", lib, "))", sep="")))
    } else {
        if (is.null(github_path)) {
            source("http://bioconductor.org/biocLite.R")
            biocLite(character(), ask=FALSE) # update dependencies, if any.
            eval(parse(text=paste("biocLite('", lib, "')", sep="")))
            if (verbose) {
                message(sprintf("Loading %s", lib))
            }
            eval(parse(text=paste("suppressPackageStartupMessages(require(", lib, "))", sep="")))
            ## eval(parse(text=paste("install.packages('", lib, "')", sep="")))
        } else {
            library(devtools)
            install_github(github_path)
        }
    }
}

autoloads_ontology = function(...) {
    require.auto("clusterProfiler", ...)
    require.auto("GO.db")
    require.auto("DOSE")
    require.auto("goseq")
    require.auto("GOstats")
    require.auto("GSEABase")
    require.auto("KEGGREST")    
    require.auto("pathview")
    require.auto("RamiGO")
    require.auto("topGO")    
}

autoloads_genome = function(...) {
    require.auto("biomaRt", ...)    
    require.auto("BSgenome")
    require.auto("BSgenome.Lmajor.friedlin", "elsayed-lab/BSgenome.Lmajor.friedlin")
    require.auto("genomeIntervals")    
    require.auto("rtracklayer")
}

autoloads_deseq = function(...) {
    require.auto("preprocessCore")
    require.auto("cbcbSEQ", "kokrah/cbcbSEQ")
    require.auto("DESeq2")
    require.auto("DESeq")
    require.auto("edgeR")
    require.auto("sva")    
}

autoloads_graphs = function(...) {
    require.auto("Cairo")
    require.auto("directlabels")
    require.auto("ggplot2")
    require.auto("googleVis")
    require.auto("gplots")
    require.auto("grid")
    require.auto("gridExtra")
    require.auto("RColorBrewer")
    require.auto("Rgraphviz")    
}

autoloads_helpers = function(...) {
    require.auto("MASS", ...)
    require.auto("mgcv")
    require.auto("Matrix")
    require.auto("devtools", ...)
    require.auto("BiocParallel")
    register(MulticoreParam(4))    
    require.auto("data.table")
    require.auto("gtools")
    require.auto("hash")    
    require.auto("knitcitations")
    require.auto("knitr")
    require.auto("knitrBootstrap", "jimhester/knitrBootstrap")
    require.auto("methods")
    require.auto("plyr")
    require.auto("RCurl")
    require.auto("reshape")
    require.auto("rjson")
    require.auto("rmarkdown")
    require.auto("roxygen2")
    require.auto("testthat")
    options(java.parameters = "-Xmx4g")  ## used for xlconnect
    require.auto("XLConnect")
    require.auto("xtable")    
}

autoloads_stats = function(...) {
    require.auto("multtest", ...)
    require.auto("qvalue")
    require.auto("robust")    
}

autoloads_misc = function(...) {
    require.auto("motifRG", ...)
    require.auto("Rsamtools")
    require.auto("scales")
    require.auto("seqinr")    
}

#' Automatic loading of stuff I use
#'
#' @return NULL currently
#' @seealso \code{\link{biocLite}} and \code{\link{install.packages}}
#' @export
autoloads_all = function(...) {
    ## I added the ... to the first entry
    ## So that it will do update.packages() if update=TRUE
    ## I don't need it after the first I think.
    autoloads_helpers(...)
    autoloads_ontology(...)
    autoloads_genome(...)
    autoloads_deseq(...)
    autoloads_graphs(...)
    autoloads_stats(...)
    autoloads_misc(...)
    ##cite_options(tooltip=TRUE)
    ##cleanbib()
    options(gvis.plot.tag="chart")
    mainfont = "Helvetica"
    ##Cairo()
    CairoFonts(regular = paste(mainfont, "style=Regular", sep = ":"),
               bold = paste(mainfont, "style=Bold", sep = ":"),
               italic = paste("SimSun", "style=Regular", sep = ":"), 
               bolditalic = paste(mainfont, "style=Bold Italic,BoldItalic", sep = ":"))
    pdf = CairoPDF
    png = CairoPNG
    x11 = CairoX11
    svg = CairoSVG
}
