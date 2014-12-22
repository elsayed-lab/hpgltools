## Time-stamp: "Mon Dec  1 11:48:37 2014 Ashton Trey Belew (abelew@gmail.com)"
## autoloads.R contains some short-cuts I wrote for myself to make
## installing/maintaining packages/dependencies easier
## 
## The following are a few reminders that I keep forgetting:
## clean_environment = environment()
## devtools::document() ## Scan for roxygen2 comments and generate .Rd files
## Reminder about roxygen comments:
##  The first sentence is the title.  The second paragraph is the
##  description.  The third paragraph are details.
##  Paragraphs are separated by single #' lines
## Use load_all() to reload a library which is being currently edited.

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
require.auto = function(lib, github_path=NULL, verbose=TRUE, update=FALSE) {
    if (isTRUE(update)) {
        update.packages(ask=FALSE)
    }
    if (isTRUE(lib %in% .packages(all.available=TRUE))) {
        if (verbose) {
            print(sprintf("Loading %s", lib))
        }
        eval(parse(text=paste("suppressPackageStartupMessages(require(", lib, "))", sep="")))
    } else {
        if (is.null(github_path)) {
            source("http://bioconductor.org/biocLite.R")
            biocLite(character(), ask=FALSE) # update dependencies, if any.
            eval(parse(text=paste("biocLite('", lib, "')", sep="")))
            if (verbose) {
                print(sprintf("Loading %s", lib))
            }
            eval(parse(text=paste("suppressPackageStartupMessages(require(", lib, "))", sep="")))
            ## eval(parse(text=paste("install.packages('", lib, "')", sep="")))
        } else {
            library(devtools)
            install_github(github_path)
        }
    }
}

#' Automatic loading of stuff I use
#'
#' @return NULL currently
#' @seealso \code{\link{biocLite}} and \code{\link{install.packages}}
#' @export
autoloads = function(...) {
    ## I added the ... to the first entry
    ## So that it will do update.packages() if update=TRUE
    ## I don't need it after the first I think.
    require.auto("biomaRt", ...)
    require.auto("BSgenome")
    require.auto("BSgenome.Lmajor.friedlin")
    require.auto("Cairo")
    require.auto("cbcbSEQ")
    require.auto("clusterProfiler")
    require.auto("data.table")
    require.auto("DESeq2")
    require.auto("DESeq")
    require.auto("devtools")
    require.auto("directlabels")
    require.auto("DOSE")
    require.auto("edgeR")
    require.auto("genomeIntervals")
    require.auto("ggplot2")
    require.auto("GO.db")
    require.auto("googleVis")
    require.auto("goseq")
    require.auto("gplots")
    require.auto("gtools")
    require.auto("gridExtra")
    require.auto("hash")
    require.auto("KEGGREST")
    require.auto("knitcitations")
    require.auto("knitr")
    require.auto("knitrBootstrap", "jimhester/knitrBootstrap")
    require.auto("methods")
    require.auto("motifRG")
    require.auto("multtest")
    require.auto("pathview")
    require.auto("plyr")
    require.auto("qvalue")
    require.auto("RamiGO")
    require.auto("RColorBrewer")
    require.auto("reshape")
    require.auto("Rgraphviz")
    require.auto("rjson")
    require.auto("rmarkdown")
    require.auto("robust")
    require.auto("roxygen2")
    require.auto("Rsamtools")
    require.auto("rtracklayer")
    require.auto("scales")
    require.auto("seqinr")
    require.auto("sva")
    require.auto("testthat")             
    require.auto("topGO")
    options(java.parameters = "-Xmx4g")  ## used for xlconnect
    require.auto("XLConnect")
    require.auto("xtable")
    ##cite_options(tooltip=TRUE)
    ##cleanbib()
    options(gvis.plot.tag="chart")
    mainfont = "Helvetica"
    ##Cairo()
    CairoFonts(regular = paste(mainfont, "style=Regular", sep = ":"),
               bold = paste(mainfont, "style=Bold", sep = ":"),
               italic = paste("SimSun", "style=Regular", sep = ":"), 
               bolditalic = paste(mainfont, "style=Bold Italic,BoldItalic", sep = ":"))
    workbook = loadWorkbook("excel/workbook.xls", create = TRUE)
    opts_chunk$set(fig.width=800/192, fig.height=800/192, dpi=192, dev="png", bootstrap.thumbnail.size=12)

    theme_set(theme_bw())
    pdf = CairoPDF
    png = CairoPNG
    x11 = CairoX11
    svg = CairoSVG
}
