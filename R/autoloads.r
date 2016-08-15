#' Grab a copy of all bioconductor packages and install them by type
#'
#' This uses jsonlite to get a copy of all bioconductor packages by name and then iterates through
#' them with BiocInstaller to install all of them.  It performs a sleep between each installation in
#' an attempt to avoid being obnoxious.  As a result, it will of a necessity take forever.
#'
#' @param release  Bioconductor release to use, should probably be adjusted to automatically find
#'     it.
#' @param mirror  Bioconductor mirror to use.
#' @param base  Base directory on the mirror to download from.
#' @param type  Type in the tree to use (software, annotation, or experiment)
#' @param suppressUpdates  For BiocLite(), don't update?
#' @param suppressAutoUpdate  For BiocLite(), don't update?
#' @param force  Install if already installed?
#' @return a number of packages installed
#' @export
bioc_all <- function(release="3.4", mirror="bioc.ism.ac.jp", base="packages", type="software",
    suppress_updates=TRUE, suppress_auto=TRUE, force=FALSE) {
    dl_url <- paste0("https://", mirror, "/", base, "/json/", release, "/tree.json")
    ## message(paste0("DL: ", dl_url))
    ## https://bioc.ism.ac.jp/packages/json/3.3/tree.json
    pkg_list <- jsonlite::fromJSON(dl_url)
    pkg_names <- pkg_list[["attr"]][["packageList"]]
    software_names <- strsplit(x=pkg_names[[1]], split=",")[[1]]
    annotation_names <- strsplit(x=pkg_names[[2]], split=",")[[1]]
    experiment_names <- strsplit(x=pkg_names[[3]], split=",")[[1]]
    ret <- 0
    attempt <- function(pkg, update=suppress_updates, auto=suppress_auto, forceme=force) {
        message(paste0("Installing: ", pkg))
        if (isTRUE(forceme)) {
            installedp <- sm(try(BiocInstaller::biocLite(pkg, ask=FALSE,
                                                         suppressUpdates=update,
                                                         suppressAutoUpdate=auto)))
            if (class(installedp) != "try-error") {
                ret <- 1
            } else {
                ret <- 0
            }
        } else {
            if (isTRUE(pkg %in% .packages(all.available=TRUE))) {
                message(paste0("Package ", pkg,  " is already installed."))
                return(1)
            } else {
                installedp <- sm(try(BiocInstaller::biocLite(pkg, ask=FALSE,
                                                             suppressUpdates=update,
                                                             suppressAutoUpdate=auto)))
                if (class(installedp) != "try-error") {
                    ret <- 1
                } else {
                    ret <- 0
                }
            }
        }
        Sys.sleep(20)
        return(ret)
    } ## End attempt
    if (type == "software") {
        for (pkg in software_names) {
            ret <- ret + attempt(pkg)
        }
    } else if (type == "annotation") {
        for (pkg in annotation_names) {
            ret <- ret + attempt(pkg)
        }
    } else if (type == "experiment") {
        for (pkg in experiment_names) {
            ret <- ret + attempt(pkg)
        }
    } else {
        soft_num <- bioc_all(release=release, mirror=mirror, base=base, type="software")
        ret <- ret + soft_num
        anno_num <- bioc_all(release=release, mirror=mirror, base=base, type="annotation")
        ret <- ret + anno_num
        expt_num <- bioc_all(release=release, mirror=mirror, base=base, type="experiment")
        ret <- ret + expt_num
    }
    return(ret)
}

#' Automatic loading and/or installing of packages.
#'
#' Load a library, install it first if necessary.
#'
#' This was taken from:
#' http://sbamin.com/2012/11/05/tips-for-working-in-r-automatically-install-missing-package/
#'
#' @param lib String name of a library to check/install.
#' @param update Update packages?
#' @return 0 or 1, whether a package was installed or not.
#' @seealso \link[BiocInstaller]{biocLite} \link{install.packages}
#' @examples
#' \dontrun{
#' require.auto("ggplot2")
#' }
#' @export
require.auto <- function(lib, update=FALSE) {
    count <- 0
    local({r <- getOption("repos")
           r["CRAN"] <- "http://cran.r-project.org"
           options(repos=r)
       })
    if (isTRUE(update)) {
        utils::update.packages(ask=FALSE)
    }
    github_path <- NULL
    ## If there is a / in the library's name, assume it is a github path
    split_lib <- strsplit(x=lib, split="/")[[1]]
    if (length(split_lib) == 2) {
        github_path <- lib
        lib <- split_lib[[2]]
    }
    if (!isTRUE(lib %in% .packages(all.available=TRUE))) {
        ##eval(parse(text=paste("suppressPackageStartupMessages(require(", lib, "))", sep="")))
        if (is.null(github_path)) {
            source("http://bioconductor.org/biocLite.R")
            ##biocLite(character(), ask=FALSE) # update dependencies, if any.
            eval(parse(text=paste("biocLite('", lib, "')", sep="")))
            count <- 1
            ##eval(parse(text=paste("suppressPackageStartupMessages(require(", lib, "))", sep="")))
            ## eval(parse(text=paste("install.packages('", lib, "')", sep="")))
        } else {
            ret <- try(devtools::install_github(github_path))
            if (class(ret) != "try-error") {
                count <- 1
            }
        }
    }
    return(count)
}

autoloads_github <- function() {
    count <- 0
    count <- count + require.auto("seandavi/GEOquery")
    count <- count + require.auto("kasperdanielhansen/minfi")
}

autoloads_ontology <- function() {
    count <- 0
    count <- count + require.auto("clusterProfiler")
    count <- count + require.auto("GO.db")
    count <- count + require.auto("DOSE")
    count <- count + require.auto("gProfileR")
    count <- count + require.auto("goseq")
    count <- count + require.auto("GOstats")
    count <- count + require.auto("GSEABase")
    count <- count + require.auto("KEGGREST")
    count <- count + require.auto("pathview")
    count <- count + require.auto("RamiGO")
    count <- count + require.auto("topGO")
    return(count)
}

autoloads_genome <- function() {
    count <- 0
    count <- count + require.auto("biomaRt")
    count <- count + require.auto("BSgenome")
    count <- count + require.auto("genomeIntervals")
    count <- count + require.auto("rtracklayer")
    count <- count + require.auto("OrganismDbi")
    count <- count + require.auto("AnnotationHub")
    count <- count + require.auto("AnnotationDbi")
    return(count)
}

autoloads_elsayedlab <- function() {
    count <- 0
    count <- count + require.auto("TxDb.TcruziCLBrener.tritryp24.genes",
    "elsayed-lab/TxDb.TcruziCLBrener.tritryp24.genes")
    count <- count + require.auto("TxDb.TcruziCLBrenerEsmer.tritryp24.genes",
    "elsayed-lab/TxDb.TcruziCLBrenerEsmer.tritryp24.genes")
    count <- count + require.auto("TxDb.TcruziCLBrenerNonEsmer.tritryp9.genes",
    "elsayed-lab/TxDb.TcruziCLBrenerNonEsmer.tritryp9.genes")
    count <- count + require.auto("TxDb.LmajorFriedlin.tritryp9.genes",
    "elsayed-lab/TxDb.LmajorFriedlin.tritryp9.genes")
    count <- count + require.auto("BSgenome.Lmajor.friedlin",
    "elsayed-lab/BSgenome.Lmajor.friedlin")
    count <- count + require.auto("BSgenome.Tcruzi.clbrener",
    "elsayed-lab/BSgenome.Tcruzi.clbrener")
    count <- count + require.auto("BSgenome.Tcruzi.esmeraldo",
    "elsayed-lab/BSgenome.Tcruzi.esmeraldo")
    count <- count + require.auto("BSgenome.Tcruzi.nonesmeraldo",
    "elsayed-lab/BSgenome.Tcruzi.nonesmeraldo")
    count <- count + require.auto("org.TcCLB.nonesmer.tritryp.db",
    "elsayed-lab/org.TcCLB.nonesmer.tritryp.db")
    count <- count + require.auto("org.TcCLB.esmer.tritryp.db",
    "elsayed-lab/org.TcCLB.esmer.tritryp.db")
    count <- count + require.auto("org.TcCLB.tritryp.db", "elsayed-lab/org.TcCLB.tritryp.db")
    count <- count + require.auto("org.LmjF.tritryp.db", "elsayed-lab/org.LmjF.tritryp.db")
    count <- count + require.auto("Trypanosoma.cruzi.CLBrener",
    "elsayed-lab/Trypanosoma.cruzi.CLBrener")
    count <- count + require.auto("Trypanosoma.cruzi.CLBrener.Esmeraldo",
    "elsayed-lab/Trypanosoma.cruzi.CLBrener.Esmeraldo")
    count <- count + require.auto("Leishmania.major.Friedlin",
                                  "elsayed-lab/Leishmania.major.Friedlin")
    return(count)
}

autoloads_deseq <- function() {
    count <- 0
    count <- count + require.auto("affy")
    count <- count + require.auto("preprocessCore")
    count <- count + require.auto("DESeq2")
    count <- count + require.auto("DESeq")
    count <- count + require.auto("edgeR")
    count <- count + require.auto("limma")
    count <- count + require.auto("RUVSeq")
    count <- count + require.auto("sva")
    count <- count + require.auto("survJamda")
    count <- count + require.auto("pasilla")
    count <- count + require.auto("preprocessCore")
    return(count)
}

autoloads_graphs <- function() {
    count <- 0
    count <- count + require.auto("Cairo")
    count <- count + require.auto("corrplot")
    count <- count + require.auto("directlabels")
    count <- count + require.auto("ggplot2")
    count <- count + require.auto("googleVis")
    count <- count + require.auto("gplots")
    count <- count + require.auto("ggrepel")
    count <- count + require.auto("grid")
    count <- count + require.auto("gridExtra")
    count <- count + require.auto("RColorBrewer")
    count <- count + require.auto("Rgraphviz")
    return(count)
}

autoloads_helpers <- function() {
    count <- 0
    count <- count + require.auto("MASS")
    count <- count + require.auto("mgcv")
    count <- count + require.auto("Matrix")
    count <- count + require.auto("matrixStats")
    count <- count + require.auto("devtools")
    count <- count + require.auto("dplyr")
    count <- count + require.auto("BiocParallel")
    count <- count + require.auto("data.table")
    count <- count + require.auto("ffpe")
    count <- count + require.auto("gtools")
    count <- count + require.auto("hash")
    count <- count + require.auto("knitcitations")
    count <- count + require.auto("knitr")
    count <- count + require.auto("methods")
    count <- count + require.auto("plyr")
    count <- count + require.auto("RCurl")
    count <- count + require.auto("reshape")
    count <- count + require.auto("rjson")
    count <- count + require.auto("rmarkdown")
    count <- count + require.auto("roxygen2")
    count <- count + require.auto("testthat")
    count <- count + require.auto("tools")
    count <- count + require.auto("openxlsx")
    count <- count + require.auto("xtable")
    return(count)
}

autoloads_stats <- function() {
    count <- 0
    count <- count + require.auto("multtest")
    count <- count + require.auto("qvalue")
    count <- count + require.auto("robust")
    return(count)
}

autoloads_misc <- function() {
    count <- 0
    count <- count + require.auto("corpcor")
    count <- count + require.auto("Rsamtools")
    count <- count + require.auto("scales")
    count <- count + require.auto("seqinr")
    count <- count + require.auto("ReactomePA")
    count <- count + require.auto("mygene")
    return(count)
}

autoloads_motif <- function() {
    count <- 0
    count <- count + require.auto("rGADEM")
###    require.auto("MotIV")
###    require.auto("motifRG")
###    require.auto("motifStack")
    return(count)
}

#' Automatic installation of stuff I use.
#'
#' hpgltools uses too many packages, ~65 at last count, I use autoloads_all() to make sure they are
#' all installed.
#'
#' @param update Update installed packages?
#' @return The number of installed packages
#' @seealso \link[BiocInstaller]{biocLite} \link{install.packages}
#' @export
autoloads_all <- function(update=FALSE) {
    helpers <- autoloads_helpers()
    misc <- autoloads_misc()
    genome <- autoloads_genome()
    graphs <- autoloads_graphs()
    stats <- autoloads_stats()
    deseq <- autoloads_deseq()
    ontology <- autoloads_ontology()
    motif <- autoloads_motif()
    if (isTRUE(update)) {
        utils::update.packages(ask=FALSE)
    }
    packages_installed <- helpers + misc + genome + graphs + stats + deseq + ontology + motif
    message(paste0("autoloads_all() installed ", packages_installed, " packages."))
    return(packages_installed)
}

## EOF
