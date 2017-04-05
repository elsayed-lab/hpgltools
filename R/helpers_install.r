#' Grab a copy of all bioconductor packages and install them by type
#'
#' This uses jsonlite to get a copy of all bioconductor packages by name and then iterates through
#' them with BiocInstaller to install all of them.  It performs a sleep between each installation in
#' an attempt to avoid being obnoxious.  As a result, it will of a necessity take forever.
#'
#' @param release  Bioconductor release to use, should probably be adjusted to automatically find
#'  it.
#' @param mirror  Bioconductor mirror to use.
#' @param base  Base directory on the mirror to download from.
#' @param type  Type in the tree to use (software or annotation)
#' @param suppress_updates  For BiocLite(), don't update?
#' @param suppress_auto  For BiocLite(), don't update?
#' @param force  Install if already installed?
#' @return a number of packages installed
#' @seealso \pkg{BiocInstaller}
#' @examples
#' \dontrun{
#'  go_get_some_coffee_this_will_take_a_while <- bioc_all()
#' }
#' @export
bioc_all <- function(release="3.4", mirror="bioconductor.statistik.tu-dortmund.de", base="packages", type="software",
    suppress_updates=TRUE, suppress_auto=TRUE, force=FALSE) {
    dl_url <- paste0("https://", mirror, "/", base, "/json/", release, "/tree.json")
    ## message(paste0("DL: ", dl_url))
    ## dl_url <- "https://bioc.ism.ac.jp/packages/json/3.3/tree.json"
    suc <- c()
    fail <- c()
    alr <- c()
    pkg_list <- jsonlite::fromJSON(dl_url)
    pkg_names <- pkg_list[["attr"]][["packageList"]]
    software_names <- strsplit(x=pkg_names[[1]], split=",")[[1]]
    annotation_names <- strsplit(x=pkg_names[[2]], split=",")[[1]]
    ## experiment_names <- strsplit(x=pkg_names[[3]], split=",")[[1]]
    ## It appears that with bioconductor release 3.4, experiment has been folded into annotation.
    installed <- list(succeeded=c(), failed=c(), already=c())
    attempt <- function(pkg, update=suppress_updates, auto=suppress_auto, forceme=force,
                        state=list(succeeded=c(), failed=c(), already=c())) {
        sleep <- 10
        suc <- state[["succeeded"]]
        fail <- state[["failed"]]
        alr <- state[["already"]]
        message(paste0("Installing: ", pkg))
        if (isTRUE(forceme)) {
            installedp <- sm(try(BiocInstaller::biocLite(pkg, ask=FALSE,
                                                         suppressUpdates=update,
                                                         suppressAutoUpdate=auto)))
            if (class(installedp) == "try-error") {
                fail <- append(fail, pkg)
            } else {
                suc <- append(suc, pkg)
            }
        } else {
            if (isTRUE(pkg %in% .packages(all.available=TRUE))) {
                message(paste0("Package ", pkg,  " is already installed."))
                alr <- append(alr, pkg)
                sleep <- 0
            } else {
                installedp <- try(sm(BiocInstaller::biocLite(pkg, ask=FALSE,
                                                             suppressUpdates=update,
                                                             suppressAutoUpdate=auto)))
                if (class(installedp) == "try-error") {
                    fail <- append(fail, pkg)
                } else {
                    suc <- append(suc, pkg)
                }
            }
        }
        Sys.sleep(sleep)
        ret <- list(
            "succeeded" = suc,
            "failed" = fail,
            "already" = alr)
        return(ret)
    } ## End attempt
    if (type == "software") {
        for (pkg in software_names) {
            installed <- attempt(pkg, state=installed)
        }
    } else if (type == "annotation") {
        for (pkg in annotation_names) {
            installed <- attempt(pkg, state=installed)
        }
    } else {
        software_installed <- bioc_all(release=release, mirror=mirror, base=base, type="software")
        annotation_installed <- bioc_all(release=release, mirror=mirror, base=base, type="annotation")
        installed <- list(
            "software" = software_installed,
            "annotation" = annotation_installed)
    }
    return(installed)
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
#' @seealso \pkg{BiocInstaller}
#'  \code{\link[BiocInstaller]{biocLite}} \code{\link{install.packages}}
#' @examples
#' \dontrun{
#'  require.auto("ggplot2")
#' }
#' @export
require.auto <- function(lib, update=FALSE) {
    count <- 0
    local({
        r <- getOption("repos")
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

## EOF
