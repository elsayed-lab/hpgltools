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
bioc_all <- function(release="3.5", mirror="bioconductor.statistik.tu-dortmund.de", base="packages", type="software",
                     suppress_updates=TRUE, suppress_auto=TRUE, force=FALSE) {
  dl_url <- paste0("https://", mirror, "/", base, "/json/", release, "/tree.json")
  ## message(paste0("DL: ", dl_url))
  ## dl_url <- "https://bioc.ism.ac.jp/packages/json/3.3/tree.json"
  suc <- c()
  ## Sadly, biocLite() does not give different returns for a successfully/failed install.
  ## Instead it always returns a character of the package(s) asked to install.
  ## That seems dumb to me, but what do I know.
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
        install_directory <- paste0(.libPaths()[1], "/", pkg)
        if (file.exists(install_directory)) {
          suc <- append(suc, pkg)
        } else {
          fail <- append(fail, pkg)
        }
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
          install_directory <- paste0(.libPaths()[1], "/", pkg)
          if (file.exists(install_directory)) {
            suc <- append(suc, pkg)
          } else {
            fail <- append(fail, pkg)
          }
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

#'  Install the set of local packrat packages so everyone may use them!
#' 
#' @export
install_packrat_globally <- function() {
  packrat_installed <- packrat::status()
  packrat_location <- packrat:::getProjectDir()
  packrat_src <- packrat:::srcDir(packrat_location)
  paths <- sm(packrat::packrat_mode())
  num_installed <- nrow(packrat_installed)
  globally_installed <- installed.packages()

  newly_installed <- 0
  message(paste0("Going to check on/install ", num_installed, " packages."))
  for (c in 1:num_installed) {
    pkg_name <- packrat_installed[c, "package"]
    pkg_ver <- packrat_installed[c, "packrat.version"]
    if (is.na(pkg_ver)) {
      next
    }
    packrat_package_path <- paste0(packrat_src, "/",
                                   pkg_name, "/",
                                   pkg_name, "_", pkg_ver, ".tar.gz")
    if (pkg_name %in% rownames(globally_installed)) {
      global_version <- globally_installed[pkg_name, "Version"]
      if (global_version == pkg_ver) {
        message(paste0("Package: ", pkg_name, " is globally installed as the same version."))
      } else {
        message(paste0("Package: ", pkg_name, " is globally installed as version: ",
                       global_version, "; packrat has version ", pkg_ver, "."))
        inst <- try(devtools::install_url(paste0("file://", packrat_package_path)))
        if (class(inst) != "try-error") {
          newly_installed <- newly_installed + 1
        }
      }
    } else {
      message(paste0("Package: ", pkg_name, " is not installed."))
      inst <- try(devtools::install_url(paste0("file://", packrat_package_path)))
      if (class(inst) != "try-error") {
        newly_installed <- newly_installed + 1
      }
    }
  } ## End of the for loop
  paths <- sm(packrat::packrat_mode())
  message(paste0("Installed ", newly_installed, " packages."))
  return(newly_installed)
}

## EOF
