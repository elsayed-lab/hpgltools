## helpers_misc.r: Some random helper functions.

#' Make a backup of an existing file with n revisions, like VMS!
#'
#' Sometimes I just want to kick myself for overwriting important files and then
#' I remember using VMS and wish modern computers were a little more like it.
#'
#' @param backup_file Filename to backup.
#' @param backups How many revisions?
backup_file <- function(backup_file, backups = 4) {
  if (file.exists(backup_file)) {
    for (i in seq(from = backups, to = 1)) {
      j <- i + 1
      i <- sprintf("%02d", i)
      j <- sprintf("%02d", j)
      test <- glue("{backup_file}.{i}")
      new <- glue("{backup_file}.{j}")
      if (file.exists(test)) {
        file.rename(test, new)
      }
    }
    newfile <- glue("{backup_file}.{i}")
    message("Renaming ", backup_file, " to ", newfile, ".")
    file.copy(backup_file, newfile)
  } else {
    message("The file does not yet exist.")
  }
}

#' Grab a copy of all bioconductor packages and install them by type
#'
#' This uses jsonlite to get a copy of all bioconductor packages by name and
#' then iterates through them with BiocManager to install all of them.  It
#' performs a sleep between each installation in an attempt to avoid being
#' obnoxious.  As a result, it will of a necessity take forever.
#'
#' @param release Bioconductor release to use, should probably be adjusted to
#'  automatically find it.
#' @param mirror Bioconductor mirror to use.
#' @param base Base directory on the mirror to download from.
#' @param type Type in the tree to use (software or annotation)
#' @param suppress_updates For BiocLite(), don't update?
#' @param suppress_auto For BiocLite(), don't update?
#' @param force Install if already installed?
#' @return a number of packages installed
#' @seealso [BiocManager] [jsonlite]
#' @examples
#' \dontrun{
#'  go_get_some_coffee_this_will_take_a_while <- bioc_all()
#' }
#' @export
bioc_all <- function(release = NULL,
                     mirror = "bioconductor.statistik.tu-dortmund.de",
                     base = "packages", type = "software",
                     suppress_updates = TRUE, suppress_auto = TRUE, force = FALSE) {
  if (is.null(release)) {
    release <- as.character(BiocManager::version())
  }
  dl_url <- glue("https://{mirror}/{base}/json/{release}/tree.json")
  ## dl_url <- "https://bioc.ism.ac.jp/packages/json/3.3/tree.json"
  suc <- c()
  ## Sadly, biocLite() does not give different returns for a successfully/failed
  ## install. Instead it always returns a character of the package(s) requested.
  ## That seems dumb to me, but what do I know.
  fail <- c()
  alr <- c()
  ## Soemthing is weird with the libcurl on the cluster...
  test <- download.file(url = dl_url, destfile = "test.json")
  pkg_list <- jsonlite::fromJSON("test.json")
  pkg_names <- pkg_list[["attr"]][["packageList"]]
  software_names <- strsplit(x = pkg_names[[1]], split = ",")[[1]]
  annotation_names <- strsplit(x = pkg_names[[2]], split = ",")[[1]]
  ## experiment_names <- strsplit(x = pkg_names[[3]], split = ",")[[1]]
  ## It appears that with bioconductor release 3.4,
  ## experiment has been folded into annotation.
  installed <- list(succeeded = c(), failed = c(), already = c())
  attempt <- function(pkg, update = suppress_updates, auto = suppress_auto,
                      forceme = force,
                      state = list(succeeded = c(),
                                   failed = c(),
                                   already = c())) {
    sleep <- 10
    suc <- state[["succeeded"]]
    fail <- state[["failed"]]
    alr <- state[["already"]]
    message("Installing: ", pkg)
    if (isTRUE(forceme)) {
      installedp <- sm(try(BiocManager::install(pkg, ask = FALSE,
                                                suppressUpdates = update,
                                                suppressAutoUpdate = auto)))
      if (class(installedp) == "try-error") {
        fail <- append(fail, pkg)
      } else {
        install_directory <- glue("{.libPaths()[1]}/{pkg}")
        if (file.exists(install_directory)) {
          suc <- append(suc, pkg)
        } else {
          fail <- append(fail, pkg)
        }
      }
    } else {
      if (isTRUE(pkg %in% .packages(all.available = TRUE))) {
        message("Package ", pkg,  " is already installed.")
        alr <- append(alr, pkg)
        sleep <- 0
      } else {
        installedp <- try(sm(BiocManager::install(pkg, ask = FALSE,
                                                  suppressUpdates = update,
                                                  suppressAutoUpdate = auto)))
        if (class(installedp) == "try-error") {
          fail <- append(fail, pkg)
        } else {
          install_directory <- glue("{.libPaths()[1]}/{pkg}")
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
      installed <- attempt(pkg, state = installed)
    }
  } else if (type == "annotation") {
    for (pkg in annotation_names) {
      installed <- attempt(pkg, state = installed)
    }
  } else {
    software_installed <- bioc_all(release = release,
                                   mirror = mirror, base = base,
                                   type = "software")
    annotation_installed <- bioc_all(release = release,
                                     mirror = mirror, base = base,
                                     type = "annotation")
    installed <- list(
      "software" = software_installed,
      "annotation" = annotation_installed)
  }
  return(installed)
}

#' Clear an R session, this is probably unwise given what I have read about R.
#'
#' @param keepers List of namespaces to leave alone (unimplemented).
#' @param depth Cheesy forloop of attempts to remove packages stops after this
#'  many tries.
#' @return A spring-fresh R session, hopefully.
#' @seealso [R.utils]
#' @export
clear_session <- function(keepers = NULL, depth = 10) {
  ## Partially taken from:
  ## https://stackoverflow.com/questions/7505547/detach-all-packages-while-working-in-r
  basic_packages <- c("package:stats", "package:graphics", "package:grDevices",
                      "package:utils", "package:datasets", "package:methods",
                      "package:base")
  package_list <- search()[ifelse(unlist(gregexpr("package:", search()))==1, TRUE, FALSE)]
  package_list <- setdiff(package_list, basic_packages)
  result <- R.utils::gcDLLs()
  if (length(package_list) > 0) {
    for (package in package_list) {
      detach(package, character.only = TRUE)
    }
  }
  tt <- sm(rm(list = ls(all.names = TRUE), envir = globalenv()))
}

#' Similarity measure which combines elements from Pearson correlation and
#' Euclidean distance.
#'
#' Here is Keith's summary:
#' Where the cor returns the Pearson correlation matrix for the input matrix,
#' and the dist function returns the Euclidean distance matrix for the input
#' matrix. The LHS of the equation is simply the sign of the correlation
#' function, which serves to preserve the sign of the interaction. The RHS
#' combines the Pearson correlation and the log inverse Euclidean distance with
#' equal weights. The result is a number in the range from -1 to 1 where values
#' close to -1 indicate a strong negative correlation and values close to 1
#' indicate a strong positive corelation.  While the Pearson correlation and
#' Euclidean distance each contribute equally in the above equation, one could
#' also assign tuning parameters to each of the metrics to allow for unequal
#' contributions.
#'
#' @param data Matrix of data
#' @param cor_method Which correlation method to use?
#' @param dist_method Which distance method to use?
#' @param cor_weight 0-1 weight of the correlation, the distance weight
#'   will be 1-cor_weight.
#' @param ... extra arguments for cor/dist
#' @author Keigth Hughitt
#' @return Matrix of the correlation-modified distances of the original matrix.
#' @export
cordist <- function(data, cor_method = "pearson", dist_method = "euclidean",
                    cor_weight = 0.5, ...) {
  cor_matrix  <- hpgl_cor(t(data), method = cor_method, ...)
  dist_matrix <- as.matrix(dist(data, method = dist_method, diag = TRUE,
                                upper = TRUE))
  dist_matrix <- log1p(dist_matrix)
  dist_matrix <- 1 - (dist_matrix / max(dist_matrix))

  if (cor_weight > 1 | cor_weight < 0) {
    stop("Invalid correlation weight.")
  }
  dist_weight <- 1 - cor_weight
  result <- (cor_weight * abs(cor_matrix)) + (dist_weight * dist_matrix)
  result <- sign(cor_matrix) * result
  return(result)
}

#' Get the current git commit for hpgltools
#'
#' One might reasonably ask about this function: "Why?"  I invoke this function
#' at the end of my various knitr documents so that if necessary I can do a
#' > git reset <commit id> and get back to the exact state of my code.
#'
#' @param gitdir Directory containing the git repository.
#' @export
get_git_commit <- function(gitdir = "~/hpgltools") {
  cmdline <- glue("cd {gitdir} && git log -1 2>&1 | grep 'commit' | awk '{{print $2}}'")
  commit_result <- system(cmdline, intern = TRUE)
  cmdline <- glue("cd {gitdir} && git log -1 2>&1 | grep 'Date' | perl -pe 's/^Date:\\s+//g'")
  date_result <- system(cmdline, intern = TRUE)
  result <- glue("{date_result}: {commit_result}")
  message("If you wish to reproduce this exact build of hpgltools, invoke the following:")
  message("> git clone http://github.com/abelew/hpgltools.git")
  message("> git reset ", commit_result)
  return(result)
}

#' Create a contrast matrix suitable for MSstats and similar tools.
#'
#' I rather like makeContrasts() from limma.  I troubled me to have to manually
#' create a contrast matrix when using MSstats.  It turns out it troubled me for
#' good reason because I managed to reverse the terms and end up with the
#' opposite contrasts of what I intended.  Ergo this function.
#'
#' Feed make_simplified_contrast_matrix() a series of numerators and denominators names
#' after the conditions of interest in an experiment and it returns a contrast
#' matrix in a format acceptable to MSstats.
#'
#' @param numerators Character list of conditions which are the numerators of a
#'  series of a/b comparisons.
#' @param denominators Character list of conditions which are the denominators of a
#'  series of a/b comparisons.
#' @return Contrast matrix suitable for use in tools like MSstats.
#' @seealso [MSstats]
#' @export
make_simplified_contrast_matrix <- function(numerators, denominators) {
  if (length(numerators) != length(denominators)) {
    stop("Need a constant number of numerators and denominators.")
  }
  conditions <- unique(c(numerators, denominators))
  contrasts <- matrix(nrow = length(numerators), ncol = length(conditions))
  colnames(contrasts) <- conditions
  rownames(contrasts) <- numerators
  for (col in (colnames(contrasts))) {
    contrasts[, col] <- 0
  }
  for (n in seq_along(numerators)) {
    num <- numerators[[n]]
    den <- denominators[[n]]
    cont <- glue("{num}_vs_{den}")
    contrasts[n, num] <- 1
    contrasts[n, den] <- -1
    rownames(contrasts)[n] <- cont
  }
  return(contrasts)
}

#' Implement the arescan function in R
#'
#' This function was taken almost verbatim from AREScore() in SeqTools
#' Available at: https://github.com/lianos/seqtools.git
#' At least on my computer I could not make that implementation work
#' So I rewrapped its apply() calls and am now hoping to extend its logic
#' a little to make it more sensitive and get rid of some of the spurious
#' parameters or at least make them more transparent.
#'
#' @param x DNA/RNA StringSet containing the UTR sequences of interest
#' @param basal I dunno.
#' @param overlapping default = 1.5
#' @param d1.3 default = 0.75  These parameter names are so stupid, lets be realistic
#' @param d4.6 default = 0.4
#' @param d7.9 default = 0.2
#' @param within.AU default = 0.3
#' @param aub.min.length default = 10
#' @param aub.p.to.start default = 0.8
#' @param aub.p.to.end default = 0.55
#' @return a DataFrame of scores
#' @seealso [IRanges] [Biostrings] [GenomicRanges]
#' @examples
#' \dontrun{
#'  ## Extract all the genes from my genome, pull a static region 120nt following the stop
#'  ## and test them for potential ARE sequences.
#'  ## FIXME: There may be an error in this example, another version I have
#'  ## handles the +/- strand genes separately, I need to return to this and check
#'  ## if it is providing the 5' UTR for 1/2 the genome, which would be
#'  ## unfortunate -- but the logic for testing remains the same.
#'  are_candidates <- hpgl_arescore(genome)
#'  utr_genes <- subset(lmajor_annotations, type == 'gene')
#'  threep <- GenomicRanges::GRanges(seqnames = Rle(utr_genes[,1]),
#'                                   ranges = IRanges(utr_genes[,3], end=(utr_genes[,3] + 120)),
#'                                   strand = Rle(utr_genes[,5]),
#'                                   name = Rle(utr_genes[,10]))
#'  threep_seqstrings <- Biostrings::getSeq(lm, threep)
#'  are_test <- hpgltools::hpgl_arescore(x = threep_seqstrings)
#'  are_genes <- rownames(are_test[ which(are_test$score > 0), ])
#' }
#' @export
hpgl_arescore <- function(x, basal = 1, overlapping = 1.5, d1.3=0.75, d4.6=0.4,
                          d7.9=0.2, within.AU = 0.3, aub.min.length = 10, aub.p.to.start = 0.8,
                          aub.p.to.end = 0.55) {
  xtype <- match.arg(substr(class(x), 1, 3), c("DNA", "RNA"))
  if (xtype == "DNA") {
    pentamer <- "ATTTA"
    overmer <- "ATTTATTTA"
  } else {
    pentamer <- "AUUUA"
    overmer <- "AUUUAUUUA"
  }
  x <- as(x, "DNAStringSet")
  pmatches <- Biostrings::vmatchPattern(pentamer, x)
  omatches <- Biostrings::vmatchPattern(overmer, x)
  basal.score <- S4Vectors::elementNROWS(pmatches) * basal
  over.score <- S4Vectors::elementNROWS(omatches) * overlapping
  no.cluster <- data.frame(d1.3 = 0, d4.6 = 0, d7.9 = 0)
  clust <- lapply(pmatches, function(m) {
    if (length(m) < 2) {
      return(no.cluster)
    }
    wg <- BiocGenerics::width(IRanges::gaps(m))
    data.frame(d1.3=sum(wg <= 3), d4.6=sum(wg >= 4 & wg <= 6), d7.9=sum(wg >= 7 & wg <= 9))
  })
  clust <- do.call(rbind, clust)
  dscores <- clust$d1.3 * d1.3 + clust$d4.6 * d4.6 + clust$d7.9 *  d7.9
  ## require.auto("Biostrings")
  au.blocks <- my_identifyAUBlocks(x, min.length = aub.min.length,
                                   p.to.start = aub.p.to.start, p.to.end = aub.p.to.end)
  aub.score <- sum(IRanges::countOverlaps(pmatches, au.blocks) * within.AU)
  score <- basal.score + over.score + dscores + aub.score
  ans <- S4Vectors::DataFrame(score = score,
                              n.pentamer = S4Vectors::elementNROWS(pmatches),
                              n.overmer = S4Vectors::elementNROWS(omatches),
                              au.blocks = au.blocks,
                              n.au.blocks = S4Vectors::elementNROWS(au.blocks))
  cbind(ans, S4Vectors::DataFrame(clust))
}

#' Wrap cor() to include robust correlations.
#'
#' Take covRob's robust correlation coefficient and add it to the set of
#' correlations available when one calls cor().  I should reimplement this using
#' S4.
#'
#' @param df Data frame to test.
#' @param method Correlation method to use. Includes pearson, spearman, kendal, robust.
#' @param ... Other options to pass to stats::cor().
#' @return Some fun correlation statistics.
#' @seealso [robust]
#' @examples
#' \dontrun{
#'  hpgl_cor(df = df)
#'  hpgl_cor(df = df, method = "robust")
#' }
#' @export
hpgl_cor <- function(df, method = "pearson", ...) {
  if (method == "robust") {
    robust_cov <- robust::covRob(df, corr = TRUE)
    correlation <- robust_cov[["cov"]]
  } else if (method == "cordist") {
    correlation <- cordist(df, ...)
  } else {
    correlation <- stats::cor(df, method = method, ...)
  }
  return(correlation)
}

#' Wrap p.adjust to add IHW adjustments as an option.
#'
#' IHW and apeglm are the two new toys I found, this adds the former as a way to
#' adjust p-values.
#'
#' @param data Column or table containing values to adjust.
#' @param pvalue_column Name of the column in a table containing the p-values.
#' @param mean_column Name of the column in a table containing the mean count
#'   values to weight.
#' @param method p adjustment method to apply.
#' @param significance Passed to IHW
#' @param type Assuming a DE table, what type of DE is this?
#' @return Newly adjusted p-values using either p.adjust() or IHW.
#' @seealso [IHW]
#' @export
hpgl_padjust <- function(data, pvalue_column = "pvalue", mean_column = "base_mean",
                         method = "fdr", significance = 0.05, type = NULL) {
  if (method == "ihw") {
    result <- ihw_adjust(data, pvalue_column = pvalue_column, type = type,
                         mean_column = mean_column, significance = significance)
  } else {
    result <- p.adjust(data[[pvalue_column]], method = method)
  }
  return(result)
}

#' Because I am not smart enough to remember t()
#'
#' It seems to me there should be a function as easy for distances are there is for correlations.
#'
#' @param df data frame from which to calculate distances.
#' @param method Which distance calculation to use?
#' @param ... Extra arguments for dist.
#' @export
hpgl_dist <- function(df, method = "euclidean", ...) {
  input <- t(as.matrix(df))
  result <- as.matrix(dist(input, method = method, ...))
  return(result)
}

#' Load a backup rdata file
#'
#' I often use R over a sshfs connection, sometimes with significant latency,
#' and I want to be able to save/load my R sessions relatively quickly. Thus
#' this function uses my backup directory to load its R environment.
#'
#' @param directory Directory containing the RData.rda.xz file.
#' @param filename Filename to which to save.
#' @return a bigger global environment
#' @seealso [saveme()]
#' @examples
#' \dontrun{
#'  loadme()
#' }
#' @export
loadme <- function(directory = "savefiles", filename = "Rdata.rda.xz") {
  savefile <- glue("{getwd()}/{directory}/{filename}")
  message("Loading the savefile: ", savefile)
  load_string <- glue("load('{savefile}', envir = globalenv())")
  message("Command run: ", load_string)
  eval(parse(text = load_string))
}

#' Perform a get_value for delimited files
#'
#' Keith wrote this as .get_value() but functions which start with . trouble me.
#'
#' @param x Some stuff to split
#' @param delimiter The tritrypdb uses ': ' ergo the default.
#' @return A value!
local_get_value <- function(x, delimiter = ": ") {
  return(gsub(pattern = "^ ", replacement = "",
              x = tail(unlist(strsplit(x, delimiter)), n = 1), fixed = TRUE))
}

#' copy/paste the function from SeqTools and figure out where it falls on its ass.
#'
#' Yeah, I do not remember what I changed in this function.
#'
#' @param x Sequence object
#' @param min.length I dunno.
#' @param p.to.start P to start of course
#' @param p.to.end The p to end -- wtf who makes names like this?
#' @return a list of IRanges which contain a bunch of As and Us.
my_identifyAUBlocks <- function (x, min.length = 20, p.to.start = 0.8, p.to.end = 0.55) {
  xtype <- match.arg(substr(class(x), 1, 3), c("DNA", "RNA"))
  stopifnot(S4Vectors::isSingleNumber(min.length) & min.length >= 5 &  min.length <= 50)
  stopifnot(S4Vectors::isSingleNumber(p.to.start) & p.to.start >= 0.5 & p.to.start <= 0.95)
  stopifnot(S4Vectors::isSingleNumber(p.to.end) & p.to.end >= 0.2 & p.to.end <= 0.7)
  stopifnot(p.to.start > p.to.end)
  if (xtype == "DNA") {
    AU <- "AT"
  } else {
    AU <- "AU"
  }
  y <- as(x, sprintf("%sStringSet", xtype))
  widths <- BiocGenerics::width(x)

  ## I don't want to steal lianos' C from his repository.
  test_seqtools <- "SeqTools" %in% installed.packages()
  if (!isTRUE(test_seqtools)) {
    message("Installing lianos/seqtools from github to get the find_au_start_end() C function.")
    test <- devtools::install_github("lianos/seqtools/R/pkg")
  }
  lib_result <- sm(requireNamespace("SeqTools"))
  att_result <- sm(try(attachNamespace("SeqTools"), silent = TRUE))
  ## library(SeqTools)
  fun <- function(i) {
    one_seq <- x[[i]]
    au <- Biostrings::letterFrequencyInSlidingView(one_seq, min.length, AU, as.prob = TRUE)
    if (is.null(au) | nrow(au) == 0) {
      return(IRanges::IRanges())
    }
    au <- as.numeric(au)
    can.start <- au >= p.to.start
    can.end <- au <= p.to.end
    posts <- .Call("find_au_start_end", au, p.to.start, p.to.end, PACKAGE = "SeqTools")
    blocks <- IRanges::IRanges(posts$start, posts$end + min.length -  1L)
    IRanges::end(blocks) <- ifelse(IRanges::end(blocks) > widths[i],
                                   widths[i],
                                   IRanges::end(blocks))
    IRanges::reduce(blocks)
  }
  au.blocks <- lapply(1:length(x), fun)
  ret <- IRanges::IRangesList(au.blocks)
  return(ret)
}

#' Automatic loading and/or installing of packages.
#'
#' Load a library, install it first if necessary.
#'
#' This was taken from:
#' http://sbamin.com/2012/11/05/tips-for-working-in-r-automatically-install-missing-package/
#' and initially provided by Ramzi Temanni.
#'
#' @param lib String name of a library to check/install.
#' @param update Update packages?
#' @return 0 or 1, whether a package was installed or not.
#' @seealso [BiocManager] [install.packages()]
#' @examples
#' \dontrun{
#'  require.auto("ggplot2")
#' }
#' @export
please_install <- function(lib, update = FALSE) {
  count <- 0
  local({
    r <- getOption("repos")
    r["CRAN"] <- "http://cran.r-project.org"
    options(repos = r)
  })
  if (isTRUE(update)) {
    utils::update.packages(ask = FALSE)
  }
  github_path <- NULL
  ## If there is a / in the library's name, assume it is a github path
  split_lib <- strsplit(x = lib, split = "/")[[1]]
  if (length(split_lib) == 2) {
    github_path <- lib
    lib <- split_lib[[2]]
  }
  if (!isTRUE(lib %in% .packages(all.available = TRUE))) {
    if (is.null(github_path)) {
      ##eval(parse(text = paste("BiocManager::install('", lib, "', update = FALSE)", sep = "")))
      installed <- BiocManager::install(lib, update = FALSE)
      count <- 1
    } else {
      ret <- try(devtools::install_github(github_path, upgrade = "never"))
      if (class(ret) != "try-error") {
        count <- 1
      }
    }
  }
  return(count)
}

#' Send the R plotter to the computer of your choice!
#'
#' Resets the display and xauthority variables to the new computer I am using so
#' that plot() works.
#'
#' @param display DISPLAY variable to use, if NULL it looks in ~/.displays/$(host).last
#' @return Fresh plotting window to the display of your choice!
#' @export
rex <- function(display = ":0") {
  if (!is.null(dev.list())) {
    for (i in seq_along(dev.list())) {
      dev.off()
    }
  }
  home <- Sys.getenv("HOME")
  host <- Sys.info()[["nodename"]]
  if (is.null(display)) {
    display <- read.table(file.path(home, ".displays", glue("{host}.last")))[1, 1]
  }
  auth <- file.path(home, ".Xauthority")
  message("Setting display to: ", display)
  result <- Sys.setenv("DISPLAY" = display, "XAUTHORITY" = auth)
  X11(display = display)
  return(NULL)
}

#' Add a little logic to rmarkdown::render to date the final outputs as per a
#' request from Najib.
#'
#' @param file Rmd file to render.
#' @param format Chosen file format.
#' @return Final filename including the prefix rundate.
#' @seealso [rmarkdown]
#' @export
renderme <- function(file, format = "html_document", overwrite = TRUE) {
  original <- Sys.getenv("TMPDIR")
  new_tmpdir <- paste0("render_tmp_", file)
  if (file.exists(new_tmpdir)) {
    if (isTRUE(overwrite)) {
      message("The TMPDIR: ", new_tmpdir, " exists, removing and recreating it.")
      removed <- unlink(new_tmpdir, recursive = TRUE)
      recreated <- dir.create(new_tmpdir)
    } else {
      message("The TMPDIR: ", new_tmpdir, " exists, leaving it alone.")
    }
  }
  Sys.setenv("TMPDIR" = new_tmpdir)
  ret <- rmarkdown::render(file, output_format = format, envir = globalenv())
  removed <- unlink(new_tmpdir, recursive = TRUE)
  Sys.setenv("TMPDIR" = original)
  rundate <- format(Sys.Date(), format = "%Y%m%d")
  outdir <- dirname(ret)
  base <- basename(ret)
  b <- tools::file_path_sans_ext(base)
  ex <- tools::file_ext(ret)
  from <- file.path(outdir, glue::glue("{b}.{ex}"))
  to <- file.path(outdir, glue::glue("{rundate}_{b}.{ex}"))
  message("Moving ", from, " to\n", basename(to), ".")
  final <- file.rename(from, to)
  return(ret)
}

#' Make a backup rdata file for future reference
#'
#' I often use R over a sshfs connection, sometimes with significant latency, and
#' I want to be able to save/load my R sessions relatively quickly.
#' Thus this function uses pxz to compress the R session maximally and relatively fast.
#' This assumes you have pxz installed and >= 4 CPUs.
#'
#' @param directory  Directory to save the Rdata file.
#' @param backups  How many revisions?
#' @param cpus  How many cpus to use for the xz call
#' @param filename  Choose a filename.
#' @return Command string used to save the global environment.
#' @seealso [loadme()]
#' @examples
#' \dontrun{
#'  saveme()
#' }
#' @export
saveme <- function(directory = "savefiles", backups = 2, cpus = 6, filename = "Rdata.rda.xz") {
  environment()
  if (!file.exists(directory)) {
    dir.create(directory)
  }
  ##savefile <- glue("{getwd()}/{directory}/{filename}")
  savefile <- file.path(getwd(), directory, filename)
  message("The savefile is: ", savefile)
  backup_file(savefile, backups = backups)
  ## The following save strings work:
  save_string <- glue(
    "con <- pipe(paste0('pxz > {savefile}'), 'wb'); \\
    save(list = ls(all.names = TRUE, envir = globalenv()),
         envir = globalenv(), file = con, compress = FALSE); \\
    close(con)")
  message("The save string is: ", save_string)
  eval(parse(text = save_string))
}

#' Calculate a simplistic distance function of a point against two axes.
#'
#' Sillydist provides a distance of any point vs. the axes of a plot.
#' This just takes the abs(distances) of each point to the axes,
#' normalizes them against the largest point on the axes, multiplies
#' the result, and normalizes against the max of all point.
#'
#' @param firstterm X-values of the points.
#' @param secondterm Y-values of the points.
#' @param firstaxis X-value of the vertical axis.
#' @param secondaxis Y-value of the second axis.
#' @return Dataframe of the distances.
#' @seealso [ggplot2]
#' @examples
#' \dontrun{
#'  mydist <- sillydist(df[,1], df[,2], first_median, second_median)
#'  first_vs_second <- ggplot(df, aes(x = .data[["first"]], y = .data[["second"]]),
#'                                     environment = hpgl_env) +
#'    ggplot2::xlab(paste("Expression of", df_x_axis)) +
#'    ggplot2::ylab(paste("Expression of", df_y_axis)) +
#'    ggplot2::geom_vline(color = "grey", xintercept=(first_median - first_mad), size = line_size) +
#'    ggplot2::geom_vline(color = "grey", xintercept=(first_median + first_mad), size = line_size) +
#'    ggplot2::geom_vline(color = "darkgrey", xintercept = first_median, size = line_size) +
#'    ggplot2::geom_hline(color = "grey", yintercept=(second_median - second_mad), size = line_size) +
#'    ggplot2::geom_hline(color = "grey", yintercept=(second_median + second_mad), size = line_size) +
#'    ggplot2::geom_hline(color = "darkgrey", yintercept = second_median, size = line_size) +
#'    ggplot2::geom_point(colour = grDevices::hsv(mydist$dist, 1, mydist$dist),
#'                        alpha = 0.6, size = size) +
#'    ggplot2::theme(legend.position = "none")
#'  first_vs_second  ## dots get colored according to how far they are from the medians
#'  ## replace first_median, second_median with 0,0 for the axes
#' }
#' @export
sillydist <- function(firstterm, secondterm, firstaxis = 0, secondaxis = 0) {
  dataframe <- data.frame(firstterm, secondterm)
  dataframe[["x"]] <- (abs(dataframe[, 1]) - abs(firstaxis)) / abs(firstaxis)
  dataframe[["y"]] <- abs((dataframe[, 2] - secondaxis) / secondaxis)
  dataframe[["x"]] <- abs(dataframe[, 1] / max(dataframe[["x"]]))
  dataframe[["y"]] <- abs(dataframe[, 2] / max(dataframe[["y"]]))
  dataframe[["dist"]] <- abs(dataframe$x * dataframe[["y"]])
  dataframe[["dist"]] <- dataframe$dist / max(dataframe[["dist"]])
  return(dataframe)
}

#' Silence
#'
#' Some libraries/functions just won't shut up.  Ergo, silence, peasant!
#' This is a simpler silence peasant.
#'
#' @param ... Some code to shut up.
#' @param wrap  Wrap the invocation and try again if it failed?
#' @return Whatever the code would have returned.
#' @export
sm <- function(..., wrap = TRUE) {
  ret <- NULL
  output <- capture.output(type = "output", {
    if (isTRUE(wrap)) {
      ret <- try(suppressWarnings(suppressMessages(...)), silent = TRUE)
      if (class(ret)[1] == "try-error") {
        if (grepl(pattern = " there is no package called", x = ret)) {
          uninstalled <- trimws(gsub(pattern = "^.* there is no package called '(.*)'.*$",
                                     replacement = "\\1",
                                     x = ret, perl = TRUE))
          message("Going to attempt to install: ", uninstalled)
          tt <- please_install(uninstalled)
        }
        ret <- sm(..., wrap = FALSE)
      }
    } else {
      ret <- suppressWarnings(suppressMessages(...))
    }
  })
  return(ret)
}

#' Some ggplot2 stats functions have not yet implemented the new dropped_aes
#' flag, and it is driving me nuts.  Hopefully this will make it less frustrating for me.
#'
#' @param ... The set of functions to silence.
#' @export
sp <- function(...) {
  suppressWarnings(print(...))
}

#' A hopefully more robust version of tempfile().
#'
#' @param pattern Filename prefix.
#' @param suffix Filename suffix.
#' @param digits Currently I use Sys.time() with this number of digits.
#' @param body No implemented, intended to use other sources of digest()
#' @param fileext Filename extension as per tempfile().
#' @return md5 based tempfilename.
#' @export
tmpmd5file <- function(pattern = "", suffix = "", digits = 6,
                       body = NULL, fileext = "") {
  if (!grepl(pattern = "^\\.", x = fileext)) {
    pattern <- paste0(".", pattern)
  }
  op <- options(digits.secs = digits)
  body_string <- digest::digest(Sys.time())
  new <- options(op)
  outdir <- tempdir()
  if (!file.exists(outdir)) {
    created <- dir.create(outdir, recursive = TRUE)
  }
  file_string <- paste0(pattern, body_string, suffix, fileext)
  file_path <- file.path(outdir, file_string)
  return(file_path)
}

#' Remove the AsIs attribute from some data structure.
#'
#' Notably, when using some gene ontology libraries, the returned data
#' structures include information which is set to type 'AsIs' which turns out to
#' be more than slightly difficult to work with.
#'
#' @param stuff  The data from which to remove the AsIs classification.
#' @export
unAsIs <- function(stuff) {
  if("AsIs" %in% class(stuff)) {
    class(stuff) <- class(stuff)[-match("AsIs", class(stuff))]
  }
  return(stuff)
}

#' Print a model as y = mx + b just like in grade school!
#'
#' Because, why not!?
#'
#' @param lm_model Model to print from glm/lm/robustbase.
#' @param as Type to return.
#' @return a string representation of that model.
#' @export
ymxb_print <- function(lm_model, as = "glue") {
  coefficients <- summary(lm_model)[["coefficients"]]
  int <- signif(x = coefficients["(Intercept)", 1], digits = 3)
  ## Some invocations of lm() provide different row/column names.
  m <- signif(x = coefficients[2, 1], digits = 3)
  ret <- NULL
  if (as != "glue") {
    retlst <- list(
      "m" = m,
      "int" = int)
    return(retlst)
  }
  if (as.numeric(int) >= 0) {
    ret <- glue::glue("y = {m}x + {int}")
  } else {
    int <- int * -1
    ret <- glue::glue("y = {m}x - {int}")
  }
  return(ret)
}


## EOF
