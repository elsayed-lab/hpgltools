#' Clear an R session, this is probably unwise given what I have read about R.
#'
#' @param keepers  List of namespaces to leave alone (unimplemented).
#' @param depth  Cheesy forloop of attempts to remove packages stops after this many tries.
#' @return  A spring-fresh R session, hopefully.
#' @export
clear_session <- function(keepers=NULL, depth=10) {
  ## Partially taken from: https://stackoverflow.com/questions/7505547/detach-all-packages-while-working-in-r
  basic_packages <- c("package:stats", "package:graphics", "package:grDevices", "package:utils",
                      "package:datasets", "package:methods", "package:base")
  package_list <- search()[ifelse(unlist(gregexpr("package:", search()))==1, TRUE, FALSE)]
  package_list <- setdiff(package_list, basic_packages)
  result <- R.utils::gcDLLs()
  if (length(package_list) > 0) {
    for (package in package_list) {
      detach(package, character.only=TRUE)
    }
  }
  tt <- sm(rm(list=ls(all.names=TRUE), envir=globalenv()))
}

#' Get the current git commit for hpgltools
#'
#' @param gitdir  Directory containing the git repository.
#' @export
get_git_commit <- function(gitdir="/home/trey/hpgltools") {
  cmdline <- paste0("cd ", gitdir, " && git log -1 2>&1 | grep 'commit' | awk '{print $2}'")
  commit_result <- system(cmdline, intern=TRUE)
  cmdline <- paste0("cd ", gitdir, " && git log -1 2>&1 | grep 'Date' | perl -pe 's/^Date:\\s+//g'")
  date_result <- system(cmdline, intern=TRUE)
  result <- paste0(date_result, ": ", commit_result)
  message(paste0("If you wish to reproduce this exact build of hpgltools, invoke the following:"))
  message("> git clone http://github.com/abelew/hpgltools.git")
  message(paste0("> git reset ", commit_result))
  message("R> packrat::restore()")
  return(result)
}

#' Perform a get_value for delimited files
#'
#' Keith wrote this as .get_value() but functions which start with . trouble me.
#' @param x  Some stuff to split
#' @param delimiter  The tritrypdb uses ': ' ergo the default.
#' @return A value!
local_get_value <- function(x, delimiter=": ") {
  return(gsub("^ ", "", tail(unlist(strsplit(x, delimiter)), n=1), fixed=TRUE))
}

#' png() shortcut
#'
#' I hate remembering my options for png()
#'
#' @param file a filename to write
#' @param width  How wide?
#' @param height  How high?
#' @param res  The chosen resolution.
#' @param ...  Arguments passed to the image plotters.
#' @return a png/svg/eps/ps/pdf with height=width=9 inches and a high resolution
#' @export
pp <- function(file, width=9, height=9, res=180, ...) {
  ext <- tools::file_ext(file)
  res <- NULL
  if (ext == "png") {
    res <- png(filename=file, width=width, height=height, units="in", res=res, ...)
  } else if (ext == "svg") {
    res <- svg(filename=file)
  } else if (ext == "ps") {
    res <- postscript(file=file, width=width, height=height, units="in", ...)
  } else if (ext == "eps") {
    res <- cairo_ps(filename=file, width=width, height=height, ...)
  } else if (ext == "pdf") {
    res <- pdf(file=file, ...)
  } else {
    message("Defaulting to tiff.")
    res <- tiff(filename=file, width=width, height=height, units="in", ...)
  }
  return(res)
}

#' Silence, m...
#'
#' Some libraries/functions just won't shut up.  Ergo, silence, peasant!
#' This is a simpler silence peasant.
#'
#' @param ... Some code to shut up.
#' @return Whatever the code would have returned.
#' @export
sm <- function(...) {
  ret <- NULL
  output <- capture.output(type="output", {
    ret <- suppressWarnings(suppressMessages(...))
  })
  return(ret)
}

#' Make a knitr report with some defaults set a priori.
#'
#' I keep forgetting to set appropriate options for knitr.  This tries to set them.
#'
#' @param name Name the document!
#' @param type Html or pdf reports?
#' @return Dated report file.
#' @seealso \pkg{knitr} \pkg{rmarkdown} \pkg{knitrBootstrap}
#' @export
make_report <- function(name="report", type="pdf") {
  knitr::opts_knit$set(
                     progress = TRUE,
                     verbose = TRUE,
                     width = 90,
                     echo = TRUE)
  knitr::opts_chunk$set(
                      error = TRUE,
                      fig.width = 8,
                      fig.height = 8,
                      dpi = 96)
  options(digits = 4,
          stringsAsFactors = FALSE,
          knitr.duplicate.label = "allow")
  ggplot2::theme_set(ggplot2::theme_bw(base_size=10))
  set.seed(1)
  output_date <- format(Sys.time(), "%Y%m%d-%H%M")
  input_filename <- name
  ## In case I add .rmd on the end.
  input_filename <- gsub("\\.rmd", "", input_filename, perl=TRUE)
  input_filename <- gsub("\\.Rmd", "", input_filename, perl=TRUE)
  input_filename <- paste0(input_filename, ".Rmd")
  if (type == "html") {
    output_filename <- paste0(name, "-", output_date, ".html")
    output_format <- "html_document"
    rmarkdown::render(output_filename, output_format)
  } else {
    output_filename <- paste0(name, "-", output_date, ".pdf")
    output_format <- "pdf_document"
  }
  message(paste0("About to run: render(input=", input_filename, ", output_file=",
                 output_filename, " and output_format=", output_format))
  result <- try(rmarkdown::render(
                             "input" = input_filename,
                             "output_file" = output_filename,
                             "output_format" = output_format), silent=TRUE)

  return(result)
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
#' Note that I did this two months ago and haven't touched it since...
#'
#' @param x DNA/RNA StringSet containing the UTR sequences of interest
#' @param basal I dunno.
#' @param overlapping default=1.5
#' @param d1.3 default=0.75  These parameter names are so stupid, lets be realistic
#' @param d4.6 default=0.4
#' @param d7.9 default=0.2
#' @param within.AU default=0.3
#' @param aub.min.length default=10
#' @param aub.p.to.start default=0.8
#' @param aub.p.to.end default=0.55
#' @return a DataFrame of scores
#' @seealso \pkg{IRanges} \pkg{Biostrings}
#' @examples
#' \dontrun{
#'  ## Extract all the genes from my genome, pull a static region 120nt following the stop
#'  ## and test them for potential ARE sequences.
#'  ## FIXME: There may be an error in this example, another version I have handles the +/- strand
#'  ## genes separately, I need to return to this and check if it is providing the 5' UTR for 1/2
#'  ## the genome, which would be unfortunate -- but the logic for testing remains the same.
#'  are_candidates <- hpgl_arescore(genome)
#'  utr_genes <- subset(lmajor_annotations, type == 'gene')
#'  threep <- GenomicRanges::GRanges(seqnames=Rle(utr_genes[,1]),
#'                                   ranges=IRanges(utr_genes[,3], end=(utr_genes[,3] + 120)),
#'                                   strand=Rle(utr_genes[,5]),
#'                                   name=Rle(utr_genes[,10]))
#'  threep_seqstrings <- Biostrings::getSeq(lm, threep)
#'  are_test <- hpgltools:::hpgl_arescore(x=threep_seqstrings)
#'  are_genes <- rownames(are_test[ which(are_test$score > 0), ])
#' }
#' @export
hpgl_arescore <- function (x, basal=1, overlapping=1.5, d1.3=0.75, d4.6=0.4,
                           d7.9=0.2, within.AU=0.3, aub.min.length=10, aub.p.to.start=0.8,
                           aub.p.to.end=0.55) {
  ## The seqtools package I am using is called in R 'SeqTools' (note the capital S T)
  ## However, the repository I want for it is 'seqtools'
  ## Ergo my stupid require.auto() will be confused by definition because it assumes equivalent names
  ##if (isTRUE('SeqTools' %in% .packages(all.available=TRUE))) {
  ##    library('SeqTools')
  ##} else {
  ##    require.auto("lianos/seqtools/R/pkg")
  ##    library('SeqTools')
  ##}
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
  basal.score <- S4Vectors::elementLengths(pmatches) * basal
  over.score <- S4Vectors::elementLengths(omatches) * overlapping
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
  au.blocks <- my_identifyAUBlocks(x, aub.min.length, aub.p.to.start, aub.p.to.end)
  aub.score <- sum(IRanges::countOverlaps(pmatches, au.blocks) * within.AU)
  score <- basal.score + over.score + dscores + aub.score
  ans <- S4Vectors::DataFrame(score=score,
                              n.pentamer=S4Vectors::elementLengths(pmatches),
                              n.overmer=S4Vectors::elementLengths(omatches),
                              au.blocks=au.blocks,
                              n.au.blocks=S4Vectors::elementLengths(au.blocks))
  cbind(ans, S4Vectors::DataFrame(clust))
}

#' copy/paste the function from SeqTools and figure out where it falls on its ass.
#'
#' Yeah, I do not remember what I changed in this function.
#'
#' @param x Sequence object
#' @param min.length I dunno.
#' @param p.to.start P to start of course
#' @param p.to.end The p to end -- wtf who makes names like this?
#'
#' @return a list of IRanges which contain a bunch of As and Us.
my_identifyAUBlocks <- function (x, min.length=20, p.to.start=0.8, p.to.end=0.55) {
  xtype = match.arg(substr(class(x), 1, 3), c("DNA", "RNA"))
  stopifnot(S4Vectors::isSingleNumber(min.length) && min.length >= 5 &&  min.length <= 50)
  stopifnot(S4Vectors::isSingleNumber(p.to.start) && p.to.start >= 0.5 && p.to.start <= 0.95)
  stopifnot(S4Vectors::isSingleNumber(p.to.end) && p.to.end >= 0.2 && p.to.end <= 0.7)
  stopifnot(p.to.start > p.to.end)
  if (xtype == "DNA") {
    AU <- "AT"
  } else {
    AU <- "AU"
  }
  y <- as(x, sprintf("%sStringSet", xtype))

  widths <- BiocGenerics::width(x)
  fun <- function(i) {
    one_seq <- x[[i]]
    au <- Biostrings::letterFrequencyInSlidingView(one_seq, min.length, AU, as.prob=TRUE)
    if (is.null(au) | nrow(au) == 0) {
      return(IRanges::IRanges())
    }
    au <- as.numeric(au)
    can.start <- au >= p.to.start
    can.end <- au <= p.to.end
    posts <- .Call("find_au_start_end", au, p.to.start, p.to.end, PACKAGE="SeqTools")
    blocks <- IRanges::IRanges(posts$start, posts$end + min.length -  1L)
    stats::end(blocks) <- ifelse(stats::end(blocks) > widths[i], widths[i], stats::end(blocks))
    IRanges::reduce(blocks)
  }
  au.blocks = lapply(1:length(x), fun)
  ret <- IRanges::IRangesList(au.blocks)
  return(ret)
}

#' Wrap cor() to include robust correlations.
#'
#' Take covRob's robust correlation coefficient and add it to the set of correlations available when
#' one calls cor().
#'
#' @param df Data frame to test.
#' @param method Correlation method to use. Includes pearson, spearman, kendal, robust.
#' @param ... Other options to pass to stats::cor().
#' @return Some fun correlation statistics.
#' @seealso \pkg{robust}
#'  \code{\link{cor}} \code{\link{cov}} \code{\link[robust]{covRob}}
#' @examples
#' \dontrun{
#'  hpgl_cor(df=df)
#'  hpgl_cor(df=df, method="robust")
#' }
#' @export
hpgl_cor <- function(df, method="pearson", ...) {
  if (method == "robust") {
    robust_cov <- robust::covRob(df, corr=TRUE)
    correlation <- robust_cov[["cov"]]
  } else {
    correlation <- stats::cor(df, method=method, ...)
  }
  return(correlation)
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
#' @seealso \pkg{ggplot2}
#' @examples
#' \dontrun{
#'  mydist <- sillydist(df[,1], df[,2], first_median, second_median)
#'  first_vs_second <- ggplot2::ggplot(df, ggplot2::aes_string(x="first", y="second"),
#'                                     environment=hpgl_env) +
#'    ggplot2::xlab(paste("Expression of", df_x_axis)) +
#'    ggplot2::ylab(paste("Expression of", df_y_axis)) +
#'    ggplot2::geom_vline(color="grey", xintercept=(first_median - first_mad), size=line_size) +
#'    ggplot2::geom_vline(color="grey", xintercept=(first_median + first_mad), size=line_size) +
#'    ggplot2::geom_vline(color="darkgrey", xintercept=first_median, size=line_size) +
#'    ggplot2::geom_hline(color="grey", yintercept=(second_median - second_mad), size=line_size) +
#'    ggplot2::geom_hline(color="grey", yintercept=(second_median + second_mad), size=line_size) +
#'    ggplot2::geom_hline(color="darkgrey", yintercept=second_median, size=line_size) +
#'    ggplot2::geom_point(colour=grDevices::hsv(mydist$dist, 1, mydist$dist),
#'                        alpha=0.6, size=size) +
#'    ggplot2::theme(legend.position="none")
#'  first_vs_second  ## dots get colored according to how far they are from the medians
#'  ## replace first_median, second_median with 0,0 for the axes
#' }
#' @export
sillydist <- function(firstterm, secondterm, firstaxis=0, secondaxis=0) {
  dataframe <- data.frame(firstterm, secondterm)
  dataframe[["x"]] <- (abs(dataframe[, 1]) - abs(firstaxis)) / abs(firstaxis)
  dataframe[["y"]] <- abs((dataframe[, 2] - secondaxis) / secondaxis)
  dataframe[["x"]] <- abs(dataframe[, 1] / max(dataframe[["x"]]))
  dataframe[["y"]] <- abs(dataframe[, 2] / max(dataframe[["y"]]))
  dataframe[["dist"]] <- abs(dataframe$x * dataframe[["y"]])
  dataframe[["dist"]] <- dataframe$dist / max(dataframe[["dist"]])
  return(dataframe)
}

#' Make a backup of an existing file with n revisions, like VMS!
#'
#' Sometimes I just want to kick myself for overwriting important files and then I remember using
#' VMS and wish modern computers were a little more like it.
#'
#' @param backup_file Filename to backup.
#' @param backups How many revisions?
backup_file <- function(backup_file, backups=4) {
  if (file.exists(backup_file)) {
    for (i in backups:01) {
      j <- i + 1
      i <- sprintf("%02d", i)
      j <- sprintf("%02d", j)
      test <- paste0(backup_file, ".", i)
      new <- paste0(backup_file, ".", j)
      if (file.exists(test)) {
        file.rename(test, new)
      }
    }
    newfile <- paste0(backup_file, ".", i)
    message(paste0("Renaming ", backup_file, " to ", newfile, "."))
    file.copy(backup_file, newfile)
  } else {
    message("The file does not yet exist.")
  }
}

#' Load a backup rdata file
#'
#' I often use R over a sshfs connection, sometimes with significant latency, and I want to be able
#' to save/load my R sessions relatively quickly. Thus this function uses my backup directory to
#' load its R environment.
#'
#' @param directory Directory containing the RData.rda.xz file.
#' @param filename  Filename to which to save.
#' @return a bigger global environment
#' @seealso \code{\link{saveme}} \code{\link{load}} \code{\link{save}}
#' @examples
#' \dontrun{
#'  loadme()
#' }
#' @export
loadme <- function(directory="savefiles", filename="Rdata.rda.xz") {
  savefile <- paste0(getwd(), "/", directory, "/", filename)
  message(paste0("Loading the savefile: ", savefile))
  load_string <- paste0("load('", savefile, "', envir=globalenv())")
  message(paste0("Command run: ", load_string))
  eval(parse(text=load_string))
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
#' @seealso \code{\link{save}} \code{\link{pipe}}
#' @examples
#' \dontrun{
#'  saveme()
#' }
#' @export
saveme <- function(directory="savefiles", backups=2, cpus=6, filename="Rdata.rda.xz") {
  environment()
  if (!file.exists(directory)) {
    dir.create(directory)
  }
  savefile <- paste0(getwd(), "/", directory, "/", filename)
  message(paste0("The savefile is: ", savefile))
  backup_file(savefile, backups=backups)
  ## The following save strings work:
  save_string <- paste0("con <- pipe(paste0('pxz -T", cpus, " > ",
                        savefile,
                        "'), 'wb');\n",
                        "save(list=ls(all.names=TRUE, envir=globalenv()), envir=globalenv(), file=con, compress=FALSE);\n",
                        "close(con)")
  message(paste0("The save string is: ", save_string))
  eval(parse(text=save_string))
}

#' Print a model as y = mx + b just like in grade school!
#'
#' Because, why not!?
#'
#' @param model Model to print from glm/lm/robustbase.
#' @return a string representation of that model.
#' @export
ymxb_print <- function(model) {
  intercept <- round(coefficients(model)[1], 2)
  x_name <- names(coefficients(model)[-1])
  slope <- round(coefficients(model)[-1], 2)
  ret <- paste0("y = ", slope, "*", x_name, " + ", intercept)
  message(ret)
  return(ret)
}

#' Resets the display and xauthority variables to the new computer I am using so that plot() works.
#'
#' This function assumes a line in the .profile which writes the DISPLAY variable to
#' ${HOME}/.displays/$(hostname).last
#'
#' @param display DISPLAY variable to use,  if NULL it looks in ~/.displays/$(host).last
#'
#' @export
rex <- function(display=":0") {
  home <- Sys.getenv("HOME")
  host <- Sys.info()[["nodename"]]
  if (is.null(display)) {
    display <- read.table(paste0(home, "/.displays/", host, ".last"))[1, 1]
  }
  auth <- paste0(home, "/.Xauthority")
  message(paste0("Setting display to: ", display))
  result <- Sys.setenv("DISPLAY" = display, "XAUTHORITY" = auth)
  X11(display=display)
  return(NULL)
}

## EOF
