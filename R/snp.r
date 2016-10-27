snp_add_file <- function(sample, input_dir="preprocessing/outputs", file_suffix="_parsed_ratio.txt") {
    tmp_dt <- as.data.table(read.table(paste0(input_dir, "/", sample, file_suffix)))
    rownames(tmp_dt) <- tmp_dt[["V1"]]
    tmp_dt <- tmp_dt[-1]
    colnames(tmp_dt) <- c("rownames", sample)
    snp_dt <- merge(snp_dt, tmp_dt, by="rownames", all.x=TRUE, all.y=TRUE)
    rownames(snp_dt) <- snp_dt[["Row.names"]]
    snp_dt <- snp_dt[-1]
    return(snp_dt)
}

#' Gather snp information for an expt
#'
#' I have some initial code for working with snps, but it seems that it will be getting more use, so
#' make it testable etc.
#'
#' @param expt an expressionset from which to extract information.
#' @return some stuff
#' @export
expt_snp <- function(expt,
                     input_dir="preprocessing/outputs",
                     file_suffix="_parsed_ratio.txt",
                     bam_suffix="_accepted_paired.bam",
                     tolower=TRUE) {
    expt <- parasite_expt
    samples <- rownames(Biobase::pData(expt$expressionset))
    if (isTRUE(tolower)) {
        samples <- tolower(samples)
    }
    sample <- samples[1]
    snp_dt <- as.data.table(read.table(paste0(input_dir, "/", sample, file_suffix)))
    rownames(snp_dt) <- snp_dt[["V1"]]
    snp_dt <- snp_dt[-1]
    colnames(snp_dt) <- c("rownames", sample)
    for (sample_num in 2:length(samples)) {
        sample <- samples[sample_num]
        message(paste0("Merging: ", sample))
        snp_dt <- snp_add_file(sample, input_dir=input_dir, file_suffix=file_suffix)
    }

    ## Now I have a data table of rownames and percentages
    ## Next I need to cross reference these against the coverage by position, ergo I must split out the rownames
    snp_dt[, c("species", "chromosome", "position", "original", "new") := tstrsplit(rownames, "_", fixed=TRUE)]

    snp_file_list <- function(row, the_samples=samples, indir=input_dir, suffix=bam_suffix) {
        filenames <- c()
        for (sample_num in 1:length(the_samples)) {
            sample <- the_samples[sample_num]
            filename <- paste0(indir, "/", sample, suffix)
            filenames <- append(filenames, filename)
        }
        return(filenames)
    }
    pileup_files <- file_list(samples)

    pileup_info <- Rsamtools::PileupFiles(pileup_files)
    ## Taken directly from the Rsamtools manual
    queries <- paste0(snp_dt[["species"]], "_",
                      snp_dt[["chromosome"]], ":",
                      as.numeric(snp_dt[["position"]]) - 1, "-",
                      as.numeric(snp_dt[["position"]]) + 1)
    which <- GenomicRanges::GRanges(queries)
    what <- "seq"
    pileups_params <- Rsamtools::ApplyPileupsParam(which=which, what=what)

    snp_calc_coverage <- function(x) {
        ## information at each pile-up position
        qme <- function(y) {
            y <- y[c("A", "C", "G", "T"), , drop=FALSE]
            y <- y + 1L
            result <- colSums(y)
            return(result)
        }
        info <- apply(x[["seq"]], 2, qme)
        retlist <- list(seqnames=x[["seqnames"]], pos=x[["pos"]], info=info)
        return(retlist)
    }
    message("Gathering coverage for all bamfiles, this takes some time.")
    coverage_result <- Rsamtools::applyPileups(pileup_info, snp_calc_coverage, param=pileups_params)
    ## result is a list of n elements where n is the number of rows in snp_dt
    ## Each element of result is in turn a list containing the following slots:
    ##  seqnames (chromosome), pos (position(s)), info (coverage by file)
    ## The piece of information we want to put into snp_dt is therefore:
    ## coverage_list <- result[[snp_dt_row]][[info]][2, ]
    ## coverage_list is in turn a character list named by filename (which begins with the sample ID)
    ## We will therefore extract the hpglID from it and the coverage for every sample
   names(coverage_result) <- snp_dt[["rownames"]]

    ## Now extract from the rather strange coverage_result data the coverage by position/sample
    new_dt <- NULL
    snp_extract_coverage <- function(element) {
        row <- element[["info"]][2, ]
        names(row) <- samples
        new_dt <<- rbind(new_dt, row)
    }
    ## unused_var <- try(lapply(coverage_result, snp_extract_coverage))
    unused_var <- try(multicore::mclapply(coverage_result, snp_extract_coverage))
    if (class(unused_var) == "try-error") {
        message("There was an error when creating the data table of coverage, hopefully the data is salvageable.")
    }
    new_dt <- as.data.table(new_dt)
    new_dt[["rownames"]] <- snp_dt[["rownames"]]

}
