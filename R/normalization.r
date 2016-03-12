## Time-stamp: <Thu Mar 10 16:56:21 2016 Ashton Trey Belew (abelew@gmail.com)>

## Note to self, @title and @description are not needed in roxygen
## comments, the first separate #' is the title, the second the
## description, the third is the long form description.  Then add in
## the @param @return @seealso @export etc...

## My root question, to which I think I have a small idea about the answer:
## Two of the very many paths to toptable()/toptags():
##  1.  normalize data -> model(~0 + condition + batch) -> limma
###     Including batch in the model loses some power, but improves the 'truth' of the result
##  2.  normalize data -> batch correction(factor(batch)) ->  model(~0 + condition) -> limma
###     Power lost before the model, also improves the 'truth' of the result.
##  Why is #1 better than #2?
### More well understood and conservative.
##  Why have we nonetheless done #2 in a few instances?  (not only because we learned that first)

#' A function suggested by Hector Corrada Bravo and Kwame Okrah for batch removal
#'
#' During a lab meeting, the following function was suggested as a quick and dirty batch removal tool
#'
#' @param normalized_counts  a data frame of log2cpm counts
#' @param model  a balanced experimental model containing condition and batch factors
#' @return a dataframe of residuals after subtracting batch from the model
#' @seealso \link[limma]{voom} \link[limma]{lmFit}
#' @examples
#' \dontrun{
#' newdata <- cbcb_batch_effect(counts, expt_model)
#' }
#' @export
cbcb_batch_effect <- function(normalized_counts, model) {
    ## model = model.matrix(~ condition + batch)
    voomed <- hpgl_voom(normalized_counts, model)
    voomed_fit <- limma::lmFit(voomed)
    modified_model <- model
    modified_model <- modified_model[, grep("batch", colnames(modified_model))] <- 0 ## Drop batch from the model
    new_data <- tcrossprod(voomed_fit$coefficient, modified_model) + residuals(voomed_fit, normalized_counts)
    return(new_data)
}

#'   Perform different batch corrections using limma, sva, ruvg, and cbcbSEQ.
#'   I found this note which is the clearest explanation of what happens with batch effect data:
#'   https://support.bioconductor.org/p/76099/
#'   Just to be clear, there's an important difference between removing a batch effect and modelling a
#'   batch effect. Including the batch in your design formula will model the batch effect in the
#'   regression step, which means that the raw data are not modified (so the batch effect is not
#'   removed), but instead the regression will estimate the size of the batch effect and subtract it
#'   out when performing all other tests. In addition, the model's residual degrees of freedom will
#'   be reduced appropriately to reflect the fact that some degrees of freedom were "spent"
#'   modelling the batch effects. This is the preferred approach for any method that is capable of
#'   using it (this includes DESeq2). You would only remove the batch effect (e.g. using limma's
#'   removeBatchEffect function) if you were going to do some kind of downstream analysis that can't
#'   model the batch effects, such as training a classifier.
#'   I don't have experience with ComBat, but I would expect that you run it on log-transformed CPM
#'   values, while DESeq2 expects raw counts as input. I couldn't tell you how to properly use the
#'   two methods together.
#'
#' @param count_table  a matrix of (pseudo)counts.
#' @param design  a model matrix defining the experimental conditions/batches/etc
#' @param batch   a string describing the method to try to remove the batch effect (or FALSE to leave it alone, TRUE uses limma)
#' @param batch1   the column in the design table describing the presumed covariant to remove.
#' @param batch2   the column in the design table describing the second covariant to remove (only used by limma at the moment).
#' @param noscale   used for combatmod, when true it removes the scaling parameter from the invocation of the modified combat.
#' @param ... more options for you!
#' @return The 'batch corrected' count table and new library size.  Please remember that the library size which comes out of this
#' may not be what you want for voom/limma and would therefore lead to spurious differential expression values.
#' @seealso \pkg{limma} \pkg{edgeR} \pkg{RUVSeq} \pkg{sva} \pkg{cbcbSEQ}
#' @examples
#' \dontrun{
#' limma_batch <- batch_counts(table, design, batch1='batch', batch2='strain')
#' sva_batch <- batch_counts(table, design, batch='sva')
#' }
#' @export
batch_counts <- function(count_table, design, batch=TRUE, batch1='batch', batch2=NULL , noscale=TRUE, ...) {
    batches <- as.factor(design[, batch1])
    conditions <- as.factor(design[, "condition"])

    num_low <- sum(count_table < 1 & count_table > 0)
    if (num_low > 0) {
        message(paste0("batch_counts: Before batch correction, ", num_low, " entries 0<x<1."))
    }
    num_zero <- sum(count_table <= 0)
    if (num_zero > 0) {
        message(paste0("batch_counts: Before batch correction, ", num_zero, " entries are >= 0."))
    }
    if (isTRUE(batch)) {
        batch <- "limma"
    }
    if (batch == "limma") {
        if (is.null(batch2)) {
            ## A reminder of removeBatchEffect usage
            ## adjusted_batchdonor = removeBatchEffect(data, batch=as.factor(as.character(des$donor)), batch2=as.factor(as.character(des$batch)))
            message("batch_counts: Using limma's removeBatchEffect to remove batch effect.")
            count_table <- limma::removeBatchEffect(count_table, batch=batches)
        } else {
            batches2 <- as.factor(design[, batch2])
            count_table <- limma::removeBatchEffect(count_table, batch=batches, batch2=batches2)
        }
    } else if (batch == 'limmaresid') {
        message("batch_counts: Using residuals of limma's lmfit to remove batch effect.")
        batch_model <- model.matrix(~batches)
        batch_voom <- limma::voom(data.frame(count_table), batch_model, normalize.method="quantile", plot=FALSE)
        batch_fit <- limma::lmFit(batch_voom, design=batch_model)
        count_table <- residuals(batch_fit, batch_voom$E)
    } else if (batch == "combatmod") {
        ## normalized_data = hpgl_combatMod(dat=data.frame(counts), batch=batches, mod=conditions, noScale=noscale, ...)
        message("batch_counts: Using a modified cbcbSEQ combatMod for batch correction.")
        count_table <- hpgl_combatMod(dat=data.frame(count_table), batch=batches, mod=conditions, noScale=noscale, ...)
    } else if (batch == "sva") {
        message("batch_counts: Using sva::fsva for batch correction.")
        df <- data.frame(count_table)
        mtrx <- as.matrix(df)
        conditional_model <- model.matrix(~conditions, data=df)
        null_model <- conditional_model[,1]
        num_surrogates <- 0
        be_surrogates <- sva::num.sv(mtrx, conditional_model, method="be")
        leek_surrogates <- sva::num.sv(mtrx, conditional_model, method="leek")
        if (be_surrogates >= 1) {
            num_surrogates <- be_surrogates
        } else {
            num_surrogates <- leek_surrogates
        }
        sva_object <- sva::sva(mtrx, conditional_model, null_model, n.sv=num_surrogates)
        ## mod_sv = cbind(conditional_model, sva_object$sv)
        fsva_result <- sva::fsva(mtrx, conditional_model, sva_object, newdat=mtrx, method="exact")
        ## new_expt$conditional_model = conditional_model
        ## new_expt$null_model = null_model
        ## new_expt$num_surrogates = num_surrogates
        ## new_expt$sva_object = sva_object
        ## new_expt$mod_sv = mod_sv
        ## new_expt$fsva_result = fsva_result
        count_table <- fsva_result$db
    } else if (batch == 'combat_noprior') {
        message("batch_counts: Using sva::combat without a prior for batch correction.")
        count_table <- sva::ComBat(count_table, batches, mod=conditions, par.prior=FALSE, prior.plots=FALSE)
    } else if (batch == 'combat') {
        message("batch_counts: Using sva::combat with a prior for batch correction.")
        count_table <- sva::ComBat(count_table, batches, mod=conditions, par.prior=TRUE, prior.plots=TRUE)
    } else if (batch == "svaseq") {
        message("batch_counts: Using sva::svaseq for batch correction.")
        message("Note to self:  If you feed svaseq a data frame you will get an error like:")
        message("data %*% (Id - mod %*% blah blah requires numeric/complex arguments.")
        df <- data.frame(count_table)
        mtrx <- as.matrix(df)
        conditional_model <- model.matrix(~conditions, data=df)
        null_model <- conditional_model[,1]
        num_surrogates <- sva::num.sv(mtrx, conditional_model)
        svaseq_result <- sva::svaseq(mtrx, conditional_model, null_model, n.sv=num_surrogates)
        plot(svaseq_result$sv, pch=19, col="blue")
        ## The following was taken from: https://www.biostars.org/p/121489/
        X <- cbind(conditional_model, svaseq_result$sv)
        Hat <- solve(t(X) %*% X) %*% t(X)
        beta <- (Hat %*% t(mtrx))
        P <- ncol(conditional_model)
        count_table <- mtrx - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
    } else if (batch == "ruvg") {
        message("Using RUVSeq and edgeR for batch correction (similar to lmfit residuals.")
        ## Adapted from: http://jtleek.com/svaseq/simulateData.html -- but not quite correct yet
        ## As it stands I do not think this does anything useful
        ##require.auto("RUVSeq")
        conditional_model <- model.matrix(~conditions, data=df)
        y <- edgeR::DGEList(counts=count_table, group=conditions)
        y <- edgeR::calcNormFactors(y, method="upperquartile")
        y <- edgeR::estimateGLMCommonDisp(y, conditional_model)
        y <- edgeR::estimateGLMTagwiseDisp(y, conditional_model)
        fit <- edgeR::glmFit(y, conditional_model)
        lrt <- edgeR::glmLRT(fit, coef=2)
        controls <- rank(lrt$table$LR) <= 400
        batch_ruv_emp <- RUVSeq::RUVg(count_table, controls, k=1)$W
        X <- cbind(conditional_model, batch_ruv_emp)
        Hat <- solve(t(X) %*% X) %*% t(X)
        beta <- (Hat %*% t(mtrx))
        P <- ncol(conditional_model)
        count_table <- mtrx - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
    } else {
        message("Did not recognize the batch correction, leaving the table alone.")
        message("Recognized batch corrections include: 'limma', 'combatmod', 'sva',")
        message("limmaresid, combat_noprior, combat, svaseq, and ruvg.")
    }
    num_low <- sum(count_table < 0)
    if (num_low > 0) {
        message(paste0("The number of elements which are < 0 after batch correction is: ", num_low))
        count_table[count_table < 0] <- 0
    }
    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Filter low-count genes from a data set.
#'
#' This was a function written by Kwame Okrah and perhaps also Laura Dillon to remove low-count genes.  It drops genes based on a threshold and number of samples.
#'
#' @param count_table  a data frame of (pseudo)counts by sample.
#' @param threshold    lower threshold of counts for each gene.
#' @param min_samples    minimum number of samples
#' @param verbose   if set to true, prints number of genes removed and remaining.
#' @return dataframe of counts without the low-count genes
#' @seealso \link[cbcbSEQ]{log2CPM} which this uses to decide what to keep
#' @examples
#' \dontrun{
#' filtered_table <- cbcb_filter_counts(count_table)
#' }
#' @export
cbcb_filter_counts <- function(count_table, threshold=2, min_samples=2, verbose=FALSE) {
    ## I think having a log2cpm here is kind of weird, because the next step in processing is to cpm the data.
    ##cpms = 2^log2CPM(counts, lib.size=lib.size)$y
    ## cpms = 2^hpgl_log2cpm(counts)
    num_before <- nrow(count_table)

    if (class(count_table) == 'ExpressionSet') {
        keep <- rowSums(Biobase::exprs(count_table) > threshold) >= min_samples
    } else {
        keep <- rowSums(count_table > threshold) >= min_samples
    }

    count_table <- count_table[keep,]

    if (verbose) {
        message(sprintf("Removing %d low-count genes (%d remaining).",
                      num_before - nrow(count_table), nrow(count_table)))
    }

    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Perform a cpm/rpkm/whatever transformation of a count table.
#'
#' I should probably tell it to also handle a simple df/vector/list of gene lengths, but I haven't.
#' cp_seq_m is a cpm conversion of the data followed by a rp-ish
#' conversion which normalizes by the number of the given oligo.  By
#' default this oligo is 'TA' because it was used for tnseq which
#' should be normalized by the number of possible transposition sites
#' by mariner.  It could, however, be used to normalize by the number
#' of methionines, for example -- if one wanted to do such a thing.
#'
#' @param data A matrix of count data
#' @param convert   A type of conversion to perform: edgecpm/cpm/rpkm/cp_seq_m
#' @param annotations   a set of gff annotations are needed if using rpkm so we can get gene lengths.
#' @param fasta   a fasta for rpkmish
#' @param pattern   for cp_seq_m counts
#' @param entry_type  used to acquire gene lengths
#' @param ... more options
#' @return dataframe of cpm/rpkm/whatever(counts)
#' @seealso \pkg{edgeR} \pkg{Biobase} \code{\link[edgeR]{cpm}}
#' @examples
#' \dontrun{
#'  converted_table = convert_counts(count_table, convert='edgecpm')
#' }
#' @export
convert_counts <- function(data, convert="raw", annotations=NULL, fasta=NULL, pattern='TA', entry_type='gene', ...) {
    data_class <- class(data)[1]
    if (data_class == 'expt') {
        count_table <- Biobase::exprs(data$expressionset)
    } else if (data_class == 'ExpressionSet') {
        count_table <- Biobase::exprs(data)
    } else if (data_class == 'matrix' | data_class == 'data.frame') {
        count_table <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }
    if (convert == "edgecpm") {
        requireNamespace("edgeR")
        count_table <- edgeR::cpm(count_table)
    } else if (convert == "cpm") {
        lib_size <- colSums(count_table)
        ## count_table = t(t((count_table$counts + 0.5) / (lib_size + 1)) * 1e+06)
        transposed <- t(count_table + 0.5)
        cp_counts <- transposed / (lib_size + 1)
        cpm_counts <- t(cp_counts * 1e+06)
        count_table <- cpm_counts
    } else if (convert == "rpkm") {
        if (is.null(annotations)) {
            stop("RPKM conversion requires gene lengths.")
        }
        count_table <- hpgl_rpkm(count_table, annotations=annotations)
    } else if (convert == "cp_seq_m") {
        counts <- edgeR::cpm(count_table)
        ## count_table = divide_seq(counts, ...)
        count_table <- divide_seq(counts, fasta=fasta, gff=annotations, pattern=pattern, entry_type=entry_type)
    }
    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

#'   Express a data frame of counts as reads per pattern per
#' million(library).
#'
#' @param counts read count matrix
#' @param pattern pattern to search against.  Defaults to 'TA'
#' @param fasta a fasta genome to search
#' @param gff the gff set of annotations to define start/ends of genes.
#' @param entry_type which type of gff entry to search against.  Defaults to 'gene'.
#' @return The 'RPseqM' counts
#' @seealso \code{\link[Rsamtools]{FaFile}} \code{\link[edgeR]{rpkm}}
#' @examples
#' \dontrun{
#' cptam <- divide_seq(cont_table, fasta="mgas_5005.fasta.xz", gff="mgas_5005.gff.xz")
#' }
#' @export
divide_seq <- function(counts, pattern="TA", fasta="testme.fasta", gff="testme.gff", entry_type="gene") {
    if (!file.exists(fasta)) {
        compressed_fasta <- paste0(fasta, '.xz')
        system(paste0("xz -d ", compressed_fasta))
    }
    raw_seq <- try(Rsamtools::FaFile(fasta))
    if (class(raw_seq)[1] == 'try-error') {
        stop(paste0("There was a problem reading: ", fasta))
    }
    gff_entries <- gff2irange(gff)
    ## print(head(gff_entries))
    ##    cds_entries = subset(gff_entries, type==entry_type)
    found_entries <- (gff_entries$type == entry_type)
    if (sum(found_entries) == 0) {
        message(paste0("There were no found entries of type: ", entry_type, "."))
        message("Going to try locus_tag, and failing that, mRNA.")
        locus_entries <- (gff_entries$type == 'locus_tag')
        mrna_entries <- (gff_entries$type == 'mRNA')
        if ((sum(locus_entries) > sum(mrna_entries)) & sum(locus_entries) > 100) {
            found_entries <- locus_entries
        } else if (sum(mrna_entries) > 100) {
            found_entries <- mrna_entries
        } else {
            stop("Unable to find any entries of type locus_tag nor mrna.")
        }
    }
    ##cds_entries = subset(gff_entries, type==entry_type)
    cds_entries <- gff_entries[found_entries,]
    names(cds_entries) <- make.names(cds_entries$locus_tag, unique=TRUE)
    cds_seq <- Biostrings::getSeq(raw_seq, cds_entries)
    names(cds_seq) <- cds_entries$locus_tag
    dict <- Biostrings::PDict(pattern, max.mismatch=0)
    result <- Biostrings::vcountPDict(dict, cds_seq)
    num_tas <- data.frame(name=names(cds_seq), tas=as.data.frame(t(result)))
    rownames(num_tas) <- make.names(num_tas$name, unique=TRUE)
    colnames(num_tas) <- c("name","TAs")
    num_tas$TAs <- num_tas$TAs + 1
    factor <- median(num_tas$TAs)
    num_tas$TAs <- num_tas$TAs / factor
    merged_tas <- merge(counts, num_tas, by="row.names", all.x=TRUE)
    rownames(merged_tas) <- merged_tas$Row.names
    merged_tas <- merged_tas[-1]
    merged_tas <- subset(merged_tas, select=-c("name"))  ## Two different ways of removing columns...
    merged_tas <- merged_tas / merged_tas$TAs
    merged_tas <- merged_tas[, !(colnames(merged_tas) %in% c("TAs"))]  ## Here is another!
    return(merged_tas)
}

#' Filter low-count genes from a data set using genefilter's pOverA()
#'
#' I keep thinking this function is pofa... oh well.
#'
#' @param count_table  input data frame of counts by sample
#' @param p   a minimum proportion of each gene's counts/sample to be greater than a minimum(A)
#' @param A   the minimum number of counts in the above proportion
#' @param verbose   If set to true, prints number of genes removed / remaining
#' @return dataframe of counts without the low-count genes
#' @seealso \pkg{genefilter} \code{\link[genefilter]{pOverA}} which this uses to decide what to keep
#' @examples
#' \dontrun{
#'  filtered_table = genefilter_pofa_counts(count_table)
#' }
#' @export
genefilter_pofa_counts <- function(count_table, p=0.01, A=100, verbose=TRUE) {
    ## genefilter has functions to work with expressionsets directly, but I think I will work merely with tables in this.
    num_before <- nrow(count_table)

    if (class(count_table) == 'ExpressionSet') {
        counts <- Biobase::exprs(count_table)
    }
    test <- genefilter::pOverA(p=p, A=A)
    filter_list <- genefilter::filterfun(test)
    answer <- genefilter::genefilter(count_table, filter_list)
    count_table <- count_table[answer,]

    if (isTRUE(verbose)) {
        removed <- num_before - nrow(count_table)
        message(paste0("Removing ", removed, " low-count genes (", nrow(count_table), " remaining)."))
    }
    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Filter genes from a dataset outside a range of variance
#'
#' @param count_table  input data frame of counts by sample
#' @param cv_min   a minimum coefficient of variance
#' @param cv_max   guess
#' @param verbose   If set to true, prints number of genes removed / remaining
#' @return dataframe of counts without the low-count genes
#' @seealso \pkg{genefilter} \code{\link[genefilter]{kOverA}} which this uses to decide what to keep
#' @examples
#' \dontrun{
#' filtered_table = genefilter_kofa_counts(count_table)
#' }
#' @export
genefilter_cv_counts <- function(count_table, cv_min=0.01, cv_max=1000, verbose=FALSE) {
    ## genefilter has functions to work with expressionsets directly, but I think I will work merely with tables in this.
    num_before <- nrow(count_table)

    if (class(count_table) == 'ExpressionSet') {
        counts <- Biobase::exprs(count_table)
    }
    test <- genefilter::cv(cv_min, cv_max)
    filter_list <- genefilter::filterfun(test)
    answer <- genefilter::genefilter(count_table, filter_list)
    count_table <- count_table[answer,]

    if (verbose) {
        message(sprintf("Removing %d low-count genes (%d remaining).",
                      num_before - nrow(count_table), nrow(count_table)))
    }
    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Filter low-count genes from a data set using genefilter's kOverA()
#'
#' @param count_table input data frame of counts by sample
#' @param k   a minimum number of samples to have >A counts
#' @param A   the minimum number of counts for each gene's sample in kOverA()
#' @param verbose   If set to true, prints number of genes removed / remaining
#' @return dataframe of counts without the low-count genes
#' @seealso \pkg{genefilter} \code{\link[genefilter]{kOverA}} which this uses to decide what to keep
#' @examples
#' \dontrun{
#'  filtered_table = genefilter_kofa_counts(count_table)
#' }
#' @export
genefilter_kofa_counts <- function(count_table, k=1, A=1, verbose=FALSE) {
    ## genefilter has functions to work with expressionsets directly, but I think I will work merely with tables in this.
    num_before <- nrow(count_table)

    if (class(count_table) == 'ExpressionSet') {
        counts <- Biobase::exprs(count_table)
    }
    test <- genefilter::kOverA(k=k, A=A)
    filter_list <- genefilter::filterfun(test)
    answer <- genefilter::genefilter(count_table, filter_list)
    count_table <- count_table[answer,]

    if (verbose) {
        message(sprintf("Removing %d low-count genes (%d remaining).",
                      num_before - nrow(count_table), nrow(count_table)))
    }
    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Use a modified version of combat on some data
#' This is a hack of Kwame's combatMod to make it not fail on corner-cases.
#'
#' @param dat a df to modify
#' @param batch a factor of batches
#' @param mod a factor of conditions
#' @param noScale the normal 'scale' option squishes the data too much, so this defaults to TRUE
#' @param prior.plots print out prior plots? FALSE
#' @return a df of batch corrected data
#' @seealso \pkg{sva} \code{\link[sva]{ComBat}}
#' @examples
#' \dontrun{
#' df_new = hpgl_combatMod(df, batches, model)
#' }
#' @export
hpgl_combatMod <- function(dat, batch, mod, noScale=TRUE, prior.plots=FALSE) {
    par.prior <- TRUE
    numCovs <- NULL
    mod <- cbind(mod, batch)
    check <- apply(mod, 2, function(x) all(x == 1))
    mod <- as.matrix(mod[, !check])
    colnames(mod)[ncol(mod)] <- "Batch"
    if (sum(check) > 0 & !is.null(numCovs)) {
        numCovs <- numCovs - 1
    }
    ##    design <- sva:::design.mat(mod, numCov = numCovs)
    ## require.auto("survJamda")
    design <- survJamda::design.mat(mod)
    batches <- survJamda::list.batch(mod)
    n.batch <- length(batches)
    n.batches <- sapply(batches, length)
    n.array <- sum(n.batches)
    NAs <- any(is.na(dat))
    if (NAs) {
        message(paste0("Found ", sum(is.na(dat)), " missing data values."))
    }
    message("Standardizing data across genes\n")
    if (!NAs) {
        B.hat <- solve(t(design) %*% design) %*% t(design) %*% t(as.matrix(dat))
    } else {
        B.hat <- apply(dat, 1, Beta.NA, design)
    }
    grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]
    if (!NAs) {
        var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array, n.array)
    }
    else {
        var.pooled <- apply(dat - t(design %*% B.hat), 1, var, na.rm = T)
    }
    stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
    if (!is.null(design)) {
        tmp <- design
        tmp[, c(1:n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% B.hat)
    }
    s.data <- (dat - stand.mean) / (sqrt(var.pooled) %*% t(rep(1, n.array)))
    if (noScale) {
        m.data <- dat - stand.mean
        mse <- ((dat - t(design %*% B.hat))^2) %*% rep(1/(n.array - ncol(design)), n.array)
        hld <- NULL
        bayesdata <- dat
        for (k in 1:n.batch) {
            message(paste0("Fitting 'shrunk' batch ", k, " effects."))
            sel <- batches[[k]]
            gammaMLE <- rowMeans(m.data[, sel])
            mprior <- mean(gammaMLE, na.rm = TRUE)
            vprior <- var(gammaMLE, na.rm = TRUE)
            prop <- vprior / (mse / (length(sel)) + vprior)
            gammaPost <- prop * gammaMLE + (1 - prop) * mprior
            for (i in sel) {
                bayesdata[, i] <- bayesdata[, i] - gammaPost
            }
            stats <- data.frame(gammaPost=gammaPost, gammaMLE=gammaMLE, prop=prop)
            hld[[paste("Batch", k, sep=".")]] <- list(stats=stats, indices=sel, mprior=mprior, vprior=vprior)
        }
        message("Adjusting data for batch effects.")
        return(bayesdata)
    } else {
        message("Fitting L/S model and finding priors.")
        batch.design <- design[, 1:n.batch]
        if (!NAs) {
            gamma.hat <- solve(t(batch.design) %*% batch.design) %*% t(batch.design) %*% t(as.matrix(s.data))
        } else {
            gamma.hat <- apply(s.data, 1, Beta.NA, batch.design)
        }
        delta.hat <- NULL
        for (i in batches) {
            delta.hat <- rbind(delta.hat, apply(s.data[, i], 1, var, na.rm = T))
        }
        gamma.bar <- apply(gamma.hat, 1, mean)
        t2 <- apply(gamma.hat, 1, var)
        a.prior <- apply(delta.hat, 1, sva:::aprior)
        b.prior <- apply(delta.hat, 1, sva:::bprior)
        if (prior.plots & par.prior) {
            par(mfrow = c(2, 2))
            tmp <- density(gamma.hat[1, ])
            plot(tmp, type="l", main="Density Plot")
            xx <- seq(min(tmp$x), max(tmp$x), length = 100)
            lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
            stats::qqnorm(gamma.hat[1, ])
            stats::qqline(gamma.hat[1, ], col = 2)
            tmp <- stats::density(delta.hat[1, ])
            invgam <- 1 / stats::rgamma(ncol(delta.hat), a.prior[1], b.prior[1])
            tmp1 <- stats::density(invgam)
            plot(tmp, typ="l", main="Density Plot", ylim=c(0, max(tmp$y, tmp1$y)))
            lines(tmp1, col = 2)
            stats::qqplot(delta.hat[1, ], invgam, xlab="Sample Quantiles", ylab="Theoretical Quantiles")
            lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
            title("Q-Q Plot")
        }
        gamma.star <- delta.star <- NULL
        if (par.prior) {
            message("Finding parametric adjustments.")
            for (i in 1:n.batch) {
                temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                                     delta.hat[i, ], gamma.bar[i],
                                     t2[i], a.prior[i], b.prior[i])
                gamma.star <- rbind(gamma.star, temp[1, ])
                delta.star <- rbind(delta.star, temp[2, ])
            }
        } else {
            message("Finding nonparametric adjustments.")
            for (i in 1:n.batch) {
                temp <- sva:::int.prior(as.matrix(s.data[, batches[[i]]]), gamma.hat[i, ], delta.hat[i, ])
                gamma.star <- rbind(gamma.star, temp[1, ])
                delta.star <- rbind(delta.star, temp[2, ])
            }
        }
        message("Adjusting the Data.")
        bayesdata <- s.data
        j <- 1
        for (i in batches) {
            bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i,] %*% gamma.star)) /
                (sqrt(delta.star[j, ]) %*% t(rep(1, n.batches[j])))
            j <- j + 1
        }
        bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, n.array)))) + stand.mean
        return(bayesdata)
    }
}

#' Converts count matrix to log2 counts-per-million reads.
#'
#' Based on the method used by limma as described in the Law et al. (2014) voom
#' paper.
#'
#' @param counts read count matrix
#' @param lib.size  library size
#' @return log2-CPM read count matrix
#' @seealso \pkg{cbcbSEQ} \pkg{edgeR}
#' @examples
#' \dontrun{
#' l2cpm <- hpgl_log2cpm(counts)
#' }
#' @export
hpgl_log2cpm <- function(counts, lib.size=NULL) {
    if (is.null(lib.size)) {
        lib.size <- colSums(counts)
    }
    transposed_adjust <- t(counts + 0.5)
    cpm <- (transposed_adjust / (lib.size + 1)) * 1e+06
    l2cpm <- t(log2(cpm))
    return(l2cpm)
}

#' Normalize a dataframe/expt, express it, and/or transform it
#'
#' Sometime soon I am going to elipsis all these variables
#'
#' @param data some data
#' @param design   design dataframe must come with it
#' @param transform   defines whether to log(2|10) transform the
#' data. Defaults to raw.
#' @param norm   specify the normalization strategy.  Defaults to
#' raw.  This makes use of DESeq/EdgeR to provide: RLE, upperquartile,
#' size-factor, or tmm normalization.  I tend to like quantile, but there are
#' definitely corner-case scenarios for all strategies.
#' @param convert   defines the output type which may be raw, cpm,
#' rpkm, or cp_seq_m.  Defaults to raw.
#' @param batch   batch correction method to try out
#' @param batch1  column from design to get batch info
#' @param batch2   a second covariate to try
#' @param filter_low   choose whether to low-count filter the data.
#' @param annotations   is used for rpkm or sequence normalizations to
#' extract the lengths of sequences for normalization
#' @param entry_type   default gff entry to cull from
#' @param fasta   fasta genome for rpkm
#' @param verbose  talk
#' @param thresh   threshold for low count filtering
#' @param min_samples   minimum samples for low count filtering
#' @param noscale   used by combatmod
#' @param p   for povera genefilter
#' @param A   for povera genefilter
#' @param k   for kovera genefilter
#' @param cv_min   for genefilter cv
#' @param cv_max   for genefilter cv
#' @param ... I should put all those other options here
#' @return edgeR's DGEList expression of a count table.  This seems to
#' me to be the easiest to deal with.
#' @seealso \link[edgeR]{cpm} \link[edgeR]{rpkm}
#' \link{hpgl_rpkm} \link[cbcbSEQ]{filterCounts} \link[DESeq2]{DESeqDataSetFromMatrix}
#' \link[DESeq]{estimateSizeFactors} \link[edgeR]{DGEList} \link[edgeR]{calcNormFactors}
#' @export
#' @examples
#' \dontrun{
#' df_raw = hpgl_norm(expt=expt)  ## Only performs low-count filtering
#' df_raw = hpgl_norm(df=a_df, design=a_design) ## Same, but using a df
#' df_ql2rpkm = hpgl_norm(expt=expt, norm='quant', transform='log2',
#'                        convert='rpkm')  ## Quantile, log2, rpkm
#' count_table = df_ql2rpkm$counts
#' }
hpgl_norm <- function(data, design=NULL, transform="raw", norm="raw",
                      convert="raw", batch="raw", batch1="batch", batch2=NULL,
                      filter_low=FALSE, annotations=NULL, entry_type="gene",
                      fasta=NULL, verbose=FALSE, thresh=2, min_samples=2,
                      noscale=TRUE, p=0.01, A=1, k=1, cv_min=0.01, cv_max=1000, ...) {
    lowfilter_performed <- FALSE
    norm_performed <- "raw"
    convert_performed <- "raw"
    transform_performed <- "raw"
    batch_performed <- "raw"
    data_class <- class(data)[1]
    if (data_class == 'expt') {
        design <- data$design
        data <- Biobase::exprs(data$expressionset)
    } else if (data_class == 'ExpressionSet') {
        data <- Biobase::exprs(data)
    } else if (data_class == 'list') {
        data <- data$count_table
        if (is.null(data)) {
            stop("The list provided contains no count_table.")
        }
    } else if (data_class == 'matrix' | data_class == 'data.frame') {
        data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }
    count_table <- as.matrix(data)
    expt_design <- design
    raw_libsize <- colSums(count_table)
    original_counts <- list(libsize=raw_libsize, counts=count_table)

    if (verbose) {
        message("This function performs normalization in a static order: low-count filter, normalization, batch, conversion, transform")
        message("These steps may be mixed/matched with the following functions: lowfilter_counts, normalize_counts, batch_counts, convert_counts, transform_counts")
    }

    ## Step 1: Low-count filtering
    lowfiltered_counts <- NULL
    if (filter_low != FALSE) {
        if (verbose) {
            message(paste0("Performing low-count filter with: ", filter_low))
        }
        count_table <- lowfilter_counts(count_table, type=filter_low, p=p, A=A, k=k, cv_min=cv_min, cv_max=cv_max, thresh=2, min_samples=2)
        ##count_table = lowfilter_counts(count_table, type=filter_low)
        lowfilter_performed <- filter_low
    }

    ## Step 2: Normalization
    ## This section handles the various normalization strategies
    ## If nothing is chosen, then the filtering is considered sufficient
    normalized_counts <- NULL
    if (norm != "raw") {
        if (verbose) {
            message(paste0("Applying normalization: ", norm))
        }
        if (is.null(expt_design)) {
            message("The experimental design is null.  Some normalizations will therefore fail.")
            message("If you receive an error about an object with no dimensions, that is likely why.")
        }
        normalized_counts <- normalize_counts(count_table, expt_design, norm=norm)
        count_table <- normalized_counts$count_table
        norm_performed <- norm
    }

    ## Step 3: Convert the data to (likely) cpm
    ## The following stanza handles the three possible output types
    ## cpm and rpkm are both from edgeR
    ## They have nice ways of handling the log2 which I should consider
    converted_counts <- NULL
    if (convert != "raw") {
        if (verbose) {
            message(paste0("Setting output type as: ", convert))
        }
        converted_counts <- convert_counts(count_table, convert=convert, annotations=annotations, fasta=fasta, entry_type=entry_type)
        count_table <- converted_counts$count_table
        convert_performed <- convert
    }

    ## Step 4: Batch correction
    batched_counts <- NULL
    if (batch != "raw") {
        if (verbose) {
            message(paste0("Applying: ", batch, " batch correction(raw means nothing)."))
        }
        ## batched_counts = batch_counts(count_table, batch=batch, batch1=batch1, batch2=batch2, design=design, ...)
        tmp_counts <- try(batch_counts(count_table, batch=batch, batch1=batch1, batch2=batch2, design=expt_design), silent=TRUE)
        batched_counts <- list(count_table=count_table)
        if (class(tmp_counts) == 'try-error') {
            warning("The batch_counts called failed.  Returning non-batch reduced data.")
        } else {
            batched_counts <- tmp_counts
        }
        count_table <- batched_counts$count_table
        batch_performed <- batch
    }

    ## Step 5: Transformation
    ## Finally, this considers whether to log2 the data or no
    transformed_counts <- NULL
    if (transform != "raw") {
        if (verbose) {
            message(paste0("Applying: ", transform, " transformation."))
        }
        transformed_counts <- transform_counts(count_table, transform=transform, converted=convert_performed, ...)
        count_table <- transformed_counts$count_table
        transform_performed <- transform
    }

    final_counts <- list(count_table=count_table, libsize=colSums(count_table))
    ret_list <- list(
        lowfilter_performed=lowfilter_performed, norm_performed=norm_performed, convert_performed=convert_performed,
        transform_performed=transform_performed, batch_performed=batch_performed,
        original_counts=original_counts, lowfiltered_counts=lowfiltered_counts,
        normalized_counts=normalized_counts, converted_counts=converted_counts,
        transformed_counts=transformed_counts, batched_counts=batched_counts,
        final_counts=final_counts, count_table=final_counts$count_table,
        libsize=final_counts$libsize
    )
    return(ret_list)
}

#' A hacked copy of Kwame's qsmooth/qstats code
#'
#' I made a couple small changes to Kwame's qstats() function to make
#' it not fail when on corner-cases.  I sent him a diff, but haven't
#' checked to see if it was useful yet.
#'
#' @param data count table to modify
#' @param groups factor of the experimental conditions
#' @param refType method for grouping conditions
#' @param groupLoc method for grouping groups
#' @param window a window, for looking!
#' @param verbose talky talky
#' @param groupCol column to define conditions
#' @param plot plot the quantiles?
#' @param ... more options
#' @return data a new data frame of normalized counts
#' @seealso \pkg{qsmooth}
#' @examples
#' \dontrun{
#' df <- hpgl_qshrink(data)
#' }
#' @export
hpgl_qshrink <- function(data=NULL, groups=NULL, refType="mean",
                        groupLoc="mean", window=99, verbose=FALSE,
                        groupCol=NULL, plot=TRUE, ...) {
    data <- as.matrix(data)
    if (is.null(groups)) {
        message("Groups were not provided.  Performing a simple quantile")
        message("normalization. This is probably not what you actually want!")
        count_rownames <- rownames(data)
        count_colnames <- colnames(data)
        normExprs <- preprocessCore::normalize.quantiles(as.matrix(data), copy=TRUE)
        rownames(normExprs) <- count_rownames
        colnames(normExprs) <- count_colnames
        return(normExprs)
    }
    res <- hpgl_qstats(Biobase::exprs, groups, refType=refType,
                       groupLoc=groupLoc, window=window)
    QBETAS <- res$QBETAS
    Qref <- res$Qref
    X <- res$model
    w <- res$smoothWeights
    wQBETAS <- QBETAS * (1 - w)
    wQBETAS <- X %*% t(wQBETAS)
    wQref <- Qref * w
    wQref <- matrix(rep(1, nrow(X)), ncol=1) %*% t(wQref)
    normExprs <- t(wQBETAS + wQref)
    RANKS <- t(matrixStats::colRanks(data, ties.method="average"))
    for (k in 1:ncol(normExprs)) {
        x <- normExprs[, k]
        normExprs[, k] <- x[RANKS[, k]]
    }

    aveTies = function (ranks, y) {
        tab = table(ranks)
        sel = tab > 1
        if (sum(sel) != 0) {
            ties = as.numeric(names(tab[sel]))
            for (k in ties) {
                sel = ranks==k
                y[sel] = mean(y[sel])
            }
        }
        y
    }
    normExprs <- aveTies(RANKS, normExprs)

    rownames(normExprs) <- rownames(data)
    colnames(normExprs) <- colnames(data)
    if (plot) {
        oldpar <- par(mar=c(4, 4, 1.5, 0.5))
        lq <- length(Qref)
        u <- (1:lq - 0.5)/lq
        if (length(u) > 10000) {
            sel <- sample(1:lq, 10000)
            plot(u[sel], w[sel], pch=".", main="Quantile reference weights",
                 xlab="u (norm. gene ranks)", ylab="Weight", ylim=c(0, 1), ...)
            ## plot(u[sel], w[sel], pch=".", main="Quantile reference weights",
            ##      xlab="u (norm. gene ranks)", ylab="Weight", ylim=c(0, 1))
        } else {
            plot(u, w, pch=".", main="Quantile reference weights",
                 xlab="u (norm. gene ranks)", ylab="Weight", ylim=c(0, 1), ...)
            ## plot(u, w, pch=".", main="Quantile reference weights",
            ##      xlab="u (norm. gene ranks)", ylab="Weight", ylim=c(0, 1))
        }
        abline(h=0.5, v=0.5, col="red", lty=2)
        par(oldpar)
    }
    return(normExprs)
}

#' A caller for different low-count filters
#'
#' @param count_table  some counts to filter
#' @param type   Filtering method to apply (cbcb, pofa, kofa, cv right now)
#' @param p  For pofa()
#' @param A   For pofa()
#' @param k   For kofa()
#' @param cv_min   For cv()
#' @param cv_max   For cv()
#' @param thresh   Minimum threshold across samples for cbcb
#' @param min_samples   Minimum number of samples for cbcb
#' @return a data frame of lowfiltered counts
#' @seealso \pkg{genefilter}
#' @examples
#' \dontrun{
#' new <- lowfilter_counts(old)
#' }
#' @export
lowfilter_counts <- function(count_table, type='cbcb', p=0.01, A=1, k=1,
                             cv_min=0.01, cv_max=1000, thresh=2, min_samples=2) {
    if (tolower(type) == 'povera') {
        type <- 'pofa'
    } else if (tolower(type) == 'kovera') {
        type <- 'kofa'
    }
    lowfiltered_counts <- NULL
    if (type == 'cbcb') {
        lowfiltered_counts <- cbcb_filter_counts(count_table, threshold=thresh,
                                                 min_samples=min_samples)
    } else if (type == 'pofa') {
        lowfiltered_counts <- genefilter_pofa_counts(count_table, p=p, A=A)
    } else if (type == 'kofa') {
        lowfiltered_counts <- genefilter_kofa_counts(count_table, k=k, A=A)
    } else if (type == 'cv') {
        lowfiltered_counts <- genefilter_cv_counts(count_table, cv_min=cv_min,
                                                   cv_max=cv_max)
    } else {
        message("Did not recognize the filtering argument, defaulting to cbcb's.
 Recognized filters are: 'cv', 'kofa', 'pofa', 'cbcb'
")
        lowfiltered_counts <- cbcb_filter_counts(count_table, threshold=thresh,
                                                 min_samples=min_samples)
    }
    count_table <- lowfiltered_counts$count_table
    return(count_table)
}

#' A hacked copy of Kwame's qsmooth/qstats code
#'
#' I made a couple small changes to Kwame's qstats() function to make
#' it not fail when on corner-cases.  I sent him a diff, but haven't
#' checked to see if it was useful yet.
#'
#' @param data the initial count data
#' @param groups the experimental conditions as a factor
#' @param refType  (or median) the method to separate groups
#' @param groupLoc   I don't remember
#' @param window window for basking
#' @return new data
#' @examples
#' \dontrun{
#' qstatted <- hpgl_qstats(data, conditions)
#' }
#' @export
hpgl_qstats <- function (data, groups, refType="mean",
                         groupLoc="mean", window=99) {
    ## require.auto("matrixStats")
    Q <- apply(data, 2, sort)
    if (refType == "median") {
        Qref <- matrixStats::rowMedians(Q)
    }
    if (refType == "mean") {
        Qref <- rowMeans(Q)
    }

    QBETAS <- c()
    SIGMA <- c()
    uGroups <- unique(groups)
    for (g in uGroups) {
        index <- (g == groups)
        if (sum(index) == 1) {
            message(paste0("There was only replicate of type: ", g))
            message("This will likely do terrible things to qsmooth.")
            QBETAS <- cbind(QBETAS, Q[, index])
            SIGMA <- cbind(SIGMA, 0)
        } else if (sum(index) > 1) {
            if (groupLoc == "mean") {
                QBETAS <- cbind(QBETAS, rowMeans(Q[, index]))
                SIGMA <- cbind(SIGMA, matrixStats::rowVars(Q[, g == groups]))
            } else if (groupLoc == "median") {
                QBETAS <- cbind(QBETAS, matrixStats::rowMedians(Q[, index]))
                SIGMA <- cbind(SIGMA, (matrixStats::rowMads(Q[, g == groups]))^2)
            }
        } else {
            warning(paste0("There were 0 of type: ", g))
        }
    }
    colnames(QBETAS) <- uGroups
    colnames(SIGMA) <- uGroups
    if (groupLoc == "mean") {
        TAU <- matrixStats::rowVars(QBETAS)
        SIGMA <- rowMeans(SIGMA)
    } else { ## median
        TAU <- matrixStats::rowMads(QBETAS)^2
        SIGMA <- matrixStats::rowMedians(SIGMA)
    }
    roughWeights <- SIGMA/(SIGMA + TAU)
    roughWeights[is.nan(roughWeights)] = 0 ## is this backward?
    roughWeights[SIGMA < 10^(-6) & TAU < 10^(-6)] = 1
    smoothWeights <- stats::runmed(roughWeights, k=window, endrule="constant")
    qstats_model <- model.matrix(~0 + factor(groups, levels=uGroups))
    qstats_result <- list(Q=Q, Qref=Qref, QBETAS=QBETAS, TAU=TAU,
                          SIGMA=SIGMA, roughWeights=roughWeights,
                          smoothWeights=smoothWeights, model=qstats_model)
    return(qstats_result)
}

#' Reads/(kilobase(gene) * million reads)
#'
#' Express a data frame of counts as reads per kilobase(gene) per
#' million(library).
#'
#' This function wraps EdgeR's rpkm in an attempt to make sure that
#' the required gene lengths get sent along.
#'
#' @param df a data frame of counts, alternately an edgeR DGEList
#' @param annotations containing gene lengths, defaulting to
#' 'gene_annotations'
#' @return rpkm_df a data frame of counts expressed as rpkm
#' @seealso \pkg{edgeR} and \code{\link[edgeR]{cpm}} \code{\link[edgeR]{rpkm}}
#' @examples
#' \dontrun{
#' rpkm_df = hpgl_rpkm(df, annotations=gene_annotations)
#' }
#' @export
hpgl_rpkm <- function(df, annotations=get0('gene_annotations')) {
    if (class(df) == "edgeR") {
        df <- df$counts
    }
    df_in <- as.data.frame(df[rownames(df) %in% rownames(annotations), ])
    if (dim(df_in)[1] == 0) {
        message("When the annotations and df were checked against each other
  the result was null.  Perhaps your annotation or df's rownames are not set?
  Going to attempt to use the column 'ID'.
")
        rownames(annotations) = make.names(annotations$ID, unique=TRUE)
        df_in <- as.data.frame(df[rownames(df) %in% rownames(annotations), ])
        if (dim(df_in)[1] == 0) {
            stop("The ID column failed too.")
        }
    }
    colnames(df_in) <- colnames(df)
    merged_annotations <- merge(df, annotations, by="row.names")
    rownames(merged_annotations) <- merged_annotations[,"Row.names"]
    ##rownames(df_in) = merged_annotations[,"Row.names"]
    ## Sometimes I am stupid and call it length...
    lenvec <- NULL
    if (is.null(merged_annotations$width)) {
        lenvec <- as.vector(merged_annotations$length)
    } else {
        lenvec <- as.vector(merged_annotations$width)
    }
    names(lenvec) <- rownames(merged_annotations)
    requireNamespace("edgeR")
    rpkm_df <- edgeR::rpkm(df_in, gene.length=lenvec)
    colnames(rpkm_df) <- colnames(df)
    return(rpkm_df)
}

#' Filter low-count genes from a data set using cbcbSEQ::filterCounts()
#'
#' @param count_table  input data frame of counts by sample
#' @param thresh   lower threshold of counts (default: 4)
#' @param min_samples   minimum number of samples (default: 2)
#' @param verbose   If set to true, prints number of genes removed / remaining
#' @return dataframe of counts without the low-count genes
#' @seealso \link[cbcbSEQ]{log2CPM} which this uses to decide what to keep
#' @examples
#' \dontrun{
#'  filtered_table = cbcb_lowfilter_counts(count_table)
#' }
#' @export
cbcb_lowfilter_counts <- function(count_table, thresh=2,
                                  min_samples=2, verbose=FALSE) {
    original_dim <- dim(count_table)
    count_table <- as.matrix(cbcbSEQ::filterCounts(count_table, thresh=thresh,
                                                   min_samples=min_samples))
    if (verbose) {
        following_dim <- dim(count_table)
        lost_rows <- original_dim[1] - following_dim[1]
        message(paste0("Low count filtering cost: ", lost_rows, " gene(s)."))
    }
    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

#'   Perform a simple normalization of a count table
#'
#' @param data A matrix of count data
#' @param design  A dataframe describing the experimental design
#' (conditions/batches/etc)
#' @param norm  A normalization to perform:
#' 'sf|quant|qsmooth|tmm|upperquartile|tmm|rle'
#' I keep wishy-washing on whether design is a required argument.
#' @return dataframe of normalized(counts)
#' @seealso \pkg{edgeR} \pkg{limma} \pkg{DESeq2}
#' @examples
#' \dontrun{
#' norm_table = normalize_counts(count_table, design=design, norm='qsmooth')
#' }
#' @export
normalize_counts <- function(data, design=NULL, norm="raw") {
    ## Note that checkUsage flagged my 'libsize = ' calls
    ## I set norm_libsize at the bottom of the function
    ## but perhaps instead I should be using these libsizes?
    data_class <- class(data)[1]
    if (data_class == 'expt') {
        design <- data$design
        count_table <- Biobase::exprs(data$expressionset)
    } else if (data_class == 'ExpressionSet') {
        count_table <- Biobase::exprs(data)
    } else if (data_class == 'list') {
        count_table <- data$count_table
        if (is.null(data)) {
            stop("The list provided contains no count_table.")
        }
    } else if (data_class == 'matrix' | data_class == 'data.frame') {
        ## some functions prefer matrix, so I am keeping this explicit
        count_table <- as.data.frame(data)
    } else {
        stop(paste0("You provided a class of type: ", data_class, ".
This works with: expt, ExpressionSet, data.frame, and matrices.
"))
    }
    if (norm == "sf") {
        ## Size-factored normalization is a part of DESeq
        factors <- DESeq2::estimateSizeFactorsForMatrix(count_table)
        num_rows <- dim(count_table)[1]
        sf_counts <- count_table / do.call(rbind, rep(list(factors), num_rows))
        ##sf_counts = counts / (libsizes * factors)
        count_table <- as.matrix(sf_counts)
        norm_performed <- 'sf'
    } else if (norm == 'sf2') {
        original_cols <- colnames(count_table)
        conds <- design$conditions
        if (is.null(conds)) {
            conds <- original_cols
        }
        cds <- DESeq::newCountDataSet(count_table, conditions=conds)
        factors <- BiocGenerics::estimateSizeFactors(cds)
        count_table <- BiocGenerics::counts(factors, normalized=TRUE)
        norm_performed <- 'sf2'
    } else if (norm == 'vsd') {
        original_cols <- colnames(count_table)
        conds <- design$conditions
        if (is.null(conds)) {
            conds <- original_cols
        }
        cds <- DESeq::newCountDataSet(count_table, conditions=conds)
        factors <- BiocGenerics::estimateSizeFactors(cds)
        dispersions <- BiocGenerics::estimateDispersions(factors, method='blind')
        count_table <- DESeq::getVarianceStabilizedData(dispersions)
        norm_performed <- 'vsd'
    } else if (norm == "quant") {
        # Quantile normalization (Bolstad et al., 2003)
        count_rownames <- rownames(count_table)
        count_colnames <- colnames(count_table)
        count_table <- preprocessCore::normalize.quantiles(as.matrix(count_table), copy=TRUE)
        rownames(count_table) <- count_rownames
        colnames(count_table) <- count_colnames
        norm_performed <- 'quant'
    } else if (norm == "qsmooth") {
        count_table <- qsmooth::qsmooth(count_table, groups=design$condition, plot=TRUE)
        norm_performed <- "qsmooth"
    } else if (norm == "qshrink") {
        count_table <- hpgl_qshrink(exprs=count_table, groups=design$condition,
                                    verbose=TRUE, plot=TRUE)
        norm_performed <- 'qshrink'
    } else if (norm == "qshrink_median") {
        count_table <- hpgl_qshrink(exprs=count_table, groups=design$condition,
                                    verbose=TRUE, plot=TRUE, refType="median",
                                    groupLoc="median", window=50)
        norm_performed <- 'qshrink_median'
    } else if (norm == "tmm") {
        ## TMM normalization is documented in edgeR
        ## Set up the edgeR data structure
        count_table <- edgeR::DGEList(counts=count_table)
        norms <- edgeR::calcNormFactors(count_table, method="TMM")
        ## libsizes = count_table$samples$lib.size
        factors <- norms$samples$norm.factors
        counts <- norms$counts
        tmm_counts <- counts / factors
        count_table <- as.matrix(tmm_counts)
        norm_performed <- "tmm"
    } else if (norm == "upperquartile") {
        ## Get the tmm normalization factors
        count_table <- edgeR::DGEList(counts=count_table)
        norms <- edgeR::calcNormFactors(count_table, method="upperquartile")
        ## libsizes = count_table$samples$lib.size
        factors <- norms$samples$norm.factors
        counts <- norms$counts
        tmm_counts <- counts / factors
        count_table <- as.matrix(tmm_counts)
        norm_performed <- "upperquartile"
    } else if (norm == "rle") {
        ## Get the tmm normalization factors
        count_table <- edgeR::DGEList(counts=count_table)
        norms <- edgeR::calcNormFactors(count_table, method="RLE")
        ## libsizes = count_table$samples$lib.size
        factors <- norms$samples$norm.factors
        counts <- norms$counts
        tmm_counts <- counts / factors
        count_table <- as.matrix(tmm_counts)
        norm_performed <- "rle"
    } else {
        message("Did not recognize the normalization, leaving the table alone.
  Recognized normalizations include: 'qsmooth', 'sf', 'sf2', 'vsd', 'quant',
  'tmm', 'qsmooth_median', 'upperquartile', and 'rle.'
")
        count_table <- as.matrix(count_table)
    }
    norm_libsize <- colSums(count_table)
    norm_counts <- list(count_table=count_table, libsize=norm_libsize, norm_performed=norm_performed)
    return(norm_counts)
}

#'   Replace the data of an expt with normalized data.
#'
#' @param expt   The original expt
#' @param transform   The transformation desired (raw, log2, log, log10)
#' @param norm   How to normalize the data (raw, quant, sf, upperquartile, tmm, rle)
#' @param convert   Conversion to perform (raw, cpm, rpkm, cp_seq_m)
#' @param batch   Batch effect removal tool to use (limma sva fsva ruv etc)
#' @param filter_low   Filter out low sequences (cbcb, pofa, kofa, others?)
#' @param annotations  used for rpkm, a df
#' @param fasta  fasta file for cp_seq_m counting of oligos
#' @param entry_type   for getting genelengths by feature type (rpkm or cp_seq_m)
#' @param verbose   talk?
#' @param use_original   whether to use the backup data in the expt class
#' @param batch1   experimental factor to extract first
#' @param batch2   a second factor to remove (only with limma's removebatcheffect())
#' @param thresh   for cbcb_lowfilter
#' @param min_samples   for cbcb_lowfilter
#' @param p   for genefilter's pofa
#' @param A   for genefilter's pofa
#' @param k   for genefilter's kofa
#' @param cv_min   for genefilter's cv()
#' @param cv_max  for genefilter's cv()
#' @param ... more options
#' @return a new expt object with normalized data and the original data saved as 'original_expressionset'
#' @seealso \pkg{genefilter} \pkg{cbcbSEQ} \pkg{limma} \pkg{sva} \pkg{edgeR} \pkg{DESeq2}
#' @examples
#' \dontrun{
#' normed <- normalize_expt(exp, transform='log2', norm='rle', convert='cpm',
#'                          batch='raw', filter_low='pofa')
#' normed_batch <- normalize_expt(exp, transform='log2', norm='rle', convert='cpm',
#'                                batch='sva', filter_low='pofa')
#' }
#' @export
normalize_expt <- function(expt, ## The expt class passed to the normalizer
    ## choose the normalization strategy
    transform="raw", norm="raw", convert="raw", batch="raw", filter_low=FALSE,
    ## annotations used for rpkm/cpseqm, original may be used to ensure double-normalization isn't performed.
    annotations=NULL, fasta=NULL, entry_type="gene", verbose=FALSE, use_original=FALSE,
    batch1="batch", batch2=NULL, ## extra parameters for batch correction
    thresh=2, min_samples=2, p=0.01, A=1, k=1, cv_min=0.01, cv_max=1000,  ## extra parameters for low-count filtering
    ...) {
    new_expt <- expt
    current <- expt$expressionset
    if (is.null(new_expt$original_expressionset)) {
        new_expt$original_expressionset = new_expt$expressionset
    } else {
        message("This function will replace the expt$expressionset slot with:")
        type <- ""
        if (transform != "raw") {
            type <- paste0(type, transform, '(')
        }
        if (norm != "raw") {
            type <- paste0(type, norm, '(')
        }
        if (convert != "raw") {
            type <- paste0(type, convert, '(')
        }
        if (filter_low != FALSE) {
            type <- paste0(type, 'low-filter(')
        }
        if (batch != "raw") {
            type <- paste0(type, 'batch-correct(')
        }
        type <- paste0(type, 'data')
        if (transform != 'raw') {
            type <- paste0(type, ')')
        }
        if (norm != "raw") {
            type <- paste0(type, ')')
        }
        if (convert != "raw") {
            type <- paste0(type, ')')
        }
        if (filter_low != FALSE) {
            type <- paste0(type, ')')
        }
        if (batch != "raw") {
            type <- paste0(type, ')')
        }
        message(type)
        message("It saves the current data into a slot named:
 expt$backup_expressionset. It will also save copies of each step along the way
 in expt$normalized with the corresponding libsizes. Keep the libsizes in mind
 when invoking limma.  The appropriate libsize is the non-log(cpm(normalized)).
 This is most likely kept in the slot called:
 'new_expt$normalized$normalized_counts$libsize' which is copied into
 new_expt$best_libsize
")
    }
    if (filter_low == FALSE) {
        message("Filter low is false, this should likely be set to something, good
 choices include cbcb, kofa, pofa (anything but FALSE).  If you want this to
 stay FALSE, keep in mind that if other normalizations are performed, then the
 resulting libsizes are likely to be strange (potentially negative!)
")
    }
    if (transform == "raw") {
        message("Leaving the data in its current base format, keep in mind that
 some metrics are easier to see when the data is log2 transformed, but
 EdgeR/DESeq don't like transformed data.
")
    }
    if (convert == "raw") {
        message("Leaving the data unconverted.  It is often advisable to cpm/rpkm
 the data to normalize for sampling differences, keep in mind though that rpkm
 has some annoying biases, and voom() by default does a cpm (though hpgl_voom()
 will try to detect this).
")
    }
    if (norm == "raw") {
        message("Leaving the data unnormalized.  This is necessary for DESeq, but
 EdgeR/limma might benefit from normalization.  Good choices include quantile,
 size-factor, tmm, etc.
")
    }
    if (batch == "raw") {
        message("Not correcting the count-data for batch effects.  If batch is
 included in EdgerR/limma's model, then this is probably wise; but in extreme
 batch effects this is a good parameter to play with.
")
    }
    new_expt$backup_expressionset <- new_expt$expressionset
    old_data <- Biobase::exprs(expt$original_expressionset)
    design <- expt$design
    normalized <- hpgl_norm(old_data, design=design, transform=transform,
                            norm=norm, convert=convert, batch=batch,
                            batch1=batch1, batch2=batch2,
                            filter_low=filter_low, annotations=annotations,
                            fasta=fasta, verbose=verbose, thresh=thresh,
                            min_samples=min_samples, p=p, A=A, k=k,
                            cv_min=cv_min, cv_max=cv_max, entry_type=entry_type)
    final_normalized <- normalized$final_counts
    libsizes <- final_normalized$libsize
    normalized_data <- as.matrix(final_normalized$count_table)
    Biobase::exprs(current) <- normalized_data
    new_expt$normalized <- normalized
    new_expt$norm_libsize <- libsizes
    new_expt$expressionset <- current
    new_expt$filtered <- filter_low
    new_expt$transform <- transform
    new_expt$best_libsize <- new_expt$normalized$normalized_counts$libsize
    new_expt$norm <- norm
    new_expt$convert <- convert
    new_expt$batch <- batch
    return(new_expt)
}

#' Perform a simple transformation of a count table (log2)
#'
#' the add argument is only important if the data was previously cpm'd because that does a +1, thus
#' this will avoid a double+1 on the data.
#'
#' @param count_table  A matrix of count data
#' @param transform   A type of transformation to perform: log2/log10/log
#' @param converted   Whether or not the data has been converted.
#' @param base   for other log scales
#' @param add   to avoid attempting a log(0)
#' @return dataframe of logx(counts)
#' @examples
#' \dontrun{
#' filtered_table = transform_counts(count_table, transform='log2', converted='cpm')
#' }
#' @export
transform_counts <- function(count_table, transform="raw", converted="raw",
                             base=NULL, add=0.5) {
    ## if (converted != "cpm") {
    ##     count_table = count_table + 1
    ## }
    num_zero <- sum(count_table == 0)
    num_low <- sum(count_table < 0)
    if (num_low > 0) {
        message(paste0("transform_counts: Found ", num_low, " values less than 0."))
    }
    if (num_zero > 0) {
        message(paste0("transform_counts: Found ", num_zero, " values equal to 0, adding ", add, "
 to the matrix.
"))
        count_table <- count_table + add
    }
    if (!is.null(base)) {
        count_table <- (log(count_table) / log(base))
    } else if (transform == "log2") {
        count_table <- log2(count_table)
    } else if (transform == "log10") {
        count_table <- log10(count_table)
    } else if (transform == "log") {  ## Natural log
        count_table <- log(count_table)  ## Apparently log1p does this.
    } else {
        message("Did not recognize the transformation, leaving the table.
 Recognized transformations include: 'log2', 'log10', 'log'
")
    }
    libsize <- colSums(count_table)
    counts <- list(count_table=count_table, libsize=libsize)
    return(counts)
}

## EOF
