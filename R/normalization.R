#' Express a data frame of counts as reads per kilobase(gene) per
#' million(library).
#'
#' @param df a data frame of counts, alternately an edgeR DGEList
#' @param annotations containing gene lengths, defaulting to
#' 'gene_annotations'
#' 
#' @return rpkm_df a data frame of counts expressed as rpkm
#' @seealso \code{\link{edgeR}} and \code{\link{cpm}}
#' @export
#' @examples
#' ## rpkm_df = hpgl_rpkm(df, annotations=gene_annotations)
hpgl_rpkm = function(df, annotations=gene_annotations) {
    if (class(df) == "edgeR") {
        df = df$counts
    }
    df = df[rownames(df) %in% rownames(annotations),]
    merged_annotations = merge(df, annotations, by="row.names")
    gene_lengths = merged_annotations$width
    rpkm_df = edgeR::rpkm(df, gene.length=gene_lengths)
    return(rpkm_df)
}

#' Converts count matrix to log2 counts-per-million reads.
#'
#' Based on the method used by limma as described in the Law et al. (2014) voom
#' paper.
#'
#' @param counts read count matrix
#' 
#' @return log2-CPM read count matrix
#' @export
#'
hpgl_log2cpm = function(counts, lib.size=NULL) {
    if (is.null(lib.size)) {
        lib.size = colSums(counts)
    }
    transposed_adjust = t(counts + 0.5)
    cpm = (transposed_adjust / (lib.size + 1)) * 1e+06
    l2cpm = t(log2(cpm))
    return(l2cpm)
}

#' Express a data frame of counts as reads per pattern per
#' million(library).
#'
#' @param counts read count matrix
#' @param pattern pattern to search against.  Defaults to 'TA'
#' @param fasta a fasta genome to search
#' @param gff the gff set of annotations to define start/ends of genes.
#' @param entry_type which type of gff entry to search against.  Defaults to 'gene'.
#' 
#' @return The 'RPseqM' counts
#' @export
divide_seq = function(counts, pattern="TA", fasta="testme.fasta", gff="testme.gff", entry_type="gene") {
    raw_seq = FaFile(fasta)
    gff_entries = import.gff3(gff, asRangedData=FALSE)
    cds_entries = subset(gff_entries, type==entry_type)
    names(cds_entries) = cds_entries$locus_tag
    cds_seq = getSeq(raw_seq, cds_entries)
    names(cds_seq) = cds_entries$locus_tag
    dict = PDict(pattern, max.mismatch=0)
    result = vcountPDict(dict, cds_seq)
    num_tas = data.frame(name=names(cds_seq), tas=as.data.frame(t(result)))
    colnames(num_tas) = c("name","TAs")
    num_tas$TAs = num_tas$TAs + 1
    factor = median(num_tas$TAs)
    num_tas$TAs = num_tas$TAs / factor
    merged_tas = merge(counts, num_tas, by.x="row.names", by.y="name")
    rownames(merged_tas) = merged_tas$Row.names
    merged_tas = merged_tas[-1]
    merged_tas = merged_tas / merged_tas$TAs
    merged_tas = merged_tas[, !(colnames(merged_tas) %in% c("TAs"))]
    return(merged_tas)
}

#' Filter low-count genes from a data set.
#'
#' @param df input data frame of counts by sample
#' @param threshold lower threshold of counts (default: 4)
#' @param min_samples minimum number of samples (default: 2)
#' @param verbose If set to true, prints number of genes removed / remaining
#' @return dataframe of counts without the low-count genes
#' @seealso \code{\link{log2CPM}} which this uses to decide what to keep
#' @export
#' @examples
#' ## filtered_table = filter_counts(count_table)
cbcb_filter_counts = function(count_table, threshold=2, min_samples=2, verbose=FALSE) {
    ## I think having a log2cpm here is kind of weird, because the next step in processing is to cpm the data.
    ##cpms = 2^log2CPM(counts, lib.size=lib.size)$y
    ## cpms = 2^hpgl_log2cpm(counts)
    num_before = nrow(count_table)

    if (class(count_table) == 'ExpressionSet') {
        keep = rowSums(exprs(count_table) > threshold) >= min_samples
    } else {
        keep = rowSums(count_table > threshold) >= min_samples
    }

    count_table = count_table[keep,]

    if (verbose) {
        print(sprintf("Removing %d low-count genes (%d remaining).",
                      num_before - nrow(count_table), nrow(count_table)))
    }

    libsize = colSums(count_table)
    counts = list(count_table=count_table, libsize=libsize)    
    return(counts)
}

#' Filter low-count genes from a data set using filterCounts()
#'
#' @param df input data frame of counts by sample
#' @param threshold lower threshold of counts (default: 4)
#' @param min_samples minimum number of samples (default: 2)
#' @param verbose If set to true, prints number of genes removed / remaining
#' @return dataframe of counts without the low-count genes
#' @seealso \code{\link{log2CPM}} which this uses to decide what to keep
#' @export
#' @examples
#' ## filtered_table = lowfilter_counts(count_table)
lowfilter_counts = function(count_table, thresh=2, min_samples=2, verbose=FALSE) {
    original_dim = dim(count_table)
    count_table = as.matrix(filterCounts(count_table, thresh=thresh, min_samples=min_samples))
    if (verbose) {
        following_dim = dim(count_table)
        lost_rows = original_dim[1] - following_dim[1]
        print(paste("Low count filtering cost:", lost_rows, "gene(s)."))
    }
    libsize = colSums(count_table)
    counts = list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Filter low-count genes from a data set using genefilter's pOverA()
#'
#' I keep thinking this function is pofa... oh well.
#'
#' @param counts input data frame of counts by sample
#' @param p a minimum proportion of each gene's counts/sample to be greater than a minimum(A) (defaults to 0.01)
#' @param A the minimum number of counts in the above proportion
#' @param verbose If set to true, prints number of genes removed / remaining
#' @return dataframe of counts without the low-count genes
#' @seealso \code{\link{genefilter}} \code{\link{pOverA}} which this uses to decide what to keep
#' @export
#' @examples
#' ## filtered_table = genefilter_pofa_counts(count_table)
genefilter_pofa_counts = function(count_table, p=0.01, A=100, verbose=FALSE) {
    ## genefilter has functions to work with expressionsets directly, but I think I will work merely with tables in this.
    num_before = nrow(count_table)

    if (class(count_table) == 'ExpressionSet') {
        counts = exprs(count_table)
    }
    test = pOverA(p=p, A=A)
    filter_list = filterfun(test)
    answer = genefilter(count_table, filter_list)
    count_table = count_table[answer,]

    if (verbose) {
        print(sprintf("Removing %d low-count genes (%d remaining).",
                      num_before - nrow(counts), nrow(counts)))
    }
    libsize = colSums(count_table)
    counts = list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Filter low-count genes from a data set using genefilter's kOverA()
#'
#' @param counts input data frame of counts by sample
#' @param k a minimum number of samples to have >A counts
#' @param A the minimum number of counts for each gene's sample in kOverA()
#' @param verbose If set to true, prints number of genes removed / remaining
#' @return dataframe of counts without the low-count genes
#' @seealso \code{\link{genefilter}} \code{\link{kOverA}} which this uses to decide what to keep
#' @export
#' @examples
#' ## filtered_table = genefilter_kofa_counts(count_table)
genefilter_kofa_counts = function(count_table, k=1, A=1, verbose=FALSE) {
    ## genefilter has functions to work with expressionsets directly, but I think I will work merely with tables in this.
    num_before = nrow(count_table)

    if (class(count_table) == 'ExpressionSet') {
        counts = exprs(count_table)
    }
    test = kOverA(k=k, A=A)
    filter_list = filterfun(test)
    answer = genefilter(count_table, filter_list)
    count_table = count_table[answer,]

    if (verbose) {
        print(sprintf("Removing %d low-count genes (%d remaining).",
                      num_before - nrow(count_table), nrow(count_table)))
    }
    libsize = colSums(count_table)
    counts = list(count_table=count_table, libsize=libsize)
    return(counts)
}

#' Filter genes from a dataset outside a range of variance
#'
#' @param counts input data frame of counts by sample
#' @param cv_min a minimum coefficient of variance
#' @param cv_max guess
#' @param verbose If set to true, prints number of genes removed / remaining
#' @return dataframe of counts without the low-count genes
#' @seealso \code{\link{genefilter}} \code{\link{kOverA}} which this uses to decide what to keep
#' @export
#' @examples
#' ## filtered_table = genefilter_kofa_counts(count_table)
genefilter_cv_counts = function(count_table, cv_min=0.01, cv_max=1000, verbose=FALSE) {
    ## genefilter has functions to work with expressionsets directly, but I think I will work merely with tables in this.
    num_before = nrow(count_table)

    if (class(count_table) == 'ExpressionSet') {
        counts = exprs(count_table)
    }
    test = cv(cv_min, cv_max)
    filter_list = filterfun(test)
    answer = genefilter(count_table, filter_list)
    count_table = count_table[answer,]

    if (verbose) {
        print(sprintf("Removing %d low-count genes (%d remaining).",
                      num_before - nrow(count_table), nrow(count_table)))
    }
    libsize = colSums(count_table)
    counts = list(count_table=count_table, libsize=libsize)    
    return(counts)
}


testme = function() {
    ## Trying out genefilter
    library("pasilla")
    data("pasillaGenes")
    library("DESeq")
    cds = estimateSizeFactors(pasillaGenes)
    cds = estimateDispersions(cds)
    fit1 = fitNbinomGLMs(cds, count ~ type + condition)
    fit0 = fitNbinomGLMs(cds, count ~ type)
    res = data.frame(
        filterstat=rowMeans(counts(cds)),
        pvalue=nbinomGLMTest(fit1, fit0),
        row.names=featureNames(cds))

    dat = counts(cds)
        
    sfun <- coxfilter(1, 1, .05)
    ffun <- filterfun(sfun)
    l2dat = log2(dat + 1)
    l2dat = l2dat[complete.cases(l2dat),]
    which <- genefilter(l2dat, ffun)
    tt = dat[which,]
    dim(tt)

    cvfun <- cv(.5,2.5)
    ffun <- filterfun(cvfun)
    which <- genefilter(dat, ffun)
    tt = dat[which,]
    dim(tt)

    ## coefficient of variation across samples must be between 0.01 and 1000.
    cvfun = cv(0.01, 1000)
    ffun = filterfun(cvfun)
    which = genefilter(dat, ffun)
    tt = dat[which,]
    dim(tt)
    
    f1 <- kOverA(k=1, A=1)
    flist <- filterfun(f1)
    ans <- genefilter(dat, flist)
    tt = dat[ans,]
    dim(tt)

    f1 = pOverA(p=0.1, A=10) ## if > p(0.1) proportion elements in a row are > A(10), then it passes to TRUE
    flist <- filterfun(f1)
    ans <- genefilter(dat, flist)
    tt = dat[ans,]
    dim(tt)    

    cov_fun = cov(l2dat)
    af <- Anova(cov_fun, .1)
    flist = filterfun(af)
    ans = genefilter(l2dat, flist)
    tt = dat[ans,]
    dim(tt)
    
    ## nsFilter(cds) ## needs an expressionset
}    
    


#' Replace the data of an expt with normalized data
#'
#' @param expt=expt The original expt
#' @param transform="log2" The transformation desired
#'
#' @return a new expt object with normalized data and the original data saved as 'original_expressionset'
#' @export
normalize_expt = function(expt, ## The expt class passed to the normalizer
    transform="raw", norm="raw", convert="raw", batch="raw", filter_low=FALSE, ## choose the normalization strategy
    annotations=NULL, verbose=FALSE, use_original=FALSE, ## annotations used for rpkm/cpseqm, original may be used to ensure double-normalization isn't performed.
    batch1="batch", batch2=NULL, ## extra parameters for batch correction
    thresh=2, min_samples=2, p=0.01, A=1, k=1, cv_min=0.01, cv_max=1000,  ## extra parameters for low-count filtering
    ...) {
    new_expt = expt
    current = expt$expressionset
    if (is.null(new_expt$original_expressionset)) {
        new_expt$original_expressionset = new_expt$expressionset
    } else {
        print(paste("This function will replace the expt$expressionset slot with the ", transform, "(", norm, "(", convert, "))'d data.", sep=""))
        print("It saves the current data into a slot named: expt$backup_expressionset")
        print("It will also save copies of each step along the way in expt$normalized with the corresponding libsizes.")
        print("Keep the libsizes in mind when invoking limma.  The appropriate libsize is the non-log(cpm(normalized)).")
    }
    new_expt$backup_expressionset = new_expt$expressionset
    old_data = exprs(expt$original_expressionset)
    design = expt$design
    normalized = hpgl_norm(df=old_data, design=design, transform=transform, norm=norm, convert=convert, batch=batch, batch1=batch1, batch2=batch2, filter_low=filter_low, annotations=annotations, verbose=verbose, thresh=thresh, min_samples=min_samples, p=p, A=A, k=k, cv_min=cv_min, cv_max=cv_max)
    final_normalized = normalized$final_counts
    libsizes = final_normalized$libsize
    normalized_data = as.matrix(final_normalized$count_table)
    exprs(current) = normalized_data
    new_expt$normalized = normalized
    new_expt$norm_libsize = libsizes
    new_expt$expressionset = current
    new_expt$filtered = filter_low
    new_expt$transform = transform
    new_expt$norm = norm
    new_expt$convert = convert
    new_expt$batch = batch
    return(new_expt)
}


#' Normalize a dataframe/expt, express it, and/or transform it
#'
#' @param expt=expt an expt class containing all the necessary
#' metadata
#' @param df=df alternately a dataframe of counts may be used
#' @param design=design but a design dataframe must come with it
#' @param convert defines the output type which may be raw, cpm,
#' rpkm, or cp_seq_m.  Defaults to raw.
#' @param transform defines whether to log(2|10) transform the
#' data. Defaults to raw.
#' @param norm specify the normalization strategy.  Defaults to
#' raw.  This makes use of DESeq/EdgeR to provide: RLE, upperquartile,
#' size-factor, or tmm normalization.  I tend to like quantile, but there are
#' definitely corner-case scenarios for all strategies.
#' @param filter_low choose whether to low-count filter the data.
#' Defaults to true.
#' @param annotations is used for rpkm or sequence normalizations to
#' extract the lengths of sequences for normalization
#'
#' @return edgeR's DGEList expression of a count table.  This seems to
#' me to be the easiest to deal with.
#' @seealso \code{\link{cpm}}, \code{\link{rpkm}},
#' \code{\link{hpgl_rpkm}}, \code{\link{filterCounts}},
#' \code{\link{DESeqDataSetFromMatrix}},
#' \code{\link{estimateSizeFactors}}, \code{\link{DGEList}},
#' \code{\link{qNorm}}, \code{\link{calcNormFactors}}
#' 
#' @export
#' @examples
#' ## df_raw = hpgl_norm(expt=expt)  ## Only performs low-count filtering
#' ## df_raw = hpgl_norm(df=a_df, design=a_design) ## Same, but using a df
#' ## df_ql2rpkm = hpgl_norm(expt=expt, norm='quant', transform='log2', convert='rpkm')  ## Quantile, log2, rpkm
#' ## count_table = df_ql2rpkm$counts
###                                                 raw|log2|log10   sf|quant|etc  cpm|rpkm|cbcbcpm
hpgl_norm = function(df=NULL, expt=NULL, design=NULL, transform="raw", norm="raw", convert="raw", batch="raw", batch1="batch", batch2=NULL, filter_low=TRUE, annotations=NULL, verbose=FALSE, thresh=2, min_samples=2, noscale=TRUE, p=0.01, A=1, k=1, cv_min=0.01, cv_max=1000, ...) {
    lowfilter_performed = FALSE
    norm_performed = "raw"    
    convert_performed = "raw"
    transform_performed = "raw"
    batch_performed = "raw"
    if (is.null(expt) & is.null(df)) {
        stop("This needs either: an expt object containing metadata; or a df, design, and colors")
    }
    if (is.null(expt)) {
        if (verbose) {
            print("expt is null, using df and design.")
        }
        count_table = df
        expt_design = design
        column_data = colnames(count_table)
    } else if (is.null(df)) {
        count_table = as.matrix(exprs(expt$expressionset))
        expt_design = expt$design
        column_data = expt$columns
    } else {
        stop("Both df and expt are defined, choose one.")
    }

    raw_libsize = colSums(count_table)
    original_counts = list(libsize=raw_libsize, counts=count_table)
    
    ## Step 1: Perform a low count filter
    lowfiltered_counts = NULL
    if (filter_low != "FALSE") {
        if (verbose) {
            print(paste0("Filtering low counts with: ", filter_low))
        }
        if (tolower(filter_low) == "povera") {
            filter_low = "pofa"
        } else if (tolower(filter_low) == "kovera") {
            filter_low = "kofa"
        }
        if (filter_low == "cbcb") {
            lowfiltered_counts = cbcb_filter_counts(count_table, thresh=thresh, min_samples=min_samples)
            count_table = lowfiltered_counts$count_table
            lowfilter_performed = "cbcb"
        } else if (filter_low == "pofa") {
            lowfiltered_counts = genefilter_pofa_counts(count_table, p=p, A=A)
            count_table = lowfiltered_counts$count_table
            lowfilter_performed = "pofa"
        } else if (filter_low == "kofa") {
            lowfiltered_counts = genefilter_kofa_counts(count_table, k=k, A=A)
            count_table = lowfiltered_counts$count_table
            lowfilter_performed = "kofa"
        } else if (filter_low == "cv") {
            lowfiltered_counts = genefilter_cv_counts(count_table, cv_min=cv_min, cv_max=cv_max)
            count_table = lowfiltered_counts$count_table
            lowfilter_performed = "cv"
        } else {
            print("Did not recognize the filtering argument, defaulting to cbcb's.")
            print("Recognized filters are: 'cv', 'kofa', 'pofa', 'cbcb'")
            lowfiltered_counts = cbcb_filter_counts(count_table, thresh=thresh, min_samples=min_samples)
            count_table = lowfiltered_counts$count_table
            lowfilter_performed = "cbcb"            
        }
    }
    
    ## Step 2: Normalization
    ## This section handles the various normalization strategies
    ## If nothing is chosen, then the filtering is considered sufficient
    normalized_counts = NULL
    if (norm != "raw") {
        if (verbose) {
            print(paste("Applying normalization:", norm))
        }
        normalized_counts = normalize_counts(count_table, design, norm=norm)
        count_table = normalized_counts$count_table
        norm_performed = norm
    }
    
    ## Step 3: Convert the data to (likely) cpm
    ## The following stanza handles the three possible output types
    ## cpm and rpkm are both from edgeR
    ## They have nice ways of handling the log2 which I should consider
    converted_counts = NULL
    if (convert != "raw") {
        if (verbose) {
            print(paste("Setting output type as:", convert))
        }
        converted_counts = convert_counts(count_table, convert=convert)
        count_table = converted_counts$count_table
        convert_performed = convert
    }
        
    ## Step 4: Transformation
    ## Finally, this considers whether to log2 the data or no
    transformed_counts = NULL
    if (transform != "raw") {
        if (verbose) {
            print(paste("Applying: ", transform, " transformation.", sep=""))
        }
        transformed_counts = transform_counts(count_table, transform=transform, annotations=annotations, converted=convert_performed, ...)
        count_table = transformed_counts$count_table
        transform_performed = transform
    }
        
    ## Step 5: Batch correction
    batched_counts = NULL
    if (batch != "raw") {    
        if (verbose) {
            print(paste("Applying: ", batch, " batch correction(raw means nothing).", sep=""))
        }
        ## batched_counts = batch_counts(count_table, batch=batch, batch1=batch1, batch2=batch2, design=design, ...)
        batched_counts = batch_counts(count_table, batch=batch, batch1=batch1, batch2=batch2, design=design)        
        count_table = batched_counts$count_table
        batch_performed = batch
    }

    final_counts = list(count_table=count_table, libsize=colSums(count_table))
    ret_list = list(
        lowfilter_performed=lowfilter_performed, norm_performed=norm_performed, convert_performed=convert_performed,
        transform_performed=transform_performed, batch_performed=batch_performed,
        original_counts=original_counts,
        lowfiltered_counts=lowfiltered_counts,
        normalized_counts=normalized_counts,
        converted_counts=converted_counts,
        transformed_counts=transformed_counts,
        batched_counts=batched_counts,
        final_counts=final_counts,
        count_table=final_counts$count_table,
        libsize=final_counts$libsize
    )
    return(ret_list)
}

batch_counts = function(count_table, design, batch=batch, batch1=batch1, batch2=batch2 , noscale=TRUE, ...) {
    if (isTRUE(batch)) {
        batch = "limma"
    }
    if (batch == "limma") {
        batches1 = as.factor(design[, batch1])
        if (is.null(batch2)) {
            ## A reminder of removeBatchEffect usage
            ## adjusted_batchdonor = removeBatchEffect(data, batch=as.factor(as.character(des$donor)), batch2=as.factor(as.character(des$batch)))
            message("Using limma's removeBatchEffect to remove batch effect.")
            count_table = limma::removeBatchEffect(count_table, batch=batches1)
        } else {
            batches2 = as.factor(design[, batch2])
            count_table = limma::removeBatchEffect(count_table, batch=batches1, batch2=batches2)
        }
    } else if (batch == "combatmod") {
        ## message("Using a modified cbcbSeq combatMod for batch correction.")
        batches = as.factor(design[, batch1])
        conditions = as.factor(design[, "condition"])
        ## normalized_data = hpgl_combatMod(dat=data.frame(counts), batch=batches, mod=conditions, noScale=noscale, ...)
        count_table = cbcbSEQ::combatMod(dat=data.frame(count_table), batch=batches, mod=conditions, noScale=noscale, ...)
    } else if (batch == "sva") {
        batches = as.factor(design[, batch1])
        conditions = as.factor(design[,"condition"])
        df = data.frame(count_table)
        mtrx = as.matrix(df)
        conditional_model = model.matrix(~conditions, data=df)
        null_model = conditional_model[,1]
        num_surrogates = num.sv(mtrx, conditional_model)
        sva_object = sva(mtrx, conditional_model, null_model, n.sv=num_surrogates)
        mod_sv = cbind(conditional_model, sva_object$sv)
        fsva_result = fsva(mtrx, conditional_model, sva_object, newdat=mtrx, method="exact")
        ## new_expt$conditional_model = conditional_model
        ## new_expt$null_model = null_model
        ## new_expt$num_surrogates = num_surrogates
        ## new_expt$sva_object = sva_object
        ## new_expt$mod_sv = mod_sv
        ## new_expt$fsva_result = fsva_result
        count_table = fsva_result$db         
    } else {
        print("Did not recognize the batch correction, leaving the table alone.")
        print("Recognized batch corrections include: 'limma', 'combatmod', 'sva'")
    }
    libsize = colSums(count_table)
    counts = list(count_table=count_table, libsize=libsize)
    return(counts)    
}

transform_counts = function(count_table, transform="raw", converted="raw", ...) {
    if (converted != "cpm") {
        count_table = count_table + 1
    }
    if (transform == "log2") {
        count_table = log2(count_table)
    } else if (transform == "log10") {
        count_table = log10(count_table)
    } else if (transform == "log") {  ## Natural log
        count_table = log(count_table)  ## Apparently log1p does this.
    } else {
        print("Did not recognize the transformation, leaving the table alone.")
        print("Recognized transformations include: 'log2', 'log10', 'log'")
    }
    libsize = colSums(count_table)
    counts = list(count_table=count_table, libsize=libsize)
    return(counts)
}

convert_counts = function(count_table, convert="raw", annotations=NULL, ...) {
    if (convert == "edgecpm") {
        count_table = edgeR::cpm(count_table)
    } else if (convert == "cpm") {
        lib_size = colSums(count_table)
        ## count_table = t(t((count_table$counts + 0.5) / (lib_size + 1)) * 1e+06)
        transposed = t(count_table + 0.5)
        cp_counts = transposed / (lib_size + 1)
        cpm_counts = t(cp_counts * 1e+06)
        count_table = cpm_counts
    } else if (convert == "rpkm") {
        if (is.null(annotations)) {
            stop("RPKM conversion requires gene lengths.")
        }
        count_table = hpgltools::hpgl_rpkm(counts_table, annotations=annotations)
    } else if (convert == "cp_seq_m") {
        counts = edgeR::cpm(count_table)
        counts_table = hpgltools::divide_seq(counts, ...)
    }
    libsize = colSums(count_table)
    counts = list(count_table=count_table, libsize=libsize)
    return(counts)
}

normalize_counts = function(count_table, design, norm="raw") {
    if (norm == "sf") {
        ## Size-factored normalization is a part of DESeq
        original_cols = colnames(count_table)
        matrix = DESeq2::DESeqDataSetFromMatrix(countData=count_table, colData=design, design=~1)
        size_factor = BiocGenerics::estimateSizeFactors(matrix)
        count_table = BiocGenerics::counts(size_factor, normalized=TRUE)
        colnames(count_table) = original_cols
        norm_performed = "sf"
    } else if (norm == "quant") {
        # Quantile normalization (Bolstad et al., 2003)
        count_rownames = rownames(count_table)
        count_colnames = colnames(count_table)
        count_table = normalize.quantiles(as.matrix(count_table), copy=TRUE)
        rownames(count_table) = count_rownames
        colnames(count_table) = count_colnames
    } else if (norm == "tmm") {
        ## TMM normalization is documented in edgeR
        ## Set up the edgeR data structure
        count_table = edgeR::DGEList(counts=count_table)
        ## Get the tmm normalization factors
        norms = edgeR::calcNormFactors(count_table, method="TMM")
        ## Set up the DESeq data structure to which to apply the new factors
        deseq_matrix =  DESeq2::DESeqDataSetFromMatrix(countData=count_table, colData=design, design=~1)
        ## Apply the edgeR tmm factors to this
        sizeFactors(deseq_matrix) = norms$samples$norm.factors
        ## Get the counts out
        factored = BiocGenerics::counts(deseq_matrix, normalized=TRUE)
        ## return this to a DGEList
        count_table = as.matrix(factored)
        norm_performed = "tmm"
    } else if (norm == "upperquartile") {
        ## Get the tmm normalization factors
        count_table = edgeR::DGEList(counts=count_table)        
        norms = edgeR::calcNormFactors(count_table, method="upperquartile")
        ## Set up the DESeq data structure to which to apply the new factors
        deseq_matrix = DESeq2::DESeqDataSetFromMatrix(countData=count_table, colData=design, design=~1)
        ## Apply the edgeR tmm factors to this
        sizeFactors(deseq_matrix) = norms$samples$norm.factors
        ## Get the counts out
        factored = BiocGenerics::counts(deseq_matrix, normalized=TRUE)
        ## return this to a DGEList
        count_table = as.matrix(factored)
    } else if (norm == "rle") {
        ## Get the tmm normalization factors
        count_table = edgeR::DGEList(counts=count_table)
        norms = edgeR::calcNormFactors(count_table, method="RLE")
        ## Set up the DESeq data structure to which to apply the new factors
        deseq_matrix = DESeq2::DESeqDataSetFromMatrix(countData=count_table, colData=design, design=~1)
        ## Apply the edgeR tmm factors to this
        sizeFactors(deseq_matrix) = norms$samples$norm.factors
        ## Get the counts out
        factored = BiocGenerics::counts(deseq_matrix, normalized=TRUE)
        ## return this to a DGEList
        count_table = as.matrix(factored)
    } else {
        print("Did not recognize the normalization, leaving the table alone.")
        print("Recognized normalizations include: 'sf', 'quant', 'tmm', 'upperquartile', and 'rle'")
        count_table = as.matrix(count_table)
    }
    norm_libsize = colSums(count_table)
    norm_counts = list(count_table=count_table, libsize=norm_libsize)
    return(norm_counts)
}


## My root question, to which I think I have a small idea about the answer:
## Two of the very many paths to toptable()/toptags():
##  1.  normalize data -> model(~0 + condition + batch) -> limma
###     Including batch in the model loses some power, but improves the 'truth' of the result
##  2.  normalize data -> batch correction(factor(batch)) ->  model(~0 + condition) -> limma
###     Power lost before the model, also improves the 'truth' of the result.
##  Why is #1 better than #2?
### More well understood and conservative.
##  Why have we nonetheless done #2 in a few instances?  (not only because we learned that first)

#' Use a modified version of combat on some data
#'
#' @param dat a df to modify
#' @param batch a factor of batches
#' @param mod a factor of conditions
#' @param noScale the normal 'scale' option squishes the data too much, so this defaults to TRUE
#' @param prior.plots print out prior plots? FALSE
#'
#' @return a df of batch corrected data
#' @seealso \code{\link{sva}}, \code{\link{combat}},
#' 
#' @export
#' @examples
#' ## df_new = hpgl_combatMod(df, batches, model)
hpgl_combatMod = function (dat, batch, mod, noScale=TRUE, prior.plots=FALSE) {
    par.prior = TRUE
    numCovs = NULL
    mod = cbind(mod, batch)
    check = apply(mod, 2, function(x) all(x == 1))
    mod = as.matrix(mod[, !check])
    colnames(mod)[ncol(mod)] = "Batch"
    if (sum(check) > 0 & !is.null(numCovs)) 
        numCovs = numCovs - 1
    design <- sva:::design.mat(mod, numCov = numCovs)
    batches <- sva:::list.batch(mod)
    n.batch <- length(batches)
    n.batches <- sapply(batches, length)
    n.array <- sum(n.batches)
    NAs = any(is.na(dat))
    if (NAs) {
        cat(c("Found", sum(is.na(dat)), "Missing Data Values\n"), 
            sep = " ")
    }
    cat("Standardizing Data across genes\n")
    if (!NAs) {
        B.hat <- solve(t(design) %*% design) %*% t(design) %*% 
            t(as.matrix(dat))
    }
    else {
        B.hat = apply(dat, 1, Beta.NA, design)
    }
    grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]
    if (!NAs) {
        var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array, 
            n.array)
    }
    else {
        var.pooled <- apply(dat - t(design %*% B.hat), 1, var, 
            na.rm = T)
    }
    stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
    if (!is.null(design)) {
        tmp <- design
        tmp[, c(1:n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% B.hat)
    }
    s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, 
        n.array)))
    if (noScale) {
        m.data <- dat - stand.mean
        mse <- ((dat - t(design %*% B.hat))^2) %*% rep(1/(n.array - 
            ncol(design)), n.array)
        hld <- NULL
        bayesdata <- dat
        for (k in 1:n.batch) {
            cat(paste("Fitting 'shrunk' batch ", k, " effects\n", 
                sep = ""))
            sel <- batches[[k]]
            gammaMLE <- rowMeans(m.data[, sel])
            mprior <- mean(gammaMLE, na.rm = TRUE)
            vprior <- var(gammaMLE, na.rm = TRUE)
            prop <- vprior/(mse/(length(sel)) + vprior)
            gammaPost <- prop * gammaMLE + (1 - prop) * mprior
            for (i in sel) {
                bayesdata[, i] <- bayesdata[, i] - gammaPost
            }
            stats <- data.frame(gammaPost = gammaPost, gammaMLE = gammaMLE, 
                prop = prop)
            hld[[paste("Batch", k, sep = ".")]] <- list(stats = stats, 
                indices = sel, mprior = mprior, vprior = vprior)
        }
        cat("Adjusting data for batch effects\n")
        return(bayesdata)
    }
    else {
        cat("Fitting L/S model and finding priors\n")
        batch.design <- design[, 1:n.batch]
        if (!NAs) {
            gamma.hat <- solve(t(batch.design) %*% batch.design) %*% 
                t(batch.design) %*% t(as.matrix(s.data))
        }
        else {
            gamma.hat = apply(s.data, 1, Beta.NA, batch.design)
        }
        delta.hat <- NULL
        for (i in batches) {
            delta.hat <- rbind(delta.hat, apply(s.data[, i], 
                1, var, na.rm = T))
        }
        gamma.bar <- apply(gamma.hat, 1, mean)
        t2 <- apply(gamma.hat, 1, var)
        a.prior <- apply(delta.hat, 1, sva:::aprior)
        b.prior <- apply(delta.hat, 1, sva:::bprior)
        if (prior.plots & par.prior) {
            par(mfrow = c(2, 2))
            tmp <- density(gamma.hat[1, ])
            plot(tmp, type = "l", main = "Density Plot")
            xx <- seq(min(tmp$x), max(tmp$x), length = 100)
            lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
            qqnorm(gamma.hat[1, ])
            qqline(gamma.hat[1, ], col = 2)
            tmp <- density(delta.hat[1, ])
            invgam <- 1/rgamma(ncol(delta.hat), a.prior[1], b.prior[1])
            tmp1 <- density(invgam)
            plot(tmp, typ = "l", main = "Density Plot", ylim = c(0, 
                max(tmp$y, tmp1$y)))
            lines(tmp1, col = 2)
            qqplot(delta.hat[1, ], invgam, xlab = "Sample Quantiles", 
                ylab = "Theoretical Quantiles")
            lines(c(0, max(invgam)), c(0, max(invgam)), col = 2)
            title("Q-Q Plot")
        }
        gamma.star <- delta.star <- NULL
        if (par.prior) {
            cat("Finding parametric adjustments\n")
            for (i in 1:n.batch) {
                temp <- sva:::it.sol(s.data[, batches[[i]]], 
                  gamma.hat[i, ], delta.hat[i, ], gamma.bar[i], 
                  t2[i], a.prior[i], b.prior[i])
                gamma.star <- rbind(gamma.star, temp[1, ])
                delta.star <- rbind(delta.star, temp[2, ])
            }
        }
        else {
            cat("Finding nonparametric adjustments\n")
            for (i in 1:n.batch) {
                temp <- sva:::int.prior(as.matrix(s.data[, batches[[i]]]), 
                  gamma.hat[i, ], delta.hat[i, ])
                gamma.star <- rbind(gamma.star, temp[1, ])
                delta.star <- rbind(delta.star, temp[2, ])
            }
        }
        cat("Adjusting the Data\n")
        bayesdata <- s.data
        j <- 1
        for (i in batches) {
            bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, 
                ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% 
                t(rep(1, n.batches[j])))
            j <- j + 1
        }
        bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, 
            n.array)))) + stand.mean
        return(bayesdata)
    }
}
