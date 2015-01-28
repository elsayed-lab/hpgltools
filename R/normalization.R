#' Express a data frame of counts as reads per killobase(gene) per
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
hpgl_log2cpm = function(counts) {
    t(log2(t(counts + 0.5) / (colSums(counts) + 1) * 1e+06))
}

divide_seq = function(counts, pattern="TA", fasta="testme.fasta", gff="testme.fasta", entry_type="gene") {
    ## Testing parameters
    ##counts=exprs(rnarpf_prometa_kexpt$expressionset)
    ##pattern="ATG"
    ##fasta="reference/lmajor_genome/TriTrypDB-6.0_LmajorFriedlin_Genome.fasta"
    ##gff="reference/lmajor_genome/TriTrypDB-6.0_LmajorFriedlin_genes.gff"
    ##entry_type="gene"
    ##
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
#' @param thresh lower threshold of counts (4 by default)
#' @param minSamples minimum number of samples
#' @return dataframe of counts without the low-count genes
#' @seealso \code{\link{log2CPM}} which this uses to decide what to keep
#' @export
#' @examples
#' ## filtered_table = filter_counts(count_table)
filter_counts = function(counts, thresh=4, minSamples=2) {
    cpms = 2^hpgl_log2cpm(counts)
    keep = rowSums(cpms > thresh) >= minSamples
    counts = counts[keep,]
    return(counts)
}

#' Replace the data of an expt with normalized data
#'
#' @param expt=expt The original expt
#' @param transform="log2" The transformation desired
#'
#' @return a new expt object with normalized data and the original data saved as 'original_expressionset'
#' @export
normalize_expt = function(expt, transform="log2", norm="quant", convert="cpm", filter_low=TRUE, annotations=NULL, verbose=FALSE, use_original=TRUE, ...) {
    new_expt = expt
    if (is.null(new_expt$original_expressionset)) {
        new_expt$original_expressionset = new_expt$expressionset
    } else {
        print("This function defaults to using the original expressionset for normalization.")
    }
    new_expt$backup_expressionset = new_expt$expressionset
    old_data = exprs(expt$original_expressionset)
    normalized_data = hpgl_norm(df=old_data, design=expt$design, transform=transform, norm=norm, convert=convert, filter_low=filter_low, annotations=annotations, verbose=verbose)
    exprs(new_expt$expressionset) = as.matrix(normalized_data$counts)
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
#' raw.  This makes use of DESeq/EdgeR/cbcbSEQ to provide: quantile,
#' RLE, upperquartile, size-factor, or tmm normalization.  I tend to
#' like quantile, but there are definitely corner-case scenarios for
#' all strategies.
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
#' df_raw = hpgl_norm(expt=expt)  ## Only performs low-count filtering
#' df_raw = hpgl_norm(df=a_df, design=a_design) ## Same, but using a df
#' df_ql2rpkm = hpgl_norm(expt=expt, norm_type='quant', filter='log2', out_type='rpkm'  ## Quantile, log2, rpkm
#' count_table = df_ql2rpkm$counts
###                                                 raw|log2|log10   sf|quant|etc  cpm|rpkm
hpgl_norm = function(df=NULL, expt=NULL, design=NULL, transform="raw", norm="raw", convert="raw", filter_low=TRUE, annotations=NULL, verbose=FALSE, ...) {
    ## Testing args
    ##df=NULL
    ##expt=rnarpf_prometa_kexpt
    ##norm_type="sf"
    ##filter="log2"
    ##out_type="cpm_seq_m"
    ##verbose=TRUE
    ## End test args
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
    if (filter_low == TRUE) {
        if (verbose) {
            print("Filtering low counts")
        }
        original_dim = dim(count_table)
        count_table = as.matrix(filter_counts(count_table))
        following_dim = dim(count_table)
        lost_rows = original_dim[1] - following_dim[1]
        if (verbose) {
            print(paste("Low count filtering cost:", lost_rows, "gene(s)."))
        }
    }
    ### This section handles the various normalization strategies
    ### If nothing is chosen, then the filtering is considered sufficient
    if (verbose) {
        print(paste("Applying normalization:", norm_type))
    }
    if (norm == "sf") {
        ## Size-factored normalization is a part of DESeq
        matrix = DESeq2::DESeqDataSetFromMatrix(countData=count_table, colData=expt_design, design=~1)
        size_factor = BiocGenerics::estimateSizeFactors(matrix)
        count_table = BiocGenerics::counts(size_factor, normalized=TRUE)
        colnames(count_table) = rownames(column_data)
        count_table = edgeR::DGEList(counts=count_table)
    } else if (norm == "quant") {
        ## Quantile normalization is documented in Kwame's cbcbSEQ
        count_table = cbcbSEQ::qNorm(count_table)
        count_table = edgeR::DGEList(counts=count_table)
    } else if (norm == "tmm") {
        ## TMM normalization is documented in edgeR
        ## Set up the edgeR data structure, this is used for TMM
        count_table = edgeR::DGEList(counts=count_table)
        ## Get the tmm normalization factors
        norms = edgeR::calcNormFactors(count_table)
        ## Set up the DESeq data structure to which to apply the new factors
        deseq_matrix =  DESeq2::DESeqDataSetFromMatrix(countData=count_table, colData=expt_design, design=~1)
        ## Apply the edgeR tmm factors to this
        sizeFactors(deseq_matrix) = norms$samples$norm.factors
        ## Get the counts out
        factored = BiocGenerics::counts(deseq_matrix, normalized=TRUE)
        ## return this to a DGEList
        count_table = edgeR::DGEList(counts=factored)        
    } else if (norm == "upperquartile") {
        ## Get the tmm normalization factors
        count_table = edgeR::DGEList(counts=count_table)        
        norms = edgeR::calcNormFactors(count_table, method="upperquartile")
        ## Set up the DESeq data structure to which to apply the new factors
        deseq_matrix = DESeq2::DESeqDataSetFromMatrix(countData=count_table, colData=expt_design, design=~1)
        ## Apply the edgeR tmm factors to this
        sizeFactors(deseq_matrix) = norms$samples$norm.factors
        ## Get the counts out
        factored = BiocGenerics::counts(deseq_matrix, normalized=TRUE)
        ## return this to a DGEList
        count_table = edgeR::DGEList(counts=factored)
    } else if (norm == "rle") {
        ## Get the tmm normalization factors
        count_table = edgeR::DGEList(counts=count_table)
        norms = edgeR::calcNormFactors(count_table, method="RLE")
        ## Set up the DESeq data structure to which to apply the new factors
        deseq_matrix = DESeq2::DESeqDataSetFromMatrix(countData=count_table, colData=expt_design, design=~1)
        ## Apply the edgeR tmm factors to this
        sizeFactors(deseq_matrix) = norms$samples$norm.factors
        ## Get the counts out
        factored = BiocGenerics::counts(deseq_matrix, normalized=TRUE)
        ## return this to a DGEList
        count_table = edgeR::DGEList(counts=factored)
    } else {
        count_table = edgeR::DGEList(counts=count_table)
    }
    ### The following stanza handles the three possible output types
    ### cpm and rpkm are both from edgeR
    ## They have nice ways of handling the log2 which I should consider
    if (verbose) {
        print(paste("Setting output type as:", out_type))
    }
    if (convert == "cpm") {
        counts = count_table$counts
        counts = edgeR::cpm(counts)
        count_table = edgeR::DGEList(counts=count_table)
    } else if (convert == "rpkm") {
        counts = count_table$counts
        if (is.null(annotations)) {
            stop("RPKM conversion requires gene lengths.")
        }
        counts = hpgltools::hpgl_rpkm(counts, annotations=annotations)
        count_table = edgeR::DGEList(counts=counts)
    } else if (convert == "cp_seq_m") {
        counts = count_table$counts
        counts = edgeR::cpm(counts)
        counts = hpgltools::divide_seq(counts, ...)
        count_table = edgeR::DGEList(counts=counts)
    } else {
        count_table = edgeR::DGEList(counts=count_table$counts)
    }
    if (verbose) {
        print(paste("Applying:", filter, "filter"))
    }
    ### Finally, this considers whether to log2 the data or no
    counts = count_table$counts
    if (transform == "log2") {
        counts = log2(counts + 1)
    } else if (transform == "log10") {
        counts = log10(counts + 1)
    } else if (transform == "log") {  ## Natural log
        counts = log(counts + 1)  ## Apparently log1p does this.
    }
    count_table = DGEList(counts=counts)
    return(count_table)
}
