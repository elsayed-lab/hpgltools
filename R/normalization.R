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
    cpm
    t(log2(t(counts + 0.5) / (colSums(counts) + 1) * 1e+06))
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
#' @param thresh lower threshold of counts (4 by default)
#' @param minSamples minimum number of samples
#' @return dataframe of counts without the low-count genes
#' @seealso \code{\link{log2CPM}} which this uses to decide what to keep
#' @export
#' @examples
#' ## filtered_table = filter_counts(count_table)
filter_counts = function(counts, thresh=2, min_samples=2) {
    ## I think having a log2cpm here is kind of weird, because the next step in processing is to cpm the data.
    ##cpms = 2^log2CPM(counts, lib.size=lib.size)$y
    ## cpms = 2^hpgl_log2cpm(counts)
    keep = rowSums(counts > thresh) >= min_samples
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
normalize_expt = function(expt, transform="log2", norm="quant", convert="cpm", filter_low=TRUE, annotations=NULL, verbose=FALSE, use_original=TRUE, thresh=2, min_samples=2, batch=NULL, batch1="batch", batch2=NULL, ...) {
    new_expt = expt
    if (is.null(new_expt$original_expressionset)) {
        new_expt$original_expressionset = new_expt$expressionset
    } else {
        print(paste("This function defaults to replacing the expt$expressionset slot with the ", transform, "(", norm, "(", convert, "))'d data.", sep=""))
        print("It saves the current data into a slot named: expt$backup_expressionset")
    }
    new_expt$backup_expressionset = new_expt$expressionset
    old_data = exprs(expt$original_expressionset)
    design = expt$design

    normalized_data = hpgl_norm(df=old_data, design=design, transform=transform, norm=norm, convert=convert, filter_low=filter_low, annotations=annotations, verbose=verbose, thresh=thresh, min_samples=min_samples)

    if (is.null(batch)) {
        exprs(new_expt$expressionset) = as.matrix(normalized_data$counts)
    } else {
        if (batch == "limma") {
            batches1 = as.factor(design[,batch1])
            if (is.null(batch2)) {
                ## A reminder of removeBatchEffect usage
                ## adjusted_batchdonor = removeBatchEffect(data, batch=as.factor(as.character(des$donor)), batch2=as.factor(as.character(des$batch)))
                message("Using limma's removeBatchEffect to remove batch effect.")
                normalized_data = removeBatchEffect(normalized_data$counts, batch=batches1)
            } else {
                batches2 = as.factor(design[,batch2])
                normalized_data = removeBatchEffect(normalized_data$counts, batch=batches1, batch2=batches2)
            }
        } else if (batch == "combatmod") {
            message("Using cbcbSeq's combatMod for batch correction.")
            batches = as.factor(design[,"batch"])
            conditions = as.factor(design[,"condition"])
            df = data.frame(normalized_data$counts)
            normalized_data = my_combatMod(dat=df, batch=batches, mod=conditions, ...)
            
        } else if (batch == "sva") {

            batches = as.factor(design[,"batch"])
            conditions = as.factor(design[,"condition"])
            df = data.frame(normalized_data$counts)
            conditional_model = model.matrix(~conditions, data=df)
            null_model = conditional_model[,1]
            num_surrogates = num.sv(as.matrix(df), conditional_model)
            sva_object = sva(as.matrix(df), conditional_model, null_model, n.sv=num_surrogates)
            mod_sv = cbind(conditional_model, sva_object$sv)
            fsva_result = fsva(as.matrix(df), conditional_model, sva_object, newdat=as.matrix(df), method="exact")
            new_expt$conditional_model = conditional_model
            new_expt$null_model = null_model
            new_expt$num_surrogates = num_surrogates
            new_expt$sva_object = sva_object
            new_expt$mod_sv = mod_sv
            new_expt$fsva_result = fsva_result
            normalized_data = fsva_result$db            
##            normalized_data = combatMod(normalized_data, batches, conditions, noScale=TRUE)
        } else {
            message("haven't implemented other batch removals, just doing limma's removeBatchEffect()")
            normalized_data = removeBatchEffect(normalized_data, batch=get(batch1))            
        }
        exprs(new_expt$expressionset) = as.matrix(normalized_data)
    } ## End the if/else batch correction
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
#' ## df_ql2rpkm = hpgl_norm(expt=expt, norm_type='quant', filter='log2', out_type='rpkm'  ## Quantile, log2, rpkm
#' ## count_table = df_ql2rpkm$counts
###                                                 raw|log2|log10   sf|quant|etc  cpm|rpkm
hpgl_norm = function(df=NULL, expt=NULL, design=NULL, transform="raw", norm="raw", convert="raw", filter_low=TRUE, annotations=NULL, verbose=FALSE, thresh=2, min_samples=2, ...) {
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

    ## Step 1: Perform a low count filter
    if (filter_low == TRUE) {
        if (verbose) {
            print("Filtering low counts")
        }
        original_dim = dim(count_table)
        count_table = as.matrix(filter_counts(count_table, thresh=thresh, min_samples=min_samples))
        if (verbose) {
            following_dim = dim(count_table)
            lost_rows = original_dim[1] - following_dim[1]
            print(paste("Low count filtering cost:", lost_rows, "gene(s)."))
        }
    }

    ## Step 2: Normalization
    ## This section handles the various normalization strategies
    ## If nothing is chosen, then the filtering is considered sufficient
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
        # Quantile normalization (Bolstad et al., 2003)
        count_rownames = rownames(count_table)
        count_colnames = colnames(count_table)
        count_table = normalize.quantiles(as.matrix(count_table))
        rownames(count_table) = count_rownames
        colnames(count_table) = count_colnames
        # Convert to a DGEList
        count_table = edgeR::DGEList(counts=count_table)
    } else if (norm == "tmm") {
        ## TMM normalization is documented in edgeR
        ## Set up the edgeR data structure
        count_table = edgeR::DGEList(counts=count_table)
        ## Get the tmm normalization factors
        norms = edgeR::calcNormFactors(count_table, method="TMM")
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
    
    ## Step 3: Convert the data to (likely) cpm
    ## The following stanza handles the three possible output types
    ## cpm and rpkm are both from edgeR
    ## They have nice ways of handling the log2 which I should consider
    if (verbose) {
        print(paste("Setting output type as:", out_type))
    }
    if (convert == "cpm") {
        counts = edgeR::cpm(count_table)
        count_table = edgeR::DGEList(counts=count_table)
    } else if (convert == "rpkm") {
        if (is.null(annotations)) {
            stop("RPKM conversion requires gene lengths.")
        }
        counts = hpgltools::hpgl_rpkm(counts_table, annotations=annotations)
        count_table = edgeR::DGEList(counts=counts)
    } else if (convert == "cp_seq_m") {
        counts = edgeR::cpm(count_table)
        counts = hpgltools::divide_seq(counts, ...)
        count_table = edgeR::DGEList(counts=counts)
    } else {
        count_table = edgeR::DGEList(counts=count_table)
    }

    ## Step 4: Transformation
    ## Finally, this considers whether to log2 the data or no
    if (verbose) {
        print(paste("Applying: ", transform, " transformation.", sep=""))
    }
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


## The original combatMod from cbcbSeq
## I intend to rewrite sections of it to make it a bit more robust to situations when we don't have perfect batch/conditions
combatMod = function (dat, batch, mod, noScale = TRUE, prior.plots = FALSE) {
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


my_combatMod = function (dat, batch, mod, noScale = TRUE, prior.plots = FALSE) {
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
