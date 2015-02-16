## Time-stamp: <Mon Feb 16 16:51:54 2015 Ashton Trey Belew (abelew@gmail.com)>
### count_tables.R contains some functions for simple count-table manipulations
### This includes reading in files, creating expressionSets, and subsets.

#' Wrap bioconductor's expressionset to include some other extraneous
#' information.  This simply calls create_experiment and then does
#' expt_subset for everything
#'
#' @param file a comma separated file describing the samples with
#' information like condition,batch,count_filename,etc
#' @param color_hash a hash which describes how to color the samples
#' 
#' @return  experiment an expressionset
#' @seealso \code{\link{pData}}, \code{\link{fData}},
#' \code{\link{exprs}}, \code{\link{hpgl_read_files}},
#' \code{\link{as.list.hash}}
#' 
#' @export
#' @examples
#' ## new_experiment = create_experiment("some_csv_file.csv", color_hash)
#' ## Dude, you need to remember that this depends on an existing data structure of
#' ## gene annotations.
create_expt = function(file, color_hash=NULL, suffix=".count.gz", header=FALSE, genes=NULL, by_type=FALSE, by_sample=FALSE, sep=",", count_dataframe=NULL, savefile="expt", ...) {
    tmp_definitions = read.csv(file=file, comment.char="#", sep=sep)
    colnames(tmp_definitions) = tolower(colnames(tmp_definitions))    
    ## Sometimes, R adds extra rows on the bottom of the data frame using this command.
    ## Thus the next line
    print("This function needs the conditions and batches to be an explicit column in the sample sheet.")
    tmp_definitions = subset(tmp_definitions, sample.id != "")
    condition_names = unique(tmp_definitions$condition)
    num_colors = length(condition_names)
    colors = colorRampPalette(brewer.pal(num_colors,"Dark2"))(num_colors)
    color_hash = hash(keys=as.character(condition_names), values=colors)
    expt_list = create_experiment(file, color_hash, suffix=suffix, header=header, genes=genes, by_type=by_type, by_sample=by_sample, count_dataframe=count_dataframe, sep=sep, low_files=low_files)
    expt = expt_list$expt
    def = expt_list$def
    new_expt = expt_subset(expt, "")
    new_expt$definitions = def
    if (savefile) {
        save(list = c("new_expt"), file=paste(savefile, ".Rdata", sep=""))
    }
    return(new_expt)
}

#' Wrap bioconductor's expressionset to include some other extraneous
#' information.
#'
#' @param file a comma separated file describing the samples with
#' information like condition,batch,count_filename,etc
#' @param color_hash a hash which describes how to color the samples
#' 
#' @return  experiment an expressionset
#' @seealso \code{\link{pData}}, \code{\link{fData}},
#' \code{\link{exprs}}, \code{\link{hpgl_read_files}},
#' \code{\link{as.list.hash}}
#' 
#' @export
#' @examples
#' ## new_experiment = create_experiment("some_csv_file.csv", color_hash)
#' ## Dude, you need to remember that this depends on an existing data structure of
#' ## gene annotations.
create_experiment = function(file, color_hash, suffix=".count.gz", header=FALSE, genes=NULL, by_type=FALSE, by_sample=FALSE, count_dataframe=NULL, sep=",", ...) {
    print("Please note that thus function assumes a specific set of columns in the sample sheet:")
    print("The most important ones are: Sample.ID, Stage, Type.")
    print("Other columns it will attempt to create by itself, but if")
    print("batch and condition are provided, that is a nice help.")
    sample_definitions = read.csv(file=file, comment.char="#", sep=sep)
    colnames(sample_definitions) = tolower(colnames(sample_definitions))
    sample_definitions = sample_definitions[grepl('(^HPGL|^hpgl)', sample_definitions$sample.id, perl=TRUE),]
    if (is.null(sample_definitions$condition)) {
        sample_definitions$condition = tolower(paste(sample_definitions$type, sample_definitions$stage, sep="_"))
        sample_definitions$batch = gsub("\\s+|\\d+|\\*", "", sample_definitions$batch, perl=TRUE)
    }
    design_colors_list = as.list.hash(color_hash)
    sample_definitions$colors = as.list(design_colors_list[as.character(sample_definitions$condition)])
    ##sample_definitions$colors = as.list(color_hash[sample_definitions$condition])
    ## The logic here is that I want by_type to be the default, but only
    ## if no one chooses either.
    if (!isTRUE(by_type) & !isTRUE(by_sample)) {
        by_type = TRUE
    }
    if (isTRUE(by_type)) {
        sample_definitions$counts = paste("processed_data/count_tables/", tolower(sample_definitions$type),
            "/", tolower(sample_definitions$stage), "/",
            sample_definitions$sample.id, suffix, sep="")
        sample_definitions$intercounts = paste("data/count_tables/", tolower(sample_definitions$type),
            "/", tolower(sample_definitions$stage), "/",
            sample_definitions$sample.id, "_inter", suffix, sep="")
    }
    if (isTRUE(by_sample)) {
        sample_definitions$counts = paste("processed_data/count_tables/", as.character(sample_definitions$sample.id), "/", as.character(sample_definitions$sample.id), suffix, sep="")
        sample_definitions$intercounts = paste("processed_data/count_tables/", as.character(sample_definitions$sample.id), "/", as.character(sample_definitions$sample.id), "_inter", suffix, sep="")
    }
    sample_definitions = as.data.frame(sample_definitions)
    rownames(sample_definitions) = sample_definitions$sample.id
    if (is.null(count_dataframe)) {
        all_count_tables = hpgl_read_files(as.character(sample_definitions$sample.id),
            as.character(sample_definitions$counts), header=header, suffix=suffix)
    } else {
        all_count_tables = count_dataframe
        colnames(all_count_tables) = rownames(sample_definitions)
    }
    all_count_matrix = data.frame(all_count_tables)
    rownames(all_count_matrix) = gsub("^exon:","", rownames(all_count_matrix))
    rownames(all_count_matrix) = make.names(gsub(":\\d+","", rownames(all_count_matrix)), unique=TRUE)    
    if (is.null(genes)) {
        gene_info = data.frame(all_count_matrix)
    } else {
        if (is.null(genes$ID)) {
            genes$ID = rownames(genes)
        }
        gene_info = genes[genes$ID %in% rownames(all_count_matrix),]
        all_count_matrix = all_count_matrix[rownames(all_count_matrix) %in% genes$ID,]                
    }
    if (is.null(sample_definitions$stage)) {
        sample_definitions$stage = "unknown"
    }
    if (is.null(sample_definitions$type)) {
        sample_definitions$type = "unknown"
    }
    if (is.null(sample_definitions$condition)) {
        sample_definitions$condition = "unknown"
    }
    if (is.null(sample_definitions$batch)) {
        sample_definitions$batch = "unknown"
    }
    if (is.null(sample_definitions$colors)) {
        sample_definitions$colors = "unknown"
    }
    if (is.null(sample_definitions$counts)) {
        sample_definitions$counts = "unknown"
    }
    if (is.null(sample_definitions$intercounts)) {
        sample_definitions$intercounts = "unknown"
    }
    metadata = new("AnnotatedDataFrame", data.frame(sample=as.character(sample_definitions$sample.id),
        stage=as.character(sample_definitions$stage),
        type=as.character(sample_definitions$type),
        condition=as.character(sample_definitions$condition),
        batch=as.character(sample_definitions$batch),
        color=as.character(sample_definitions$colors),
        counts=sample_definitions$counts,
        intercounts=sample_definitions$intercounts))
    sampleNames(metadata) = colnames(all_count_matrix)
    feature_data = new("AnnotatedDataFrame", gene_info)
    featureNames(feature_data) = rownames(all_count_matrix)
    experiment = new("ExpressionSet", exprs=all_count_matrix,
        phenoData=metadata, featureData=feature_data)
    ret = list(expt=experiment, def=sample_definitions)
    return(ret)
}

#' Extract a subset of samples following some rule(s) from an
#' experiment class
#'
#' @param expt an expt which is a home-grown class containing an
#' expressionSet, design, colors, etc.
#' @param subset a valid R expression which defines a subset of the
#' design to keep.
#' 
#' @return  metadata an expt class which contains the smaller set of
#' data
#' @seealso \code{\link{ExpressionSet}}, \code{\link{pData}},
#' \code{\link{exprs}}, and \code{\link{fData}}
#' 
#' @export
#' @examples
#' ## smaller_expt = expt_subset(big_expt, "condition=='control'")
#' ## all_expt = expt_subset(expressionset, "")  ## extracts everything
expt_subset = function(expt, subset, by_definitions=FALSE) {
    if (class(expt) == "ExpressionSet") {
        expressionset = expt
    } else if (class(expt) == "expt") {
        expressionset = expt$expressionset
    } else {
        stop("expt is neither an expt nor ExpressionSet")
    }
    if (isTRUE(by_definitions)) {
        initial_metadata = expt$definitions
    } else {
        initial_metadata = Biobase::pData(expressionset)
    }
    r_expression=paste("subset(initial_metadata,", subset, ")")
    samples = eval(parse(text=r_expression))
    ## design = data.frame(sample=samples$sample, condition=samples$condition, batch=samples$batch)
    design = as.data.frame(samples)
    ## This is to get around stupidity with respect to needing all factors to be in a DESeqDataSet
    conditions = as.factor(as.character(design$condition))
    batches = as.factor(as.character(design$batch))
    design$condition = conditions
    design$batch = batches
    samplenames = as.character(samples$sample)
    colors = as.character(samples$color)
    names = paste(conditions, batches, sep="-")
    expressionset = expressionset[, sampleNames(expressionset) %in% samplenames]
    columns = data.frame(sample=colnames(exprs(expressionset)))
    rownames(columns) = colnames(exprs(expressionset))
    metadata = list(initial_metadata=initial_metadata,
        original_expressionset=expressionset,
        expressionset=expressionset,
        samples=samples,
        design=design,
        stages=initial_metadata$stage,
        types=initial_metadata$type,
        conditions=conditions,
        batches=batches,
        samplenames=samplenames,
        colors=colors,
        names=names,
        columns=columns)
    class(metadata) = "expt"
    return(metadata)
}

#' Read a bunch of count tables and create a usable data frame from
#' them.
#'
#' @param ids a list of experimental ids
#' @param files a list of files to read
#' @param header whether or not the count tables include a header row
#'        (default: FALSE)
#' @param include_summary_rows whether HTSeq summary rows should be included
#'        (default: FALSE)
#' 
#' @return  count_table a data frame of count tables
#' @seealso \code{\link{create_experiment}}
#' 
#' @export
#' @examples
#' ## count_tables = hpgl_read_files(as.character(sample_ids), as.character(count_filenames))
hpgl_read_files = function(ids, files, header=FALSE, include_summary_rows=FALSE, suffix=NULL, ...) {
    ## load first sample
    lower_filenames = files
    dirs = dirname(lower_filenames)
    low_files = tolower(basename(files))
    if (!is.null(suffix)) {
        low_hpgl = gsub(suffix, "", basename(files))
        low_hpgl = tolower(low_hpgl)
        low_hpgl = paste(low_hpgl, suffix, sep="")
    } else {
        low_hpgl = gsub("HPGL","hpgl", basename(files))
    }
    lower_filenames = paste(dirs, low_files, sep="/")
    lowhpgl_filenames = paste(dirs, low_hpgl, sep="/")

    if (file.exists(tolower(files[1]))) {
        files = tolower(files)
    } else if (file.exists(lowhpgl_filenames[1])) {
        files = lowhpgl_filenames
    } else if (file.exists(lower_filenames[1])) {
        files = lower_filenames
    }
    ##count_table = read.table(files[1], header=header, ...)
    count_table = read.table(files[1], header=header)    
    colnames(count_table) = c("ID", ids[1])
    ## iterate over and append remaining samples
    for (table in 2:length(files)) {
        tmp_count = read.table(files[table], header=header)
        colnames(tmp_count) = c("ID", ids[table])
        count_table = merge(count_table, tmp_count, by="ID")
    }
    
    rm(tmp_count)
    ## set row and columns ids
    rownames(count_table) = count_table$ID
    count_table = count_table[-1]
    colnames(count_table) = ids

    ## remove summary fields added by HTSeq
    if (!include_summary_rows) {
        htseq_meta_rows = c('__no_feature', '__ambiguous', '__too_low_aQual',
                            '__not_aligned', '__alignment_not_unique')
        count_table = count_table[!rownames(count_table) %in% htseq_meta_rows,]
    }
    return(count_table)
}
