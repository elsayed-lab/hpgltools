## Time-stamp: <Tue Dec 16 16:35:51 2014 Ashton Trey Belew (abelew@gmail.com)>
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
#' \code{\link{exprs}}, \code{\link{my_read_files}},
#' \code{\link{as.list.hash}}
#' 
#' @export
#' @examples
#' ## new_experiment = create_experiment("some_csv_file.csv", color_hash)
#' ## Dude, you need to remember that this depends on an existing data structure of
#' ## gene annotations.
create_expt = function(file, color_hash=NULL, suffix=".count.gz", header=FALSE, genes=tooltip_data, by_type=FALSE, by_sample=FALSE) {
    ## Test params
    ##file = "all_samples.csv"
    ##suffix=".htseq.gz"
    ##header=FALSE
    ##by_type=TRUE
    ## End test params
    if (is.null(color_hash)) {
        tmp_definitions = read.csv(file=file, comment.char="#")
        ## Sometimes, R adds extra rows on the bottom of the data frame using this command.
        ## Thus the next line
        tmp_definitions = subset(tmp_definitions, Sample.ID != "")
        condition_names = unique(tmp_definitions$condition)
        num_colors = length(condition_names)        
        colors = colorRampPalette(brewer.pal(num_colors,"Dark2"))(num_colors)
        color_hash = hash(keys=as.character(condition_names), values=colors)
    }
    expt = create_experiment(file, color_hash, suffix=suffix, header=header, genes=genes, by_type=by_type, by_sample=by_sample)
    expt = expt_subset(expt, "")
    data_rows = nrow(exprs(expt$expressionset))
    data_cols = ncol(exprs(expt$expressionset))
    ## EdgeR makes toy data with this call:
    ##bcv = 0.2 ## (dispersion, 0.4 for data, 0.1 for very similar, 0.01 for replicates)
    ##data_cols = 20
    ##data_rows = 2
    ##counts <- matrix( rnbinom(40,size=1/bcv^2,mu=10), data_cols, data_rows)
    ##expt$sample_data = counts(make_exampledata(ngenes=data_rows, columns=data_cols))
    return(expt)
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
#' \code{\link{exprs}}, \code{\link{my_read_files}},
#' \code{\link{as.list.hash}}
#' 
#' @export
#' @examples
#' ## new_experiment = create_experiment("some_csv_file.csv", color_hash)
#' ## Dude, you need to remember that this depends on an existing data structure of
#' ## gene annotations.
create_experiment = function(file, color_hash, suffix=".count.gz", header=FALSE, genes=tooltip_data, by_type=FALSE, by_sample=FALSE) {
    ## Testing parameters
    ##file="human_samples.csv"
    ##color_hash=color_hash
    ##suffix="_hg19.count.gz"
    ##header=TRUE
    ## End testing parameters
    print("Please note that thus function assumes a specific set of columns in the sample sheet:")
    print("The most important ones are: Sample.ID, Stage, Type.")
    print("Other columns it will attempt to create by itself, but if")
    print("batch and condition are provided, that is a nice help.")
    sample_definitions = read.csv(file=file, comment.char="#")
    sample_definitions = sample_definitions[grepl('^HPGL', sample_definitions$Sample.ID, perl=TRUE),]
    if (is.null(sample_definitions$condition)) {
        sample_definitions$condition = tolower(paste(sample_definitions$Type, sample_definitions$Stage, sep="_"))
        sample_definitions$batch = gsub("\\s+|\\d+|\\*", "", sample_definitions$Batch, perl=TRUE)
    }
    design_colors_list = as.list.hash(color_hash)
    sample_definitions$colors = as.list(design_colors_list[ sample_definitions$condition ])
    ## The logic here is that I want by_type to be the default, but only
    ## if no one chooses either.
    if (!isTRUE(by_type) & !isTRUE(by_sample)) {
        by_type = TRUE
    }
    if (isTRUE(by_type)) {
        sample_definitions$counts = paste("processed_data/count_tables/", tolower(sample_definitions$Type),
            "/", tolower(sample_definitions$Stage), "/",
            sample_definitions$Sample.ID, suffix, sep="")
        sample_definitions$intercounts = paste("data/count_tables/", tolower(sample_definitions$Type),
            "/", tolower(sample_definitions$Stage), "/",
            sample_definitions$Sample.ID, "_inter", suffix, sep="")
    }
    if (isTRUE(by_sample)) {
        sample_definitions$counts = paste("processed_data/count_tables/", as.character(sample_definitions$Sample.ID), "/", as.character(sample_definitions$Sample.ID), suffix, sep="")
        sample_definitions$intercounts = paste("processed_data/count_tables/", as.character(sample_definitions$Sample.ID), "/", as.character(sample_definitions$Sample.ID), "_inter", suffix, sep="")
    }
    if (is.null(sample_definitions$stage)) {
        sample_definitions$stage = sample_definitions$Stage
    }
    if (is.null(sample_definitions$type)) {
        sample_definitions$type = sample_definitions$Type
    }
    if (is.null(sample_definitions$batch)) {
        sample_definitions$batch = sample_definitions$Batch
    }
    sample_definitions = as.data.frame(sample_definitions)
    rownames(sample_definitions) = sample_definitions$Sample.ID
    all_count_tables = my_read_files(as.character(sample_definitions$Sample.ID),
        as.character(sample_definitions$counts), header=header)
    all_count_matrix = as.matrix(all_count_tables)
    if (!is.null(genes)) {
        if (is.null(genes$ID)) {
            genes$ID = rownames(genes)
        }
        gene_info = genes[genes$ID %in% rownames(all_count_matrix),]
        all_count_matrix = all_count_matrix[rownames(all_count_matrix) %in% genes$ID,]                
    }
    metadata = new("AnnotatedDataFrame", data.frame(sample=as.character(sample_definitions$Sample.ID),
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
    return(experiment)
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
expt_subset = function(expt, subset) {
    if (class(expt) == "ExpressionSet") {
        expressionset = expt
    } else if (class(expt) == "expt") {
        expressionset = expt$expressionset
    } else {
        stop("expt is neither an expt nor ExpressionSet")
    }
    initial_metadata = Biobase::pData(expressionset)
    r_expression=paste("subset(initial_metadata,", subset, ")")
    samples = eval(parse(text=r_expression))
##    design = data.frame(sample=samples$sample, condition=samples$condition, batch=samples$batch)
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
#' 
#' @return  initial_count a data frame of count tables
#' @seealso \code{\link{create_experiment}}
#' 
#' @export
#' @examples
#' ## count_tables = my_read_files(as.character(sample_ids), as.character(count_filenames))
my_read_files = function(ids, files, header=FALSE) {
    initial_count = read.table(files[1])
    colnames(initial_count) = c("ID", ids[1])
    for (table in 2:length(files)) {
        tmp_count = read.table(files[table], header=header)
        colnames(tmp_count) = c("ID", ids[table])
        initial_count = merge(initial_count, tmp_count, by="ID")
    }
    rm(tmp_count)
    rownames(initial_count) = initial_count$ID
    initial_count = initial_count[-1]
    colnames(initial_count) = ids
    return(initial_count)
}
