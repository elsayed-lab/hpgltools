## Make sure I didn't introduce any stupid syntax errors.
expect_error(expect_error(library(hpgltools)))
## Make sure its friends load.
expect_error(expect_error(autoloads_all()))
## Load pasilla
expect_error(expect_error(require.auto("pasilla")))

## Load the pasilla data set
datafile = system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
## Load the counts and drop super-low counts genes
counts = read.table(datafile, header=TRUE, row.names=1)
counts = counts[rowSums(counts) > ncol(counts),]
## Set up a quick design to be used by cbcbSEQ and hpgltools
design = data.frame(row.names=colnames(counts),
    condition=c("untreated","untreated","untreated",
        "untreated","treated","treated","treated"),
    libType=c("single-end","single-end","paired-end",
        "paired-end","single-end","paired-end","paired-end"))
metadata = design
colnames(metadata) = c("condition", "batch")
metadata$Sample.id = rownames(metadata)

## Make sure it is still possible to create an expt
expect_error(expect_error(pasilla_expt = create_expt(count_dataframe=counts, meta_dataframe=metadata)))

cbcb_data = counts
hpgl_data = exprs(pasilla_expt$expressionset)
expect_equal(cbcb_data, hpgl_data)
