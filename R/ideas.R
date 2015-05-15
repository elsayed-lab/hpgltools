## Time-stamp: <Thu May 14 14:41:15 2015 Ashton Trey Belew (abelew@gmail.com)>

## Idea from Ginger/Hector
## Take the qvalues from two datasets
## Rank order them and calculate the rank order correlation of those
## lists.  Similarly, rank order the logFC and correlate

## Using log-odds plots of fdr between two sets of genes in order to
## answer the question:

## How similar are these two sets?
## As the fdr moves towards 1.0, the logodds ratio should drop to 0.
## The more before 1.0, the 'better'

## It is more difficult to discern important genes in the low power
## experiment with respect to the high power than viceversa.
## The same set of tools will not likely work for that.

## sylvio scatterplot of fold change 24 v 72 hours
## look for oddities in the shape

## scatterplot of t/b from limma's toptable
## histograms of t/b at 24/72

## Following makeContrasts(), most easily following write_limma()
## since it creates a list of all the comparisons
## Have functions which pull the pairwise B/t/fold change numbers and
## easily provide scatterplots/histograms of them to show changes in
## the distribution.


