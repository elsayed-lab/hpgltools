## The following was taken from ggplot2's ggplot2.r
## I presume it is a blanket importer cue for roxygen2 to add
## import statements to the NAMESPACE file so that when ggplot2 is used
## it will ensure that these libraries are available.
## I checked the roxygen documentation and it appears that
## imports are saved as the exclusive set, as a result repeating these
## at each function declaration serves to make explicit what each function
## requires while not (I think) adding excessive cruft to the NAMESPACE

## #' @import scales grid gtable
## #' @importFrom plyr defaults
## #' @importFrom stats setNames
## NULL
