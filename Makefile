VERSION=2017.10
export _R_CHECK_FORCE_SUGGESTS_=FALSE

all: clean roxygen reference check build test

build:
	@echo "Performing build with R CMD build hpgltools"
	R CMD build .

check:
	@echo "Performing check with R CMD check hpgltools"
	export _R_CHECK_FORCE_SUGGESTS_=FALSE && R CMD check . --no-build-vignettes

clean:
	@echo "Cleaning up"
	rm -rf hpgltools
	rm -rf ./..Rcheck &
	rm -rf hpgltools.Rcheck/
	rm -f hpgltools_${VERSION}.tar.gz
	rm -f inst/*.fai
	rm -rf vignettes/circos vignettes/pasilla_* vignettes/org.Spombe.eg.db \
      vignettes/wt_mga  vignettes/wt_mga_sig
	rm -f vignettes/*.gff vignettes/*.pdf vignettes/gene2pubmed.gz vignettes/NCBI.sqlite
	rm -rf R/.Rhistory vignettes/.Rhistory R/eupathdb R/*.rda R/*.Rdata
	for testdir in travis all_functions slow_tests; do \
	  rm -rf tests/$${testdir}/circos tests/$${testdir}/excel tests/$${testdir}/excel_test tests/$${testdir}/preprocessing tests/$${testdir}/test_gprofiler \
	    tests/$${testdir}/saved_plots tests/$${testdir}/excel_test_sig tests/$${testdir}/kegg_pathways tests/$${testdir}/pathview tests/$${testdir}/UP000* tests/$${testdir}/topgo \
	    tests/$${testdir}/pathview_in tests/$${testdir}/eupathdb tests/$${testdir}/BSgenome* tests/$${testdir}/testing_write_expt ;\
	  rm -f tests/$${testdir}/*.pdf tests/$${testdir}/*.png tests/$${testdir}/*.xlsx tests/$${testdir}/*.rda \
	    tests/$${testdir}/*.gff tests/$${testdir}/*.gb tests/$${testdir}/*.map tests/$${testdir}/*.xml \
	    tests/$${testdir}/*.Rdata tests/$${testdir}/*.json tests/$${testdir}/*.tab ;\
	done

clean_vignette:
	rm -f vignettes/*.rda vignettes/*.map vignettes/*.Rdata

deps:
	@echo "Invoking devtools::install_dev_deps()"
	R -e "suppressPackageStartupMessages(suppressMessages(source('http://bioconductor.org/biocLite.R')));\
all = as.data.frame(devtools::dev_package_deps('.', dependencies=TRUE)); needed = all[['diff']] < 0; needed = all[needed, 'package']; for (t in needed) { biocLite(t) }"

document: roxygen vignette reference

install:
	@echo "Performing R CMD INSTALL hpgltools globally."
	R CMD INSTALL .

install_bioconductor:
	R -e "library(hpgltools); bioc_all()"

prereq:
	@echo "Checking a few prerequisites."
	R -e "suppressPackageStartupMessages(suppressMessages(source('http://bioconductor.org/biocLite.R')));\
bioc_prereq <- c('R.utils', 'pasilla','testthat','roxygen2','Biobase','preprocessCore','devtools','rmarkdown','knitr','ggplot2','data.table','foreach','survival');\
for (req in bioc_prereq) { if (class(try(suppressMessages(eval(parse(text=paste0('library(', req, ')')))))) == 'try-error') { biocLite(req) } } \
## hahaha looks like lisp!"

push:
	echo "Pushing to github."
	git commit -a && git push

reference:
	@echo "Generating reference manual with R CMD Rd2pdf"
	rm -f inst/doc/reference.pdf
	R CMD Rd2pdf . -o inst/doc/reference.pdf --no-preview

roxygen:
	@echo "Generating documentation with devtools::document()"
	R -e "suppressPackageStartupMessages(devtools::document())"

suggests:
	@echo "Installing suggested packages."
	R -e "source('http://bioconductor.org/biocLite.R');\
library(desc);\
d = description\$$new(); suggests = d\$$get('Suggests');\
 suggests = gsub(pattern='\\n', replacement='', x=suggests);\
 suggests = gsub(pattern=' ', replacement='', x=suggests);\
 suggests = strsplit(x=suggests, split=',');\
 for (pkg in suggests[[1]]) { if (! pkg %in% installed.packages()) { biocLite(pkg); } else { message(paste0(pkg, ' is already installed.')) } };"

test: roxygen
	@echo "Installing hpgltools."
	R CMD INSTALL .
	@echo "Running run_tests.R"
	tests/testthat.R

update:
	R -e "source('http://bioconductor.org/biocLite.R'); biocLite(); library(BiocInstaller); biocValid()"

update_bioc:
	R -e "source('http://bioconductor.org/biocLite.R'); biocLite(); biocLite('BiocUpgrade');"

vignette:
	@echo "Building vignettes with devtools::build_vignettes()"
	R -e "devtools::build_vignettes()"

vt:	clean_vignette vignette reference install
