VERSION=2017.10
export _R_CHECK_FORCE_SUGGESTS_=FALSE

all: clean roxygen reference check build test

build:
	@echo "Performing build with R CMD build hpgltools"
	R CMD build .

check: roxygen
	@echo "Performing check with R CMD check hpgltools"
	rm -rf ./..Rcheck 2>/dev/null 1>&2
	export _R_CHECK_FORCE_SUGGESTS_=FALSE && R CMD check . --no-build-vignettes
	@rm -rf ./..Rcheck 2>/dev/null 1>&2

clean:
	@echo "Cleaning up"
	rm -rf hpgltools
	rm -rf ./..Rcheck &
	rm -rf hpgltools.Rcheck/
	rm -f hpgltools_${VERSION}.tar.gz
	rm -f inst/*.fai
	rm -rf vignettes/circos vignettes/pasilla_* vignettes/org.Spombe.eg.db \
      vignettes/wt_mga vignettes/wt_mga_sig vignettes/*_files
	rm -f vignettes/*.gff vignettes/*.pdf vignettes/gene2pubmed.gz vignettes/NCBI.sqlite \
			vignettes/*.R vignettes/*.html vignettes/*.gb vignettes/*.gz vignettes/*.Rdata \
			vignettes/*.xlsx vignettes/*.tex vignettes/*.log vignettes/*.aux vignettes/*.map \
			vignettes/*.rda
	rm -rf R/.Rhistory vignettes/.Rhistory R/EuPathDB R/*.rda R/*.Rdata
	rm -rf tests/testthat/circos tests/testthat/EuPathDB tests/testthat/excel tests/testthat/excel_test \
		tests/testthat/preprocessing tests/testthat/test_gprofiler \
		tests/testthat/saved_plots tests/testthat/excel_test_sig \
		tests/testthat/kegg_pathways tests/testthat/pathview \
		tests/testthat/UP000* tests/testthat/topgo \
 		tests/testthat/pathview_in tests/testthat/eupathdb \
		tests/testthat/BSgenome* tests/testthat/testing_write_expt ;\
	  rm -f tests/testthat/*.pdf tests/testthat/*.png tests/testthat/*.xlsx tests/testthat/*.rda \
		tests/testthat/*.gff tests/testthat/*.gb tests/testthat/*.map tests/testthat/*.xml \
		tests/testthat/*.Rdata tests/testthat/*.json tests/testthat/*.tab tests/testthat/*kgml* ;\

clean_vignette:
	rm -f vignettes/*.rda vignettes/*.map vignettes/*.Rdata inst/reference/reference.pdf

covr: install
	@echo "Invoking covr::codecov()"
	R -e "x <- covr::package_coverage('.'); covr::report(x)"

deps:
	@echo "Invoking devtools::install_dev_deps()"
	R -e "all = as.data.frame(devtools::dev_package_deps('.', dependencies=TRUE)); needed = all[['diff']] < 0; needed = all[needed, 'package']; for (t in needed) { BiocManager::install(t) }"

document: roxygen vignette reference

install: clean
	@echo "Performing R CMD INSTALL hpgltools."
	R CMD INSTALL --install-tests .

install_bioconductor:
	R -e "library(hpgltools); bioc_all()"

prereq:
	@echo "Checking a few prerequisites."
	R -e "install.packages('BiocManager', repo='http://cran.rstudio.com/')"
	R -e "bioc_prereq <- c('devtools', 'R.utils', 'pasilla','testthat','roxygen2','Biobase','preprocessCore','devtools','rmarkdown','knitr','ggplot2','data.table','foreach','survival');\
for (req in bioc_prereq) { if (class(try(suppressMessages(eval(parse(text=paste0('library(', req, ')')))))) == 'try-error') { BiocManager::install(req) } } \
## hahaha looks like lisp!"

push:
	echo "Pushing to github."
	git commit -a && git push

reference:
	@echo "Generating reference manual with R CMD Rd2pdf"
	@mkdir -p inst/doc
	R CMD Rd2pdf . -o inst/reference/reference.pdf --no-preview

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
	@mkdir -p doc
	@echo "Building vignettes with devtools::build_vignettes()"
	R -e "devtools::build_vignettes()"
	mv doc inst/doc
	cp inst/reference/* inst/doc

vt:	clean_vignette vignette reference install
