VERSION=2021.03
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
	rm -rf R/.Rhistory vignettes/.Rhistory R/EuPathDB R/*.rda R/*.Rdata R/*.gb R/*.tab
	rm -rf tests/testthat/circos tests/testthat/EuPathDB tests/testthat/excel tests/testthat/excel_test \
		tests/testthat/preprocessing tests/testthat/test_gprofiler \
		tests/testthat/saved_plots tests/testthat/excel_test_sig \
		tests/testthat/kegg_pathways tests/testthat/pathview \
		tests/testthat/UP000* tests/testthat/topgo \
		tests/testthat/pathview_in tests/testthat/eupathdb \
		tests/testthat/BSgenome* tests/testthat/testing_write_expt \
		tests/testthat/a909_sig	tests/testthat/a909_tables \
		tests/testthat/mtb_rmats tests/testthat/mtb_suppa \
		tests/testthat/.Rhistory tests/testthat/parasite
	rm -f tests/testthat/*.pdf tests/testthat/*.png tests/testthat/*.xlsx tests/testthat/*.rda \
		tests/testthat/*.gff tests/testthat/*.gb tests/testthat/*.map tests/testthat/*.xml \
		tests/testthat/*.Rdata tests/testthat/*.json tests/testthat/*.tab tests/testthat/*kgml* \
		tests/testthat/*.csv tests/testthat/*.html

clean_vignette:
	rm -f vignettes/*.rda vignettes/*.map vignettes/*.Rdata inst/reference/reference.pdf

covr: install
	@echo "Invoking covr"
	R -e "x <- covr::package_coverage('.', quiet=FALSE); covr::report(x, file='hpgltools-report.html')"
##	R -e "x <- covr::package_coverage('.', type='all', quiet=FALSE); covr::report(x, file='hpgltools-report.html')"

deps:
	@echo "Invoking devtools::install_dev_deps()"
	R -e "all = as.data.frame(devtools::dev_package_deps('.', dependencies=TRUE)); needed = all[['diff']] < 0; needed = all[needed, 'package']; BiocManager::install(needed)"

document: roxygen vignette reference

install: reference roxygen clean
	echo "Performing R CMD INSTALL hpgltools."
	R CMD INSTALL --install-tests .

install_bioconductor:
	R -e "library(hpgltools); bioc_all()"

prereq:
	@echo "Checking a few prerequisites that seem to fall between the cracks sometimes."
	R -e "if (! 'BiocManager' %in% installed.packages()) { install.packages('BiocManager', repo='http://cran.rstudio.com/') }"
	R -e "bioc_prereq <- c('devtools', 'R.utils', 'pasilla','testthat','roxygen2','Biobase','preprocessCore','devtools','rmarkdown','knitr','ggplot2','data.table','foreach','survival');\
for (req in bioc_prereq) { if (! req %in% installed.packages()) { BiocManager::install(req) } }"

push:
	echo "Pushing to github."
	git commit -a && git push

reference:
	@echo "Generating reference manual with R CMD Rd2pdf"
	mkdir -p inst/doc
	rm -f inst/doc/reference.pdf
	R CMD Rd2pdf . -o inst/doc/reference.pdf --no-preview 2>/dev/null &

roxygen:
	@echo "Generating documentation with devtools::document()"
	R -e "suppressPackageStartupMessages(devtools::document())"

test: install
	@echo "Running run_tests.R"
	R -e "library(hpgltools); library(testthat); test_local(path = '.', reporter = 'summary', stop_on_failure = FALSE)"

vignette:
	@mkdir -p doc
	@echo "Building vignettes with devtools::build_vignettes()"
	R -e "devtools::build_vignettes(install=FALSE)"
	mv doc/* inst/doc/
	cp inst/reference/* inst/doc
	cp vignettes/*.R inst/doc
	cp vignettes/*.Rmd inst/doc
	cp vignettes/*.html inst/doc

vt:	clean_vignette vignette reference install
