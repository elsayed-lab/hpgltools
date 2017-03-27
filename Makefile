VERSION=2017.01
export _R_CHECK_FORCE_SUGGESTS_=FALSE

all: clean roxygen reference check build test

install: prereq roxygen
	@echo "Performing R CMD INSTALL hpgltools"
	@R CMD INSTALL .

reference:
	@echo "Generating reference manual with R CMD Rd2pdf"
	@rm -f inst/doc/reference.pdf
	@R CMD Rd2pdf . -o inst/doc/reference.pdf --no-preview

check:
	@echo "Performing check with R CMD check hpgltools"
	@export _R_CHECK_FORCE_SUGGESTS_=FALSE && R CMD check . --no-build-vignettes

build:
	@echo "Performing build with R CMD build hpgltools"
	@R CMD build .

test: install
	@echo "Running run_tests.R"
	tests/testthat.R

roxygen:
	@echo "Generating documentation with roxygen2::roxygenize()"
	@Rscript -e "suppressPackageStartupMessages(roxygen2::roxygenize())"

vignette:
	@echo "Building vignettes with devtools::build_vignettes()"
	@Rscript -e "devtools::build_vignettes()"

document: roxygen vignette reference

clean_vignette:
	@rm -f inst/doc/* vignettes/*.rda vignettes/*.map vignettes/*.Rdata

vt:	clean_vignette vignette reference install

clean:
	rm -rf hpgltools/
	rm -rf hpgltools.Rcheck/
	rm -rf hpgltools_${VERSION}.tar.gz
	find . -type f -name '*.Rdata' -exec rm -rf {} ';' 2>/dev/null
	find . -type d -name excel -exec rm -rf {} ';' 2>/dev/null
	find . -type d -name reference -exec rm -rf {} ';' 2>/dev/null

prereq:
	@Rscript -e "suppressPackageStartupMessages(suppressMessages(source('http://bioconductor.org/biocLite.R')));\
bioc_prereq <- c('pasilla','testthat','roxygen2','Biobase','preprocessCore','devtools','rmarkdown','knitr','ggplot2','data.table','foreach','survival');\
for (req in bioc_prereq) { if (class(try(suppressMessages(eval(parse(text=paste0('library(', req, ')')))))) == 'try-error') { biocLite(req) } }; install.packages(\".\", repos=NULL, type=\"source\"); message(\"Checking Dependencies.\"); devtools::install_deps(\".\", dependencies=TRUE); message(\"Finished dependencies.\"); \
## hahaha looks like lisp!"

update_bioc:
	Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite(); biocLite('BiocUpgrade');"

update:
	Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite(); library(BiocInstaller); biocValid()"

install_bioconductor:
	Rscript -e "library(hpgltools); bioc_all()"
