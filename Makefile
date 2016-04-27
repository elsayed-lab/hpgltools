VERSION=2016.02
export _R_CHECK_FORCE_SUGGESTS_=FALSE

all: clean prereq document reference check build install test

install:
	echo "Performing R CMD INSTALL hpgltools"
	cd ../ && R CMD INSTALL hpgltools

reference:
	echo "Generating reference manual with R CMD Rd2pdf"
	rm -f inst/doc/reference.pdf
	R CMD Rd2pdf . -o inst/doc/reference.pdf

check:
	echo "Performing check with R CMD check hpgltools"
	cd ../ && export _R_CHECK_FORCE_SUGGESTS_=FALSE && R CMD check hpgltools --no-build-vignettes

build:
	echo "Performing build with R CMD build hpgltools"
	cd ../ && R CMD build hpgltools

test:
	echo "Running run_tests.R"
	./run_tests.R

roxygen:
	echo "Generating documentation with roxygen2::roxygenize()"
	rm -f NAMESPACE && Rscript -e "roxygen2::roxygenize()"

document:
	echo "Generating documentation with devtools::document()"
	rm -f NAMESPACE && Rscript -e "devtools::document()"

vignette:
	echo "Building vignettes with devtools::build_vignettes()"
	Rscript -e "devtools::build_vignettes()"

clean_vignette:
	rm -f inst/doc/*

vt:	clean_vignette vignette install

clean:
	rm -rf hpgltools/
	rm -rf hpgltools.Rcheck/
	rm -rf hpgltools_${VERSION}.tar.gz
	rm -rf $(find . -type f -name '*.rda' | grep -v 'hpgltools.rda')
	find . -type f -name '*.Rdata' -exec rm -rf {} ';' 2>/dev/null
	find . -type d -name excel -exec rm -rf {} ';' 2>/dev/null
	find . -type d -name reference -exec rm -rf {} ';' 2>/dev/null

autoloads:
	Rscript -e "library(devtools); devtools::load_all('.'); autoloads_all()"

prereq:
	Rscript -e "source('http://bioconductor.org/biocLite.R');\
pasilla = try(library('pasilla'));\
if (class(pasilla) == 'try-error') { biocLite('pasilla'); library('pasilla') };\
tt = try(library('testthat'));\
if (class(tt) == 'try-error') { biocLite('testthat'); library('testthat') };\
bb = try(library('Biobase'));\
if (class(bb) == 'try-error') { biocLite('Biobase'); library('Biobase') };\
prep = try(library('preprocessCore'));\
if (class(prep) == 'try-error') { biocLite('preprocessCore'); library('preprocessCore') };\
devtools = try(library('devtools'));\
if (class(devtools) == 'try-error') { biocLite('devtools'); library('devtools') };\
rmarkdown = try(library('rmarkdown')); \
if (class(rmarkdown) == 'try-error') {install_github('rstudio/rmarkdown'); library('rmarkdown') };\
knitrbootstrap = try(library('knitrBootstrap'));\
if (class(knitrbootstrap) == 'try-error') { install_github('jimhester/knitrBootstrap'); library('knitrBootstrap') };\
" ;
