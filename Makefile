VERSION=2016.02
install: clean document
	cd ../ && R CMD INSTALL hpgltools

all: clean test document reference check build install clean

reference:
	rm -f inst/doc/reference.pdf
	R CMD Rd2pdf . -o inst/doc/reference.pdf

check:
	cd ../ && R CMD check hpgltools --no-build-vignettes

build:
	cd ../ && R CMD BUILD hpgltools

test:
	./run_tests.R

roxygen:
	rm -f NAMESPACE && Rscript -e "roxygen2::roxygenize()"

document:
	rm -f NAMESPACE && Rscript -e "devtools::document()"

vignette:
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
seqt = try(lirary('SeqTools'));\
if (class(SeqTools)) == 'try-error') { install_github('lianos/seqtools/R/pkg'); library('SeqTools') };\
prep = try(library('preprocessCore'));\
if (class(prep) == 'try-error') { biocLite('preprocessCore'); library('preprocessCore') };\
devtools = try(library('devtools'));\
if (class(devtools) == 'try-error') { biocLite('devtools'); library('devtools') };\
rmarkdown = try(library('rmarkdown')); \
if (class(rmarkdown) == 'try-error') {install_github('rstudio/rmarkdown'); library('rmarkdown') };\
knitrbootstrap = try(library('knitrBootstrap'));\
if (class(knitrbootstrap) == 'try-error') { install_github('jimhester/knitrBootstrap'); library('knitrBootstrap') };\
cbcbSEQ = try(library('cbcbSEQ')); \
if (class(cbcbSEQ) == 'try-error') { install_github('kokrah/cbcbSEQ'); library('cbcbSEQ') };\
" ;	cd .. && R CMD INSTALL hpgltools_0.1.tar.gz ; cd hpgltools ; Rscript -e "library(hpgltools); autoloads_all(update=TRUE)"
