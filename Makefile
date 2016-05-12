VERSION=2016.02
export _R_CHECK_FORCE_SUGGESTS_=FALSE

all: clean prereq document reference check build test

install:
	echo "Performing R CMD INSTALL hpgltools"
	cd ../ && R CMD INSTALL hpgltools && cd hpgltools

reference:
	echo "Generating reference manual with R CMD Rd2pdf"
	rm -f inst/doc/reference.pdf
	R CMD Rd2pdf . -o inst/doc/reference.pdf

check:
	echo "Performing check with R CMD check hpgltools"
	cd ../ && export _R_CHECK_FORCE_SUGGESTS_=FALSE && R CMD check hpgltools --no-build-vignettes && cd hpgltools

build:
	echo "Performing build with R CMD build hpgltools"
	cd ../ && R CMD build hpgltools && cd hpgltools

test: install
	echo "Running run_tests.R"
	./run_tests.R

roxygen:
	echo "Generating documentation with roxygen2::roxygenize()"
	Rscript -e "roxygen2::roxygenize()"

document:
	echo "Generating documentation with devtools::document()"
	Rscript -e "devtools::document()"

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
	find . -type f -name '*.Rdata' -exec rm -rf {} ';' 2>/dev/null
	find . -type d -name excel -exec rm -rf {} ';' 2>/dev/null
	find . -type d -name reference -exec rm -rf {} ';' 2>/dev/null

autoloads:
	Rscript -e "library('hpgltools'); autoloads_all();"

prereq:
	Rscript -e "suppressMessages(source('http://bioconductor.org/biocLite.R'));\
bioc_prereq <- c('pasilla','testthat','roxygen2','Biobase','preprocessCore','devtools','rmarkdown','knitr');\
for (req in bioc_prereq) { if (class(try(suppressMessages(eval(parse(text=paste0('library(', req, ')')))))) == 'try-error') { biocLite(req) } };\
## hahaha looks like lisp!"

update_bioc:
	Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite(); biocLite('BiocUpgrade');"

update:
	Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite(); library(BiocInstaller); biocValid()"
