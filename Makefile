VERSION=0.1
build:
	cd ../ && R CMD build hpgltools --no-build-vignettes
	mv ../hpgltools_${VERSION}.tar.gz .
	tar xavf hpgltools_${VERSION}.tar.gz
	R CMD check hpgltools --no-build-vignettes
	R CMD Rd2pdf hpgltools && mv hpgltools.pdf hpgltools/inst/doc
	R CMD INSTALL hpgltools

inst:
	cd ../ && R CMD build hpgltools --no-build-vignettes && R CMD INSTALL hpgltools && rm hpgltools_${VERSION}.tar.gz

clean:
	rm -rf hpgltools/
	rm -rf hpgltools.Rcheck/
	rm -rf hpgltools_${VERSION}.tar.gz

install:
	Rscript -e "source('http://bioconductor.org/biocLite.R');\
pasilla = try(library('pasilla'));\
if (class(pasilla) == 'try-error') { biocLite('pasilla'); library('pasilla') };\
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
" ;	cd .. && R CMD build hpgltools ; R CMD INSTALL hpgltools_0.1.tar.gz ; cd hpgltools && R CMD INSTALL googleVis && R CMD INSTALL BSgenome.Lmajor.friedlin ; Rscript -e "require(hpgltools)" 
