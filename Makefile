VERSION=0.1
build:
	cd ../ && R CMD build myr --no-build-vignettes
	mv ../myr_${VERSION}.tar.gz .
	tar xavf myr_${VERSION}.tar.gz
	R CMD check myr --no-build-vignettes
	R CMD Rd2pdf myr && mv myr.pdf myr/inst/doc
	R CMD INSTALL myr

inst:
	cd ../ && R CMD build myr --no-build-vignettes && R CMD INSTALL myr && rm myr_${VERSION}.tar.gz

clean:
	rm -rf myr/
	rm -rf myr.Rcheck/
	rm -rf myr_${VERSION}.tar.gz

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
" ;	cd .. && R CMD build myr ; R CMD INSTALL myr_0.1.tar.gz ; cd myr && R CMD INSTALL googleVis && R CMD INSTALL BSgenome.Lmajor.friedlin ; Rscript -e "require(myr)" 
