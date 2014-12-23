# Makefile to use knitr for package vignettes
# put all PDF targets here, separated by spaces
PDFS=myr.pdf myr.html
all: $(PDFS)
clean:
	rm -rf *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ *.html
%.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$*.Rmd', output_format='pdf_document')"
%.html: %.Rmd
	Rscript -e "library('knitrBootstrap') ; rmarkdown::render('$*.Rmd', output_format='knitrBootstrap::bootstrap_document')"

