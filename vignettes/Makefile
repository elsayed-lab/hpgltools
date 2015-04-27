# Makefile to use knitr for package vignettes
# put all PDF targets here, separated by spaces
PDFS=hpgltools.pdf hpgltools.html fission_examples.pdf fission_examples.html examples.pdf examples.html

all: $(PDFS)
clean:
	rm -rf *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ *.html
%.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$*.Rmd', output_format='pdf_document')"
%.html: %.Rmd
	Rscript -e "library('knitrBootstrap') ; rmarkdown::render('$*.Rmd', output_format='knitrBootstrap::bootstrap_document')"

