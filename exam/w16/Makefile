LATEX = latex
PDFLATEX = pdflatex
BIBTEX = bibtex
RM = rm -f
CP = cp
TAR = tar
TARFLAGS = zcvf
PSTOPDF = ps2pdf
EPSTOPDF = epstopdf

REXE = R --vanilla
RCMD = $(REXE) CMD
RSCRIPT = Rscript --vanilla

default: pdf-docs

r-scripts: mt531w16.R

html-docs: mt531w16.pdf

pdf-docs: mt531w16.pdf

%.html: %.Rmd
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\",output_format=\"html_document\")"


%.tex: %.Rnw
	$(RSCRIPT) -e "library(knitr); knit(\"$*.Rnw\")"

%.pdf: %.tex
	$(PDFLATEX) $*
	-$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*

%.R: %.Rmd
	Rscript --vanilla -e "library(knitr); purl(\"$*.Rmd\",output=\"$*.R\")"

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	$(RM) -r cache figure

