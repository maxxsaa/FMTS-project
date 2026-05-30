.PHONY: all figures pdf clean

RSCRIPT ?= Rscript
LATEX   ?= pdflatex

all: pdf

figures:
	cd .. && $(RSCRIPT) run_analysis.R

pdf: figures
	$(LATEX) -interaction=nonstopmode report.tex
	$(LATEX) -interaction=nonstopmode report.tex

clean:
	rm -f report.aux report.log report.out report.toc
