.PHONY: all tables figures pdf clean

RSCRIPT ?= Rscript
LATEX   ?= pdflatex

all: pdf

figures:
	cd .. && $(RSCRIPT) timeseries_Rscript.R

tables:
	$(RSCRIPT) export_report_tables.R

pdf: tables
	$(LATEX) -interaction=nonstopmode report.tex
	$(LATEX) -interaction=nonstopmode report.tex

clean:
	rm -f report.aux report.log report.out report.toc report.pdf
