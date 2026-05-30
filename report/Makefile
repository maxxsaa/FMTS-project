.PHONY: all figures tables pdf clean

RSCRIPT ?= Rscript
LATEX   ?= pdflatex

all: pdf

figures:
	cd .. && $(RSCRIPT) run_analysis.R

tables: figures
	cp ../outputs/tables/test_accuracy_by_method.csv tables/ 2>/dev/null || true
	cp ../outputs/tables/forecast_vs_test.csv tables/ 2>/dev/null || true
	cp ../outputs/tables/seasonally_adjusted.csv tables/ 2>/dev/null || true

pdf: figures tables
	$(LATEX) -interaction=nonstopmode report.tex
	$(LATEX) -interaction=nonstopmode report.tex

clean:
	rm -f report.aux report.log report.out report.toc
