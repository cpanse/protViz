# Makefile
# Christian Panse <cp@fgcz.ethz.ch>
# 2014-11-08

MKDIR_P = mkdir -p

FILES = PTM_MarkerFinder.pdf
.PHONY: clean directories

LATEX   = pdflatex
BIBTEX  = bibtex
SVNFILES = PTM_MarkerFinder.Rnw
JSSFILES = $(SVNFILES) PTM_MarkerFinder.pdf


jss: sweave stangle
	pdflatex PTM_MarkerFinder.tex \
	&& bibtex PTM_MarkerFinder \
	&& pdflatex PTM_MarkerFinder.tex \
	&& pdflatex PTM_MarkerFinder.tex \
	&& ls -l $(JSSFILES)


sweave:
	R CMD Sweave PTM_MarkerFinder.Rnw

stangle:
	R CMD Stangle PTM_MarkerFinder.Rnw  && R --no-save < PTM_MarkerFinder.R

clean:
	@rm -fv *.aux *.log *.ps *.blg *.bbl *.txt *.out *.brf *.toc

