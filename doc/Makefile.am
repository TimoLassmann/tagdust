EXTRA_DIST = User-Manual.tex\
User-Manual.pdf \
ref.bib \
figures/figure1.pdf



LATEX=pdflatex

all:
	@echo Everything seemed to work ok.



doc: User-Manual.pdf


User-Manual.pdf: User-Manual.tex User-Manual.bbl
	$(LATEX) $<


User-Manual.bbl: User-Manual.aux ref.bib
	bibtex $<


User-Manual.aux: User-Manual.tex
	$(LATEX) $?




CLEANFILES= *.bcf *.nav *.snm *.log *.aux *.dvi *.bbl *.blg *.ilg *.toc *.lof *.lot *.idx *.ind *.ps  *~ *.out  *.run.xml

DISTCLEANFILES= *pdf






