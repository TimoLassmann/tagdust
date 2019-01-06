EXTRA_DIST = casava_6nt_arch.txt UsingTagDustwithCASA.pdf UsingTagDustwithCASA.tex



LATEX=pdflatex

all:
	@echo Everything seemed to work ok.


doc: UsingTagDustwithCASA.pdf


UsingTagDustwithCASA.pdf: UsingTagDustwithCASA.tex
	$(LATEX) $<
	$(LATEX) $<


CLEANFILES= *.bcf *.nav *.snm *.log *.aux *.dvi *.bbl *.blg *.ilg *.toc *.lof *.lot *.idx *.ind *.ps  *~ *.out  *.run.xml

DISTCLEANFILES=



