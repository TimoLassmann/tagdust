
EXTRA_DIST = $(PNG)  $(PDF) $(HTML) style.css template_modern.tex

if HAVE_PANDOC


dirs=$(top_srcdir)/reproducibility/scripts

MARKDOWNS=$(foreach dir,$(dirs), $(sort $(wildcard $(dir)/*.md)))

PDF=$(notdir $(patsubst %.md,%.pdf,  $(MARKDOWNS)))
HTML=$(notdir $(patsubst %.md,%.html,  $(MARKDOWNS)))
PNG=$(wildcard *.png)


vpath %.md  $(top_srcdir)/reproducibility/scripts
vpath %.pdf  $(top_srcdir)/reproducibility/doc/
vpath %.html $(top_srcdir)/reproducibility/doc/
vpath %.tex $(top_srcdir)/reproducibility/doc/
vpath %.css $(top_srcdir)/reproducibility/doc/

.PHONY : all 

all : $(PDF) $(HTML)

%.pdf : %.md  $(top_srcdir)/reproducibility/doc/template_modern.tex
	pandoc -f markdown  -s $< --template=$(top_srcdir)/reproducibility/doc/template_modern.tex --toc  -o $@ --chapters  --listings


%.html: %.md $(top_srcdir)/reproducibility/doc/style.css
	pandoc -f markdown -s -S --toc -c  $(top_srcdir)/reproducibility/doc/style.css  -t html5   $<  -o $@

 

all-local:
	@echo "Everything seemed to work ok."

endif

CLEANFILES=
DISTCLEANFILES=

