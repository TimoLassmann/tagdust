
EXTRA_DIST = $(PNG)  $(PDF) $(HTML) style.css template_modern.tex


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



%.pdf : %.md  $(top_srcdir)/reproducibility/doc/template_modern.tex
	pandoc -f markdown  -s $< --template=$(top_srcdir)/reproducibility/doc/template_modern.tex --toc -o $@ --chapters  --listings


%.html: %.md $(top_srcdir)/reproducibility/doc/style.css
	pandoc -f markdown -s -S --toc -c  $(top_srcdir)/reproducibility/doc/style.css  -t html5   $<  -o $@

 

all-local:
	@echo "Everything seemed to work ok."

if HAVE_PANDOC

all : $(PDF) $(HTML)

else

all : 
	@echo "Skipping pandoc on benchmark markdown file.";
endif




CLEANFILES=
DISTCLEANFILES=

