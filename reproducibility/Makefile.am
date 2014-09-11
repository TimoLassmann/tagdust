AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

if HAVE_PANDOC
DOC = doc
else
DOC =
endif


SUBDIRS = builder scripts $(DOC)


all: 
if HAVE_PANDOC
if GOT_64BIT
	@echo "${GOAHEAD}";
	@echo "GAGA";
endif
endif

