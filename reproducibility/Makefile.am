AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

EXTRA_DIST = rnaseq_pipeline_4tagdust_paper-0.8.tar.gz 

SUBDIRS = builder scripts doc bin


all: 
if HAVE_PANDOC
if GOT_64BIT
	@echo "All good";
endif
endif

