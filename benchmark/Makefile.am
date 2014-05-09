AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

EXTRA_DIST = \
README.txt \
5barread3.sh \ 
barread.sh

CLEANFILES = *.fq   *_logfile.txt *_tagdust_arch.txt *btrim*  *_fastxbarcodefile.txt *_fastx*

.PHONY: benchmark


benchmark:
	./barread.sh
	./5barread3.sh
	@echo HelloWorld



