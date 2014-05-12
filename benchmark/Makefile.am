AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

EXTRA_DIST = \
README.txt \
5barread3.sh \ 
barread.sh \
all_arch.txt

CLEANFILES = *.fq   *_logfile.txt *_tagdust_arch.txt *btrim*  *_fastxbarcodefile.txt *_fastx*

.PHONY: benchmark


benchmark:
	./barread.sh -b ../dev/EDITTAG_4nt_ed_2.txt
	./barread.sh -b ../dev/EDITTAG_6nt_ed_3.txt
	./5barread3.sh -b ../dev/EDITTAG_4nt_ed_2.txt
	./5barread3.sh -b ../dev/EDITTAG_6nt_ed_3.txt
	@echo HelloWorld



