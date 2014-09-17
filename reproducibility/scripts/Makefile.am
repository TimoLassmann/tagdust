
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}
EXTRA_DIST = run_benchmark.md

PIPEFILES= run.mk plotting.R 5barread3.sh barread.sh Figure1.tex


DISTCLEANFILES = all_arch.txt $(PIPEFILES)


CLEANFILES =   *.fq   *_logfile.txt *_tagdust_arch.txt *btrim*  *_fastxbarcodefile.txt *_fastx*

.PHONY: run


all: $(PIPEFILES)  run

run:
	$(MAKE) -f run.mk message


$(PIPEFILES): $(top_srcdir)/reproducibility/scripts/run_benchmark.md $(abs_top_builddir)/reproducibility/builder/mdlp
	$(abs_top_builddir)/reproducibility/builder/mdlp   $(top_srcdir)/reproducibility/scripts/run_benchmark.md



all-local:
	@echo "\
	\
	\
	Everything seemed to work ok.\
	\
	\
	"

