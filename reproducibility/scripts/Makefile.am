
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}
EXTRA_DIST = run_benchmark.md

PIPEFILES= run.mk

CLEANFILES = $(PIPEFILES) 

.PHONY: run


all: $(PIPEFILES)  run

run:
	$(MAKE)	-k -f run.mk


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

