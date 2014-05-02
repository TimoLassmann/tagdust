
bin_PROGRAMS = tagdust simreads evalres


evalres_SOURCES = \
malloc_macro.h \
evalresults.c \
interface.h \
interface.c \
misc.h \
misc.c \
nuc_code.h \
nuc_code.c \
io.h \
io.c

evalres_LDADD = -lm

simreads_SOURCES= \
malloc_macro.h \
tagdust2.h \
interface.h \
interface.c \
nuc_code.h \
nuc_code.c \
misc.h \
misc.c \
simulate_reads.c
simreads_CPPFLAGS = -DSIMREADS
simreads_LDADD = -lm

tagdust_SOURCES =  \
malloc_macro.h \
barcode_hmm.h \
interface.h \
misc.h \
io.h \
nuc_code.h \
tagdust2.h \
main.c \
interface.c \
nuc_code.c \
misc.c \
io.c \
barcode_hmm.c \
paired_pre_processing.c \
calibrateQ.c \
test_architectures.c

tagdust_LDADD = -lpthread -lm


TESTS_ENVIRONMENT = \
	devdir='$(top_srcdir)'/dev; \
         export devdir;
	$(VALGRIND)

check_PROGRAMS = tagdust_rtest simreads_rtest evalres_rtest

tagdust_rtest_SOURCES = $(tagdust_SOURCES)
tagdust_rtest_LDADD = $(tagdust_LDADD)
tagdust_rtest_CPPFLAGS  =  -DRTEST

simreads_rtest_SOURCES = $(simreads_SOURCES)
simreads_rtest_LDADD = $(simreads_LDADD)
simreads_rtest_CPPFLAGS  =  -DRTEST

evalres_rtest_SOURCES = $(evalres_SOURCES)
evalres_rtest_LDADD = $(evalres_LDADD)
evalres_rtest_CPPFLAGS  =  -DRTEST






