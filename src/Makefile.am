
bin_PROGRAMS = tagdust simreads evalres


evalres_SOURCES = \
evalresults.c \
interface.h \
interface.c \
misc.h \
misc.c \
io.h \
io.c

evalres_LDADD = -lm

simreads_SOURCES= \
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
exact.h \
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
exact.c \
paired_pre_processing.c \
calibrateQ.c \
test_architectures.c

tagdust_LDADD = -lpthread -lm


TESTS_ENVIRONMENT = \
	devdir='$(top_srcdir)'/dev; \
         export devdir;
	$(VALGRIND)
	
#TESTS = tagdust_test

#check_PROGRAMS = simulate_reads run_programs eval_results

#tagdust_test_SOURCES = sserdt.c sserdt.h
#tagdust_test_CPPFLAGS  =  -DITEST

	





