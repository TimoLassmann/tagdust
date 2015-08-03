bin_PROGRAMS = tagdust simreads evalres merge rename_qiime

evalres_SOURCES = \
barcode_hmm.h \
io.h \
nuc_code.h \
interface.h \
kslib.h \
misc.h \
tagdust2.h \
kslib.c \
evalresults.c \
interface.c \
misc.c \
nuc_code.c \
io.c


evalres_LDADD = -lm
evalres_CPPFLAGS = -DEVALRES

simreads_SOURCES= \
nuc_code.h \
sim.h \
interface.h \
kslib.h \
misc.h \
tagdust2.h \
kslib.c \
interface.c \
nuc_code.c \
misc.c \
simulate_reads.c

simreads_CPPFLAGS = -DSIMREADS
simreads_LDADD = -lm

tagdust_SOURCES =  \
kslib.c \
main.c \
interface.c \
nuc_code.c \
misc.c \
io.c \
barcode_hmm.c \
calibrateQ.c \
test_architectures.c

tagdust_LDADD = -lpthread -lm
tagdust_CPPFLAGS = -DTAGDUST

merge_SOURCES = \
kslib.c \
interface.c \
nuc_code.c \
io.c \
misc.c \
merge.c

merge_LDADD = -lpthread -lm
merge_CPPFLAGS = -DMERGE


rename_qiime_SOURCES = \
rename_illumina_demultiplexed_for_qiime.c \
nuc_code.c \
kslib.c \
interface.c \
misc.c \
io.c

rename_qiime_CPPFLAGS = -DRENAME
rename_qiime_LDADD = -lm

TESTS = tagdustiotest

TESTS_ENVIRONMENT = \
	devdir='$(top_srcdir)'/dev; \
         export devdir; \
	$(VALGRIND)

check_PROGRAMS = tagdust_rtest simreads_rtest evalres_rtest tagdustiotest


tagdustiotest_SOURCES = io.c misc.c interface.c kslib.c
tagdustiotest_LDADD = $(tagdust_LDADD)
tagdustiotest_CPPFLAGS  =  -DUTEST -DTAGDUST

tagdust_rtest_SOURCES = $(tagdust_SOURCES)
tagdust_rtest_LDADD = $(tagdust_LDADD)
tagdust_rtest_CPPFLAGS  =  -DRTEST -DTAGDUST

simreads_rtest_SOURCES = $(simreads_SOURCES)
simreads_rtest_LDADD = $(simreads_LDADD)
simreads_rtest_CPPFLAGS  =  -DRTEST -DSIMREADS

evalres_rtest_SOURCES = $(evalres_SOURCES)
evalres_rtest_LDADD = $(evalres_LDADD)
evalres_rtest_CPPFLAGS  =  -DRTEST -DEVALRES






