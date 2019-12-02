ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS =  @SIMD_FLAGS@ @OPENMP_CFLAGS@

CLEANFILES =

AM_LDFLAGS = @DEAD_STRIP@
AM_LDFLAGS += -Wno-undef



bin_PROGRAMS = tagdust

tagdust_SOURCES = \
tldevel.h \
tldevel.c \
main.c \
interface.h \
interface.c \
arch_lib.h \
arch_lib.c \
seq_stats.h \
seq_stats.c \
io.h \
io.c \
misc.h \
misc.c \
nuc_code.h \
nuc_code.c \
hmm.h \
init_hmm.h \
init_hmm.c \
hmm_model_bag.h \
hmm_model_bag.c \
core_hmm_functions.h \
core_hmm_functions.c \
test_arch.h \
test_arch.c \
calibrate_hmm.h \
calibrate_hmm.c \
rng.h \
rng.c

tagdust_LDADD = -lpthread -lm
# misc.h \
# misc.c \
# nuc_code.h \
# nuc_code.c \
# barcode_hmm.h \
# barcode_hmm.c \
# seq_stats.h \
# seq_stats.c \
# hmm_model_bag.h \
# hmm_model_bag.c \
# core_hmm_functions.h \
# core_hmm_functions.c \
# test_architectures.c \
# io.c


# tagdust_CPPFLAGS = -DTAGDUST

# tagdust_SOURCES =  \
# kslib.c \
# main.c \
# interface.c \
# nuc_code.c \
# misc.c \
# io.c \
# barcode_hmm.c \
# calibrateQ.c \
# test_architectures.c

# tagdust_LDADD = -lpthread -lm
# tagdust_CPPFLAGS = -DTAGDUST


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

TESTS =

TESTS_ENVIRONMENT = \
	devdir='$(top_srcdir)'/dev; \
         export devdir; \
	$(VALGRIND)

check_PROGRAMS = arch_lib_test interface_test


arch_lib_test_SOURCES = \
arch_lib.h \
arch_lib.c \
tldevel.h \
tldevel.c 

arch_lib_test_CPPFLAGS = $(AM_CPPFLAGS) -DARCH_TEST

interface_test_SOURCES = \
interface.h \
interface.c \
tldevel.h \
tldevel.c

interface_test_CPPFLAGS = $(AM_CPPFLAGS) -DINTERFACE_TEST



# tagdustiotest_SOURCES = io.c misc.c interface.c kslib.c
# tagdustiotest_LDADD = $(tagdust_LDADD)
# tagdustiotest_CPPFLAGS  =  -DUTEST -DTAGDUST

# tagdust_rtest_SOURCES = $(tagdust_SOURCES)
# tagdust_rtest_LDADD = $(tagdust_LDADD)
# tagdust_rtest_CPPFLAGS  =  -DRTEST -DTAGDUST

# simreads_rtest_SOURCES = $(simreads_SOURCES)
# simreads_rtest_LDADD = $(simreads_LDADD)
# simreads_rtest_CPPFLAGS  =  -DRTEST -DSIMREADS

# evalres_rtest_SOURCES = $(evalres_SOURCES)
# evalres_rtest_LDADD = $(evalres_LDADD)
# evalres_rtest_CPPFLAGS  =  -DRTEST -DEVALRES






