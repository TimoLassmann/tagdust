
ACLOCAL_AMFLAGS = -I m4
AM_LDFLAGS = @DEAD_STRIP@
AM_LDFLAGS += -Wno-undef

LIBS = ../libtldevel.a -lm


AM_CPPFLAGS = -I$(top_srcdir)/lib_tldevel

TESTLIST = \
tldevel_test \
rng_test \
checkpoint_test \
rbtree_test \
logsum_test \
hash_test \
minhash_test

TESTS = $(TESTLIST)

check_PROGRAMS = $(TESTLIST)

tldevel_test_SOURCES = tldevel_test.c
tldevel_test_CPPFLAGS = $(AM_CPPFLAGS)
tldevel_test_LDADD = ../libtldevel.a

rng_test_SOURCES = tlrng_test.c
rng_test_CPPFLAGS = $(AM_CPPFLAGS)
rng_test_LDADD = ../libtldevel.a

checkpoint_test_SOURCES = tlcheckpoint_test.c
checkpoint_test_CPPFLAGS = $(AM_CPPFLAGS)
checkpoint_test_LDADD = ../libtldevel.a


rbtree_test_SOURCES = tlrbtree_test.c
rbtree_test_CPPFLAGS = $(AM_CPPFLAGS)
rbtree_test_LDADD = ../libtldevel.a

logsum_test_SOURCES = tllogsum_test.c
logsum_test_CPPFLAGS = $(AM_CPPFLAGS)
logsum_test_LDADD = ../libtldevel.a


hash_test_SOURCES = tlhashtable_test.c
hash_test_CPPFLAGS = $(AM_CPPFLAGS)
hash_test_LDADD = ../libtldevel.a

minhash_test_SOURCES = tlminhash_test.c
minhash_test_CPPFLAGS = $(AM_CPPFLAGS)
minhash_test_LDADD = ../libtldevel.a

CLEANFILES = \
	MAIN_0.chk


