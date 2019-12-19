
SUBDIRS = dev

ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST= README.org autogen.sh version.sh VERSION

AM_LDFLAGS = @DEAD_STRIP@
AM_LDFLAGS += -Wno-undef


noinst_LIBRARIES = libtldevel.a

TLHASHTABLESOURCES = tlhashtable.h tlhashtable.c
TLMINHASHSOURCES = tlminhash.h tlminhash.c
TLDEVELSOURCES = tldevel.h  tldevel.c
TLRNGSOURCES = tlrng.c tlrng.h
TLLOGSUMSOURCES = tllogsum.h tllogsum.c
TLBITVECSOURCES = tlbitvec.h tlbitvec.c
TLRBTREESOURCES  = tlrbtree.h tlrbtree.c
TLCHECKPOINTSOURCES = tlcheckpoint.h tlcheckpoint.c
TLMISCSOURCES = tlmisc.h tlmisc.c
TLALPHABETSOURCES = tlalphabet.h  tlalphabet.c
TLSEQIOSOURCES = tlseqio.h tlseqio.c

libtldevel_a_SOURCES = \
$(TLHASHTABLESOURCES) \
$(TLMINHASHSOURCES) \
$(TLDEVELSOURCES) \
$(TLRNGSOURCES) \
$(TLLOGSUMSOURCES) \
$(TLBITVECSOURCES) \
$(TLRBTREESOURCES) \
$(TLCHECKPOINTSOURCES) \
$(TLMISCSOURCES) \
$(TLALPHABETSOURCES)
# !ENABLE_ZLIB actually works!
if ENABLE_ZLIB
libtldevel_a_SOURCES += $(TLSEQIOSOURCES)
endif
