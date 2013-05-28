AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}


bin_PROGRAMS = tagdust
tagdust_SOURCES = src/spline_fitting.c src/spline_fitting.h src/exact.h src/barcode_hmm.h  src/interface.h src/misc.h src/io.h src/nuc_code.h src/tagdust2.h src/main.c src/interface.c src/nuc_code.c src/misc.c src/io.c src/barcode_hmm.c  src/exact.c src/sim.c 
tagdust_LDFLAGS =  -lpthread -lm


dist_noinst_SCRIPTS = autogen.sh


