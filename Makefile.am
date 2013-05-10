AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}


bin_PROGRAMS = tagdust
tagdust_SOURCES = src/barcode_hmm.h src/dbscan.h src/interface.h src/misc.h src/pst.h src/gmm.h src/io.h src/nuc_code.h src/tagdust2.h src/main.c src/interface.c src/nuc_code.c src/misc.c src/io.c src/barcode_hmm.c src/pst.c src/gmm.c src/dbscan.c 
tagdust_LDFLAGS =  -lpthread -lm


dist_noinst_SCRIPTS = autogen.sh


