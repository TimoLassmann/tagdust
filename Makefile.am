
VFLAGS = --track-origins=yes --tool=memcheck --leak-check=yes --error-exitcode=1


AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

EXTRA_DIST = doc/User-Manual.pdf doc/latex/refman.pdf


SUBDIRS = src dev

dist_noinst_SCRIPTS = autogen.sh



rsync:
	ssh lm3 'rm -rf tmp/testing'
	rsync -ra -e ssh --exclude '/.git' --exclude '/bin' --exclude '/doc' --exclude '*.tar.gz' --exclude '/build' . lm3:tmp/testing

val: rsync
	ssh lm3 'cd tmp/testing && make clean && make check '




