




ACLOCAL_AMFLAGS = -I m4

SUBDIRS = $(TLDEVELDIR)


SUBDIRS +=  src #dev  casava_demo doc
dist_noinst_SCRIPTS = autogen.sh


count:
	line=$$(find . \( -name "*.[c,h]" -o -name "*.am" -o -name "*.org" -o -name "*.sh" \)  -type f -print0  | xargs -0 wc -l | grep total); echo	$$line
