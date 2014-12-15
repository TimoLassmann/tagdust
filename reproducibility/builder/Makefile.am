EXTRA_DIST = mdlp.md ostart.c

bin_PROGRAMS = mdlp
noinst_PROGRAMS = code_generator

# Source files for the mdlp
code_generator_SOURCES = ostart.c

# Source files for mdlp
nodist_mdlp_SOURCES = nstart.c

# generated.c is a built source that must be cleaned
CLEANFILES = nstart.c

# Build rule for generated.c
nstart.c:  code_generator$(EXEEXT)
	./code_generator$(EXEEXT) $(top_srcdir)/reproducibility/builder/mdlp.md







