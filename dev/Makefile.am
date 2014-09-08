AUTOMAKE_OPTIONS = serial-tests
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

EXTRA_DIST = \
README.txt \
EDITTAG_3nt_ed_1.txt \
EDITTAG_3nt_ed_2.txt \
EDITTAG_4nt_ed_2.txt \
EDITTAG_4nt_ed_3.txt \
EDITTAG_6nt_ed_3.txt \
EDITTAG_6nt_ed_4.txt \
bar_read_test.sh \
sanity_test.sh \
casava_test.sh \
barread1_tagdust_results_gold.txt \
barread2_tagdust_results_gold.txt \
barread_paired_tagdust_results_gold.txt \
read_paired_tagdust_results_gold.txt \
casava_read1.fastq.gz \
casava_read2.fastq.gz \
casava_read3.fastq.gz \
casava_arch.txt \
casava_out_BC_TTAGGC_READ1_gold.txt \
casava_out_BC_TTAGGC_READ2_gold.txt


TESTS_ENVIRONMENT = \
	devdir='$(top_srcdir)'/dev; \
	export devdir; \
	valparam='$(VALGRIND)'; \
	export valparam; 
	
TESTS=sanity_test.sh bar_read_test.sh casava_test.sh
TEST_EXTENSIONS= .sh
SH_LOG_COMPILER= $(BASH)

CLEANFILES = *.fq  *_btrim_pattern.txt *_logfile.txt *_tagdust_arch.txt *_results.txt  *_fastxbarcodefile.txt










