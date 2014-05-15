AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

EXTRA_DIST = \
README.txt \
5barread3.sh \ 
barread.sh \
plotting.R

CLEANFILES = *.fq   *_logfile.txt *_tagdust_arch.txt *btrim*  *_fastxbarcodefile.txt *_fastx*

DISTCLEANFILES = all_arch.txt

.PHONY: benchmark

benchmark: rplot
	R --slave --vanilla < plotting.R
	@echo HelloWorld


rplot: barread4r.tsv 5barread34r.tsv
	@echo HelloWorld

barread4r.tsv: barread_results.txt
	cat barread_results.txt | awk 'BEGIN{print "Program\tsimerror\tbarcodes\tSensitivity\tSpecificity\tPrecision\tKappa"} {if(NR>1){x = split($$1,a,"_");printf "%s\t%f\t%f\t%f\t%f\t%f\t%f\n" ,  a[1],a[2],a[3],$$2,$$3,$$4,$$5 }}' > barread4r.tsv
	
5barread34r.tsv: 5barread3_results.txt
	cat 5barread3_results.txt | awk 'BEGIN{print "Program\tsimerror\tbarcodes\tSensitivity\tSpecificity\tPrecision\tKappa"} {if(NR>1){x = split($$1,a,"_");printf "%s\t%f\t%f\t%f\t%f\t%f\t%f\n" ,  a[1],a[2],a[3],$$2,$$3,$$4,$$5 }}' > 5barread34r.tsv


barread_results.txt: all_arch.txt
	./barread.sh -b ../dev/EDITTAG_4nt_ed_2.txt -r
	./barread.sh -b ../dev/EDITTAG_6nt_ed_3.txt -r
	
5barread3_results.txt: all_arch.txt
	./5barread3.sh -b ../dev/EDITTAG_4nt_ed_2.txt -r
	./5barread3.sh -b ../dev/EDITTAG_6nt_ed_3.txt -r

all_arch.txt: 
	./barread.sh -b ../dev/EDITTAG_4nt_ed_2.txt
	./barread.sh -b ../dev/EDITTAG_6nt_ed_3.txt
	./5barread3.sh -b ../dev/EDITTAG_4nt_ed_2.txt
	./5barread3.sh -b ../dev/EDITTAG_6nt_ed_3.txt
	cat  *tagdust_arch.txt  | sort | uniq  > all_arch.txt

