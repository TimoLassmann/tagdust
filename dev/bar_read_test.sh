#!/bin/bash




echo "Running tagdust regression tests:";

printf "     %s\n"  "1) single-end barcode and read";


#make --silent clean

error=$(${valparam} ../src/simreads_rtest ${devdir}/EDITTAG_6nt_ed_4.txt -seed 42 -sim_barnum 4  -sim_readlen 20 -sim_readlen_mod 0 -sim_numseq 10000 -sim_endloss 0 -sim_random_frac 0.1 -o barread1.fq   -sim_error_rate 0.02  2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  sim_reads SUCCESS;
else
printf "%20s%10s\n"  sim_reads FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi



error=$(${valparam} ../src/tagdust_rtest -seed 42 barread1.fq -arch barread1.fq_tagdust_arch.txt -o barread1_tagdust 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  tagdust SUCCESS;
else
printf "%20s%10s\n"  tagdust FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi


error=$(${valparam} ../src/evalres_rtest -name tagdust  barread1_tagdust*.fq -o barread1_tagdust  2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  eval_run SUCCESS;
else
printf "%20s%10s\n"  eval_run FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

if cmp -s <( sort ${devdir}/barread1_tagdust_results_gold.txt ) <(sort  barread1_tagdust_results.txt); then
printf "%20s%10s\n"  results SUCCESS;
else
printf "%20s%10s\n"  results FAIL;
exit 1;
fi

#make --silent clean

printf "     %s\n"  "2) single-end, 5' and 3' linkers, barcode and read";

error=$(${valparam} ../src/simreads_rtest ${devdir}/EDITTAG_6nt_ed_4.txt -seed 42 -sim_barnum 4  -sim_5seq GGGGGGG -sim_3seq TTTTTTT  -sim_readlen 20 -sim_readlen_mod 0 -sim_numseq 10000 -sim_endloss 0 -sim_random_frac 0.1 -o barread2.fq   -sim_error_rate 0.02 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  sim_reads SUCCESS;
else
printf "%20s%10s\n"  sim_reads FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi



error=$(${valparam} ../src/tagdust_rtest  -seed 42   barread2.fq -arch barread2.fq_tagdust_arch.txt -o barread2_tagdust 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  tagdust SUCCESS;
else
printf "%20s%10s\n"  tagdust FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi


error=$(${valparam} ../src/evalres_rtest -name tagdust  barread2_tagdust*.fq -o barread2_tagdust  2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  eval_run SUCCESS;
else
printf "%20s%10s\n"  eval_run FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

if cmp -s <( sort ${devdir}/barread2_tagdust_results_gold.txt ) <(sort barread2_tagdust_results.txt); then
printf "%20s%10s\n"  results SUCCESS;
else
printf "%20s%10s\n"  results FAIL;
exit 1;
fi


#make --silent clean


printf "     %s\n"  "3) paired-end";
printf "     %s\n"  "   read1: 5' and 3' linkers and read";
printf "     %s\n"  "   read2: just read";

error=$(${valparam} ../src/simreads_rtest ${devdir}/EDITTAG_6nt_ed_4.txt -seed 42 -sim_barnum 0  -sim_5seq GGGGGGG -sim_3seq TTTTTTT  -sim_readlen 20 -sim_readlen_mod 0 -sim_numseq 10000 -sim_endloss 0 -sim_random_frac 0.1 -o barread5_read1.fq   -sim_error_rate 0.02 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  sim_reads SUCCESS;
else
printf "%20s%10s\n"  sim_reads FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

error=$(${valparam} ../src/simreads_rtest ${devdir}/EDITTAG_6nt_ed_4.txt -seed 42 -sim_barnum 0   -sim_readlen 20 -sim_readlen_mod 0 -sim_numseq 10000 -sim_endloss 0 -sim_random_frac 0.00  -o barread6_read2.fq   -sim_error_rate 0.02 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  sim_reads SUCCESS;
else
printf "%20s%10s\n"  sim_reads FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

cat barread5_read1.fq_tagdust_arch.txt barread6_read2.fq_tagdust_arch.txt > barread_combo2.fq_tagdust_arch.txt


error=$(${valparam} ../src/tagdust_rtest  -seed 42  -sim_numseq 1 barread5_read1.fq  barread6_read2.fq  -arch barread_combo2.fq_tagdust_arch.txt  -o barread_paired2_tagdust 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  tagdust SUCCESS;
else
printf "%20s%10s\n"  tagdust FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi


error=$(${valparam} ../src/evalres_rtest -name tagdust barread_paired2_tagdust_*READ1.fq -o read_paired_tagdust  2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  eval_run SUCCESS;
else
printf "%20s%10s\n"  eval_run FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

if cmp -s <( sort ${devdir}/read_paired_tagdust_results_gold.txt ) <(sort  read_paired_tagdust_results.txt); then
printf "%20s%10s\n"  results SUCCESS;
else
printf "%20s%10s\n"  results FAIL;
exit 1;
fi


#make --silent clean

printf "     %s\n"  "4) paired-end";
printf "     %s\n"  "   read1: 5' and 3' linkers, barcode and read";
printf "     %s\n"  "   read2: just read";

error=$(${valparam} ../src/simreads_rtest ${devdir}/EDITTAG_6nt_ed_4.txt -seed 42 -sim_barnum 4  -sim_5seq GGGGGGG -sim_3seq TTTTTTT  -sim_readlen 20 -sim_readlen_mod 0 -sim_numseq 10000 -sim_endloss 0 -sim_random_frac 0.1 -o barread3_read1.fq   -sim_error_rate 0.02 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  sim_reads SUCCESS;
else
printf "%20s%10s\n"  sim_reads FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

error=$(${valparam} ../src/simreads_rtest ${devdir}/EDITTAG_6nt_ed_4.txt -seed 42 -sim_barnum 0   -sim_readlen 20 -sim_readlen_mod 0 -sim_numseq 10000 -sim_endloss 0 -sim_random_frac 0.00  -o barread4_read2.fq   -sim_error_rate 0.02 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  sim_reads SUCCESS;
else
printf "%20s%10s\n"  sim_reads FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

cat barread3_read1.fq_tagdust_arch.txt barread4_read2.fq_tagdust_arch.txt > barread_combo.fq_tagdust_arch.txt


error=$(${valparam} ../src/tagdust_rtest  -seed 42  -sim_numseq 1 barread3_read1.fq  barread4_read2.fq  -arch barread_combo.fq_tagdust_arch.txt  -o barread_paired_tagdust 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  tagdust SUCCESS;
else
printf "%20s%10s\n"  tagdust FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi


error=$(${valparam} ../src/evalres_rtest -name tagdust barread_paired_tagdust_*READ1.fq -o barread_paired_tagdust  2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  eval_run SUCCESS;
else
printf "%20s%10s\n"  eval_run FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

if cmp -s <( sort ${devdir}/barread_paired_tagdust_results_gold.txt ) <(sort  barread_paired_tagdust_results.txt); then
printf "%20s%10s\n"  results SUCCESS;
else
printf "%20s%10s\n"  results FAIL;
exit 1;
fi




make --silent clean





