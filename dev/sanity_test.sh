#!/bin/bash


echo "Running tagdust sanity tests:";

printf "     %s\n"  "1) Sanity 1: what happens with no architecture.";


#make --silent clean

error=$(${valparam} ../src/simreads_rtest ${devdir}/EDITTAG_6nt_ed_4.txt -seed 42 -sim_barnum 4  -sim_readlen 20 -sim_readlen_mod 0 -sim_numseq 10000 -sim_endloss 0 -sim_random_frac 0.1 -o sanity_barread1.fq   -sim_error_rate 0.02  2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  sim_reads SUCCESS;
else
printf "%20s%10s\n"  sim_reads FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

error=$(${valparam} ../src/tagdust_rtest -seed 42 sanity_barread1.fq -o sanity1 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  tagdust SUCCESS;
else
printf "%20s%10s\n"  tagdust FAILED;
echo "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi




make --silent clean




