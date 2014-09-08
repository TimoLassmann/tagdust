#!/bin/bash

echo "Running casava regression test(s):";

printf "     %s\n"  "1) Default paired-end read with index.";

error=$(${valparam} ../src/tagdust_rtest -seed 42 -arch ${devdir}/casava_arch.txt  ${devdir}/casava_read1.fastq.gz ${devdir}/casava_read2.fastq.gz ${devdir}/casava_read3.fastq.gz   -o casava_out 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  tagdust SUCCESS;
else
printf "%20s%10s\n"  tagdust FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

printf "     %s\n"  " Comparing to old PhiX results (1)";
if cmp -s <( sort ${devdir}/casava_out_BC_TTAGGC_READ1_gold.txt) <(sort  casava_out_BC_TTAGGC_READ1.fq); then
printf "%20s%10s\n"  results SUCCESS;
else
printf "%20s%10s\n"  results FAIL;
exit 1;
fi

printf "     %s\n"  " Comparing to old PhiX results (2)";
if cmp -s <( sort ${devdir}/casava_out_BC_TTAGGC_READ2_gold.txt) <(sort  casava_out_BC_TTAGGC_READ2.fq); then
printf "%20s%10s\n"  results SUCCESS;
else
printf "%20s%10s\n"  results FAIL;
exit 1;
fi

