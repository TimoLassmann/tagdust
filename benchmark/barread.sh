#!/bin/bash




echo "Running tagdust benchmarks:";

array=(0.01 0.015 0.02 0.025 0.03)
len=${#array[*]}

numbar=(8 24 48)
lenbar=${#numbar[*]}

barcodes="EDITTAG_3nt_ed_1.txt"

make --silent clean

for (( c=0; c < $lenbar; c+=1 )); do

j=0

while [ $j -lt $len ]
do
make --silent clean
name="barread${array[$j]}.fq"

error=$(../src/simreads ../dev/$barcodes  -seed 42 -sim_barnum ${numbar[$c]}  -sim_readlen 20 -sim_readlen_mod 0 -sim_numseq 100000 -sim_endloss 0 -sim_random_frac 0.1 -o $name   -sim_error_rate ${array[$j]}  2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  sim SUCCESS;
else
printf "%20s%10s\n"  sim FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

suffix="_tagdust_arch.txt"
arch=$name$suffix
suffix="_tagdust";
outfile=$name$suffix


error=$( ../src/tagdust -t 80  -seed 42  $name -arch $arch -o $outfile )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  tagdust SUCCESS;
else
printf "%20s%10s\n"  tagdust FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi


suffix="_fastxbarcodefile.txt"
arch=$name$suffix
suffix="_fastx";
outfile=$name$suffix

error=$(  cat $name | fastx_barcode_splitter.pl --bcfile $arch --prefix $outfile -bol   2>&1)
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  fastx_split SUCCESS;
else
printf "%20s%10s\n"  fastx_split FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi







suffix="_tagdust_BC_";
outfile=$name$suffix

error=$( ../src/evalres -name tagdust_${array[$j]}_${numbar[$c]}_$barcodes  $outfile*.fq -o barread  -sim_random_frac 0.1 -sim_numseq 100000 2>&1)
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  eval_run SUCCESS;
else
printf "%20s%10s\n"  eval_run FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi


suffix="_fastxBC";
outfile=$name$suffix


error=$( ../src/evalres -name fastx_${array[$j]}_${numbar[$c]}_$barcodes $outfile* -o barread  -sim_random_frac 0.1 -sim_numseq 100000 2>&1)
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  eval_run SUCCESS;
else
printf "%20s%10s\n"  eval_run FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

let j++
let myrun++
done

done


#






