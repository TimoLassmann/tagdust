#!/bin/bash

#  make_benchmark_datasets.sh
#  
#
#  Created by timo on 16/06/2011.
#  Copyright 2011 RIKEN Yokohama Institute, Genome Exploration Research Group. All rights reserved.

DIRECTORY=
QTHRESHOLD=



function usage()
{
cat <<EOF
usage: $0 -d <directory> -q <threshold>
EOF
exit 1;
}

while getopts d:q: opt
do
case ${opt} in
d) DIRECTORY=${OPTARG};;
q) QTHRESHOLD=${OPTARG};;
*) usage;;
esac
done

if [ "${DIRECTORY}" = "" ]; then usage; fi
if [ "${QTHRESHOLD}" = "" ]; then usage; fi

#QSUB=`type qsub 2> /dev/null | sed -e "s/.* is //"`


array=(0.00 0.01 0.02 0.03 0.04 0.05)
len=${#array[*]}

numbar=(8 24 48 96)

myrun=0
j=0









for (( c=0; c < 4; c+=1 )); do
for (( i=4; i <= 8; i+=2 )); do
j=0

while [ $j -lt $len ]
do
if [ -n "${QSUB}" ]; then
echo "#!/bin/bash" > $myrun.tagdust.qsub;
echo "#\$ -cwd" >>  $myrun.tagdust.qsub;
echo "#\$ -N td$myrun" >>  $myrun.tagdust.qsub;
echo "#\$ -q osc-lm3.q@osc-lm3.gsc.riken.jp" >>   $myrun.tagdust.qsub;
echo "#\$ -e td$myrun.stderr" >>  $myrun.tagdust.qsub;
echo "#\$ -o td$myrun.stdout" >>  $myrun.tagdust.qsub;
echo "export PATH=/home/lassmann/bin:\$PATH" >>  $myrun.tagdust.qsub;
echo "/home/lassmann/bin/tagdust -sim_numseq 100000  -sim_5seq agggaggacgatgcgg   -sim_3seq gtgtcagtcacttccagcgg  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.1 -o $DIRECTORY  -sim_error_rate ${array[$j]}  -sim_InDel_frac 0.1 -sim_sequenced_len 1100  -sim_barlen $i  -sim_barnum ${numbar[$c]} " >>  $myrun.tagdust.qsub;
qsub $myrun.tagdust.qsub;
else
tagdust -sim_numseq 100000  -sim_5seq agggaggacgatgcgg   -sim_3seq gtgtcagtcacttccagcgg  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.1 -o $DIRECTORY  -sim_error_rate ${array[$j]}  -sim_InDel_frac 0.1 -sim_sequenced_len 1100  -sim_barlen $i  -sim_barnum ${numbar[$c]} -threshold $QTHRESHOLD
fi
let j++
let myrun++
done
done
done

for (( c=0; c < 4; c+=1 )); do
for (( i=4; i <= 8; i+=2 )); do
j=0

while [ $j -lt $len ]
do
if [ -n "${QSUB}" ]; then
echo "#!/bin/bash" > $myrun.tagdust.qsub;
echo "#\$ -cwd" >>  $myrun.tagdust.qsub;
echo "#\$ -N td$myrun" >>  $myrun.tagdust.qsub;
echo "#\$ -q osc-lm3.q@osc-lm3.gsc.riken.jp" >>   $myrun.tagdust.qsub;
echo "#\$ -e td$myrun.stderr" >>  $myrun.tagdust.qsub;
echo "#\$ -o td$myrun.stdout" >>  $myrun.tagdust.qsub;
echo "export PATH=/home/lassmann/bin:\$PATH" >>  $myrun.tagdust.qsub;
echo "tagdust -sim_numseq 100000  -sim_5seq agggaggacgatgcgg   -sim_3seq gtgtcagtcacttccagcgg  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.1 -o $DIRECTORY  -sim_error_rate ${array[$j]}  -sim_InDel_frac 0.5 -sim_sequenced_len 1100  -sim_barlen $i  -sim_barnum ${numbar[$c]} " >>  $myrun.tagdust.qsub;
qsub $myrun.tagdust.qsub;
else
tagdust -sim_numseq 100000  -sim_5seq agggaggacgatgcgg   -sim_3seq gtgtcagtcacttccagcgg  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.1 -o $DIRECTORY  -sim_error_rate ${array[$j]}  -sim_InDel_frac 0.5 -sim_sequenced_len 1100  -sim_barlen $i  -sim_barnum ${numbar[$c]} -threshold $QTHRESHOLD
fi
let j++
let myrun++
done
done
done






j=0

while [ $j -lt $len ]
do
if [ -n "${QSUB}" ]; then
echo "#!/bin/bash" > $myrun.tagdust.qsub;
echo "#\$ -cwd" >>  $myrun.tagdust.qsub;
echo "#\$ -N td$myrun" >>  $myrun.tagdust.qsub;
echo "#\$ -q osc-lm3.q@osc-lm3.gsc.riken.jp" >>   $myrun.tagdust.qsub;
echo "#\$ -e td$myrun.stderr" >>  $myrun.tagdust.qsub;
echo "#\$ -o td$myrun.stdout" >>  $myrun.tagdust.qsub;
echo "export PATH=/home/lassmann/bin:\$PATH" >>  $myrun.tagdust.qsub;
echo "tagdust -sim_numseq 100000  -sim_5seq agggaggacgatgcgg   -sim_3seq gtgtcagtcacttccagcgg  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.1 -o $DIRECTORY  -sim_error_rate ${array[$j]}  -sim_InDel_frac 0.1 -sim_sequenced_len 1100" >>  $myrun.tagdust.qsub;
qsub $myrun.tagdust.qsub;
else
tagdust -sim_numseq 100000  -sim_5seq agggaggacgatgcgg   -sim_3seq gtgtcagtcacttccagcgg  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.1 -o $DIRECTORY  -sim_error_rate ${array[$j]}  -sim_InDel_frac 0.1 -sim_sequenced_len 1100 -threshold $QTHRESHOLD
fi

let j++
let myrun++
done


j=0

while [ $j -lt $len ]
do
if [ -n "${QSUB}" ]; then

echo "#!/bin/bash" > $myrun.tagdust.qsub;
echo "#\$ -cwd" >>  $myrun.tagdust.qsub;
echo "#\$ -N td$myrun" >>  $myrun.tagdust.qsub;
echo "#\$ -q osc-lm3.q@osc-lm3.gsc.riken.jp" >>   $myrun.tagdust.qsub;
echo "#\$ -e td$myrun.stderr" >>  $myrun.tagdust.qsub;
echo "#\$ -o td$myrun.stdout" >>  $myrun.tagdust.qsub;
echo "export PATH=/home/lassmann/bin:\$PATH" >>  $myrun.tagdust.qsub;
echo "tagdust -sim_numseq 100000  -sim_5seq agggaggacgatgcgg   -sim_3seq gtgtcagtcacttccagcgg  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.1 -o $DIRECTORY  -sim_error_rate ${array[$j]}  -sim_InDel_frac 0.5 -sim_sequenced_len 1100" >>  $myrun.tagdust.qsub;
qsub $myrun.tagdust.qsub;
else
tagdust -sim_numseq 100000  -sim_5seq agggaggacgatgcgg   -sim_3seq gtgtcagtcacttccagcgg  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.1 -o $DIRECTORY  -sim_error_rate ${array[$j]}  -sim_InDel_frac 0.5 -sim_sequenced_len 1100 -threshold $QTHRESHOLD
fi
let j++
let myrun++
done


for (( c=0; c < 4; c+=1 )); do
for (( i=4; i <= 8; i+=2 )); do
j=0

while [ $j -lt $len ]
do
if [ -n "${QSUB}" ]; then
echo "#!/bin/bash" > $myrun.tagdust.qsub;
echo "#\$ -cwd" >>  $myrun.tagdust.qsub;
echo "#\$ -N td$myrun" >>  $myrun.tagdust.qsub;
echo "#\$ -q osc-lm3.q@osc-lm3.gsc.riken.jp" >>   $myrun.tagdust.qsub;
echo "#\$ -e td$myrun.stderr" >>  $myrun.tagdust.qsub;
echo "#\$ -o td$myrun.stdout" >>  $myrun.tagdust.qsub;
echo "export PATH=/home/lassmann/bin:\$PATH" >>  $myrun.tagdust.qsub;
echo "tagdust -sim_numseq 100000  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.1 -o $DIRECTORY  -sim_error_rate ${array[$j]}  -sim_InDel_frac 0.1 -sim_sequenced_len 1100  -sim_barlen $i  -sim_barnum ${numbar[$c]}" >>  $myrun.tagdust.qsub;
qsub $myrun.tagdust.qsub;
else
 tagdust -sim_numseq 100000  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.1 -o $DIRECTORY  -sim_error_rate ${array[$j]}  -sim_InDel_frac 0.1 -sim_sequenced_len 1100  -sim_barlen $i  -sim_barnum ${numbar[$c]} -threshold $QTHRESHOLD
fi
let j++
let myrun++
done
done
done


for (( c=0; c < 4; c+=1 )); do
for (( i=4; i <= 8; i+=2 )); do

j=0

while [ $j -lt $len ]
do
if [ -n "${QSUB}" ]; then
echo "#!/bin/bash" > $myrun.tagdust.qsub;
echo "#\$ -cwd" >>  $myrun.tagdust.qsub;
echo "#\$ -N td$myrun" >>  $myrun.tagdust.qsub;
echo "#\$ -q osc-lm3.q@osc-lm3.gsc.riken.jp" >>   $myrun.tagdust.qsub;
echo "#\$ -e td$myrun.stderr" >>  $myrun.tagdust.qsub;
echo "#\$ -o td$myrun.stdout" >>  $myrun.tagdust.qsub;
echo "export PATH=/home/lassmann/bin:\$PATH" >>  $myrun.tagdust.qsub;
echo "tagdust -sim_numseq 100000  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.1 -o $DIRECTORY  -sim_error_rate ${array[$j]}  -sim_InDel_frac 0.5 -sim_sequenced_len 1100  -sim_barlen $i  -sim_barnum ${numbar[$c]}" >>  $myrun.tagdust.qsub;
qsub $myrun.tagdust.qsub;
else

tagdust -sim_numseq 100000  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.1 -o $DIRECTORY  -sim_error_rate ${array[$j]}  -sim_InDel_frac 0.5 -sim_sequenced_len 1100  -sim_barlen $i  -sim_barnum ${numbar[$c]} -threshold $QTHRESHOLD

fi
let j++
let myrun++
done
done
done





