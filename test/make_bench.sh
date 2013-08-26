#!/bin/bash

#  make_benchmark_datasets.sh
#  
#
#  Created by timo on 16/06/2011.
#  Copyright 2011 RIKEN Yokohama Institute, Genome Exploration Research Group. All rights reserved.

array=(0.00 0.01 0.02 0.03 0.04 0.05)
len=${#array[*]}

for (( i=4; i <= 10; i+=2 )); do
j=0

while [ $j -lt $len ]
do
tagdust -simulation $i  -numbarcode $NUMBARCODE  -e ${array[$j]}  -i $INDEL -o  $base.$NUMBARCODE.$i.${array[$j]}.$INDEL;
let j++
done
done