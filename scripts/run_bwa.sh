#!/bin/bash
#
# Copyright (C) 2010,2011 Timo Lassmann <timolassmann@gmail.com>
# 
# This file is part of delve.
# 
# Delve is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Delve is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Delve.  If not, see <http://www.gnu.org/licenses/>.
#



GENOME_DIR= 

ALNDIR= 

function usage()
{
cat <<EOF
usage: $0  -g <genome file> -t <fasta directory>
EOF
exit 1;
}

while getopts g:t: opt
do
case ${opt} in
g) GENOME_DIR=${OPTARG};;
t) ALNDIR=${OPTARG};;
*) usage;;
esac
done




if [ "${GENOME_DIR}" = "" ]; then usage; fi
if [ "${ALNDIR}" = "" ]; then usage; fi



QSUB=`type qsub 2> /dev/null | sed -e "s/.* is //"`

for file in $ALNDIR/*                       
	do
		if [ -f $file ]; then 
			if [[ $file =~ _extracted.fq$ ]]; then
				if [ -f $$file.sorted.bam ]; then
					echo "$file already processed"
				else
						echo "working on $file"
#bwa   aln $GENOME_DIR $file -t 80  > tmp.sai
#bwa samse $GENOME_DIR tmp.sai $file > $file.sam

						 bwa aln -t 80  $GENOME_DIR $file  | bwa samse $GENOME_DIR - $file | samtools view -Su - | samtools sort - $file.sorted


#					fi
				fi
			fi
		fi
	done



cat *sam | awk -v thres=20  'BEGIN{for (i = 0; i <= 80; i++) {storage[$i] = 1;}}  {  if($2 != 4){ x = split($1,a,"[:,;]");  print a[9],$1 } }'


cat *sam | awk ' {  if($2 != 4){ x = split($1,a,"[:,;]");  printf "%f\t%f\n" ,$5,a[9] } }'








