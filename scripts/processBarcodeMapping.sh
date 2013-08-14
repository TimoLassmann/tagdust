#!/bin/sh

#  processBarcodeMapping.sh
#  tagdust2
#
#  Created by lassmann on 8/12/13.
#  Copyright (c) 2013 lassmann. All rights reserved.


samtools view -b -q 20 all_mapping.bam | bamToBed  -i -   | intersectBed -a stdin -b ~/db/hg19_ann/refGene_merged.bed  -s -wao | awk 'BEGIN{line = 0}{\
line++;\
if($8 != -1){\
	x = split($4,a,":");\
	y = split($10,b,",");\
	table[a[x],b[1]]++;\
}\
}\
END{\
for (var in table){\
printf "%s\t%s\t%d\n",substr(var,1,6),substr(var,6) , table[var];\
}\
}'




#chr1    6337332 6337381 GA2X:1:58:18647:5607#0/1;BC:GAGATC      30      -       chr1    6324331 6453826 ACOT7,NM_181866;ACOT7,NM_181865;ACOT7,NM_181864;ACOT7,NM_007274      0       -       49

#chr1    4191075 4191122 GA2X:5:31:7689:1176#0/1;BC:CTGAAA       39      +       .       -1      -1      .       -1      .       0
