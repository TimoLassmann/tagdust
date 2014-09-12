% TagDust2 Active Supplement
% Timo Lassmann
% 11 Sepember 2014


# Introduction

The purpose of this document is to provide a complete description and the supporting code to reproduce all figures from the TagDust2 paper.

For convenience during development all steps are executed within __makefiles__ ensuring that computationally expensive initial steps do not have to be repeated if downstream analysis steps are modified. The makefiles also keep the documentation up to date.

Finally I use the literal programming paradigm to describe and highlight the important steps. This file actually contains the code needed for the analysis. 


The package is organized as follows:

1. the configure.ac file contains the instructions to check the versions of external programs used in the pipeline. As a default if the installed version is higher than the required version a warning will be produced. Otherwise the configuration will fail. 
2. once configured, make will compile the C code used to extract code from this document.
3. After all programs are extracted ro the __src__ directory, a makefile will call the __run.mk__ makefile and start the pipelines. 

# Installation and Usage

## Files
Here is the code of the  actual scripts used in the analysis.


### File: run.mk

The makefile __run.mk__ executes the pipeline __process_run.mk__ on all target libraries. Afterwards additional makefiles collect output files and summarize the results.

~~~~{.Makefile}

.PHONY:  message benchmark

benchmark: rplot
	R --slave --vanilla < plotting.R
	rm Rplots.pdf
	@echo HelloWorld

rplot: barread_4nt_4r.tsv 5barread3_4nt_4r.tsv barread_6nt_4r.tsv 5barread3_6nt_4r.tsv
	@echo HelloWorld

barread_4nt_4r.tsv: barread_results.txt
	cat barread_results.txt | grep _4nt_ | awk 'BEGIN{print "Program\tsimerror\tbarcodes\tRecall\tSpecificity\tPrecision\tKappa"} {if(NR>0){x = split($$1,a,"_");printf "%s\t%f\t%f\t%f\t%f\t%f\t%f\n" ,  a[1],a[2],a[3],$$2,$$3,$$4,$$5 }}' > barread_4nt_4r.tsv
	
5barread3_4nt_4r.tsv: 5barread3_results.txt
	cat 5barread3_results.txt |  grep _4nt_ | awk 'BEGIN{print "Program\tsimerror\tbarcodes\tRecall\tSpecificity\tPrecision\tKappa"} {if(NR>0){x = split($$1,a,"_");printf "%s\t%f\t%f\t%f\t%f\t%f\t%f\n" ,  a[1],a[2],a[3],$$2,$$3,$$4,$$5 }}' > 5barread3_4nt_4r.tsv

barread_6nt_4r.tsv: barread_results.txt
	cat barread_results.txt | grep _6nt_ | awk 'BEGIN{print "Program\tsimerror\tbarcodes\tRecall\tSpecificity\tPrecision\tKappa"} {if(NR>0){x = split($$1,a,"_");printf "%s\t%f\t%f\t%f\t%f\t%f\t%f\n" ,  a[1],a[2],a[3],$$2,$$3,$$4,$$5 }}' > barread_6nt_4r.tsv
	
5barread3_6nt_4r.tsv: 5barread3_results.txt
	cat 5barread3_results.txt |  grep _6nt_ | awk 'BEGIN{print "Program\tsimerror\tbarcodes\tRecall\tSpecificity\tPrecision\tKappa"} {if(NR>0){x = split($$1,a,"_");printf "%s\t%f\t%f\t%f\t%f\t%f\t%f\n" ,  a[1],a[2],a[3],$$2,$$3,$$4,$$5 }}' > 5barread3_6nt_4r.tsv



barread_results.txt: all_arch.txt
	./barread.sh -b ../../dev/EDITTAG_4nt_ed_2.txt -r
	./barread.sh -b ../../dev/EDITTAG_6nt_ed_3.txt -r
	
5barread3_results.txt: all_arch.txt
	./5barread3.sh -b ../../dev/EDITTAG_4nt_ed_2.txt -r
	./5barread3.sh -b ../../dev/EDITTAG_6nt_ed_3.txt -r

all_arch.txt: 
	./barread.sh -b ../../dev/EDITTAG_4nt_ed_2.txt
	./barread.sh -b ../../dev/EDITTAG_6nt_ed_3.txt
	./5barread3.sh -b ../../dev/EDITTAG_4nt_ed_2.txt
	./5barread3.sh -b ../../dev/EDITTAG_6nt_ed_3.txt
	cat  *tagdust_arch.txt  | sort | uniq  > all_arch.txt
	



all: message

message:
	@echo run make benchmark..;

~~~~

### File: barread.sh

~~~~{.bash}
#!/bin/bash

echo "Running tagdust benchmarks:";

barcodes=
run=0

function usage()
{
cat <<EOF
usage: $0 -b  <barcode file> -r
EOF
exit 1;
}

while getopts b:r opt
do
case ${opt} in
b) barcodes=${OPTARG};;
r) run=1;;
*) usage;;
esac
done

if [ "${barcodes}" = "" ]; then usage; fi

array=(0.01 0.015 0.02 0.025 0.03)
len=${#array[*]}

numbar=(8 24 48)
lenbar=${#numbar[*]}


if [[ $run -eq 1 ]]; then
make --silent clean
fi



for (( c=0; c < $lenbar; c+=1 )); do

j=0

while [ $j -lt $len ]
do

name="barread${array[$j]}_${numbar[$c]}.fq"

error=$(../../src/simreads $barcodes  -seed 42 -sim_barnum ${numbar[$c]}  -sim_readlen 20 -sim_readlen_mod 0 -sim_numseq 100000 -sim_endloss 0 -sim_random_frac 0.1 -o $name   -sim_error_rate ${array[$j]}  2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  sim SUCCESS;
else
printf "%20s%10s\n"  sim FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

if [[ $run -eq 1 ]]; then

suffix="_tagdust_arch.txt"
arch=$name$suffix
suffix="_tagdust";
outfile=$name$suffix


error=$( ../../src/tagdust -t 80  -seed 42  $name -arch $arch -o $outfile 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  tagdust SUCCESS;
else
printf "%20s%10s\n"  tagdust FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi


suffix="_tagdustarchsel";
outfile=$name$suffix


error=$( ../../src/tagdust -t 80  -seed 42  $name -arch all_arch.txt -o $outfile 2>&1 )
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

error=$(  cat $name | ../bin/fastx_barcode_splitter.pl --bcfile $arch --prefix $outfile -bol   2>&1)
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

error=$( ../../src/evalres -name tagdust_${array[$j]}_${numbar[$c]}_$barcodes  $outfile*.fq -o barread  -sim_random_frac 0.1 -sim_numseq 100000 2>&1)
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  eval_run SUCCESS;
else
printf "%20s%10s\n"  eval_run FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi



suffix="_tagdustarchsel_BC_";
outfile=$name$suffix

error=$( ../../src/evalres -name tagdustallarch_${array[$j]}_${numbar[$c]}_$barcodes  $outfile*.fq -o barread  -sim_random_frac 0.1 -sim_numseq 100000 2>&1)
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


error=$( ../../src/evalres -name fastx_${array[$j]}_${numbar[$c]}_$barcodes $outfile* -o barread  -sim_random_frac 0.1 -sim_numseq 100000 2>&1)
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  eval_run SUCCESS;
else
printf "%20s%10s\n"  eval_run FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi
fi

let j++
let myrun++
done

done

if [[ $run -eq 1 ]]; then
make --silent clean
fi

~~~~


### File: 5barread3.sh

~~~~{.bash}
#!/bin/bash

echo "Running tagdust benchmarks:";

barcodes=
run=0

function usage()
{
cat <<EOF
usage: $0 -b  <barcode file> -r
EOF
exit 1;
}

while getopts b:r opt
do
case ${opt} in
b) barcodes=${OPTARG};;
r) run=1;;
*) usage;;
esac
done

if [ "${barcodes}" = "" ]; then usage; fi

array=(0.01 0.015 0.02 0.025 0.03)
len=${#array[*]}

numbar=(8 24 48)
lenbar=${#numbar[*]}



if [[ $run -eq 1 ]]; then
make --silent clean
fi

#barcodes="EDITTAG_4nt_ed_2.txt"

for (( c=0; c < $lenbar; c+=1 )); do

j=0

while [ $j -lt $len ]
do

name="5barread3${array[$j]}_${numbar[$c]}.fq"

error=$(../../src/simreads $barcodes  -sim_5seq agggaggacgatgcgg   -sim_3seq gtgtcagtcacttccagcgg  -seed 42 -sim_barnum ${numbar[$c]}  -sim_readlen 20 -sim_readlen_mod 0 -sim_numseq 100000 -sim_endloss 0 -sim_random_frac 0.1 -o $name   -sim_error_rate ${array[$j]}  2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  sim SUCCESS;
else
printf "%20s%10s\n"  sim FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

if [[ $run -eq 1 ]]; then
suffix="_tagdust_arch.txt"
arch=$name$suffix
suffix="_tagdust";
outfile=$name$suffix


error=$( ../../src/tagdust -t 80  -seed 42  $name -arch $arch -o $outfile 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  tagdust SUCCESS;
else
printf "%20s%10s\n"  tagdust FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi


suffix="_tagdustarchsel";
outfile=$name$suffix


error=$( ../../src/tagdust -t 80  -seed 42  $name -arch all_arch.txt -o $outfile 2>&1 )
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  tagdust SUCCESS;
else
printf "%20s%10s\n"  tagdust FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi






suffix="_btrim_pattern.txt"
arch=$name$suffix
suffix="_btrim";
outfile=$name$suffix

error=$( ../bin/btrim64 -p $arch  -t $name -o $outfile -l 18 -B  2>&1)
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  btrim SUCCESS;
else
printf "%20s%10s\n"  btrim FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi

suffix="_cutadadapt.fq";
outfile=$name$suffix
error=$( ../bin/cutadapt --discard-untrimmed -n 2 --front agggaggacgatgcgg    --adapter=gtgtcagtcacttccagcgg  $name   -o  $outfile  2>&1)
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  cutadapt SUCCESS;
else
printf "%20s%10s\n"  cutadapt FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi


suffix="_cutadadapt.fq";
cutadaptout=$name$suffix
suffix="_fastxbarcodefile.txt"
arch=$name$suffix
suffix="_fastx";
outfile=$name$suffix

error=$(  cat $cutadaptout | ../bin/fastx_barcode_splitter.pl --bcfile $arch --prefix $outfile -bol   2>&1)
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

error=$( ../bin/../src/evalres -name tagdust_${array[$j]}_${numbar[$c]}_$barcodes  $outfile*.fq -o 5barread3  -sim_random_frac 0.1 -sim_numseq 100000 2>&1)
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  eval_run SUCCESS;
else
printf "%20s%10s\n"  eval_run FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi


suffix="_tagdustarchsel_BC_";
outfile=$name$suffix

error=$( ../../src/evalres -name tagdustallarch_${array[$j]}_${numbar[$c]}_$barcodes  $outfile*.fq -o 5barread3  -sim_random_frac 0.1 -sim_numseq 100000 2>&1)
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  eval_run SUCCESS;
else
printf "%20s%10s\n"  eval_run FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi





suffix="_btrim";
outfile=$name$suffix

error=$( ../../src/evalres -name btrim_${array[$j]}_${numbar[$c]}_$barcodes $outfile* -o 5barread3  -sim_random_frac 0.1 -sim_numseq 100000 2>&1)
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


error=$( ../../src/evalres -name fastx_${array[$j]}_${numbar[$c]}_$barcodes $outfile* -o 5barread3  -sim_random_frac 0.1 -sim_numseq 100000 2>&1)
status=$?
if [[ $status -eq 0 ]]; then
printf "%20s%10s\n"  eval_run SUCCESS;
else
printf "%20s%10s\n"  eval_run FAILED;
printf "with ERROR $status and Message:\n\n$error\n\n";
exit 1;
fi


fi;

let j++
let myrun++
done

done

if [[ $run -eq 1 ]]; then

make --silent clean
fi


~~~~ 


### File: plotting.R

Here is the R script to make the plots.

~~~~{.R}

library(ggplot2)
library(gplots)
library(scales)
require("reshape")

yaxis_label_format <- function(x) {
	lab <- sprintf('%0.2f', x) # Format the strings as HH:MM:SS
}


mat = read.table("barread_4nt_4r.tsv",header =T,sep="\t")



mat$barcodes = paste(mat$barcodes, "Barcodes")
mat$barcodes  = factor(mat$barcodes, levels=c('8 Barcodes','24 Barcodes','48 Barcodes'))

s = subset(mat,mat$Program != "tagdustallarch")
s$Program  = factor(s$Program, levels=c('tagdust','fastx','btrim'))

m = melt(s, id=c("Program","simerror","barcodes"),measure.vars = c("Recall", "Precision"))
ggplot(m, aes(x = simerror,y = value ,color=Program, fill=factor(Program))) + geom_line()  +facet_grid(scales="free", variable ~ barcodes  ) +   scale_x_continuous(breaks=c(0.01,0.015,0.02,0.025,0.03), labels = expression("1","1.5","2","2.5","3")) +  xlab("Simulated Error Rate (%)") + scale_y_continuous(labels=yaxis_label_format)+ ylab("")

ggsave("barread_4nt.pdf",width=8,height=4,dpi = 300)




mat = read.table("barread_6nt_4r.tsv",header =T,sep="\t")

mat$barcodes = paste(mat$barcodes, "Barcodes")
mat$barcodes  = factor(mat$barcodes, levels=c('8 Barcodes','24 Barcodes','48 Barcodes'))

s = subset(mat,mat$Program != "tagdustallarch")
s$Program  = factor(s$Program, levels=c('tagdust','fastx','btrim'))
m = melt(s, id=c("Program","simerror","barcodes"),measure.vars = c("Recall", "Precision"))
ggplot(m, aes(x = simerror,y = value ,color=Program, fill=factor(Program))) + geom_line()  +facet_grid(scales="free", variable ~ barcodes  ) +   scale_x_continuous(breaks=c(0.01,0.015,0.02,0.025,0.03), labels = expression("1","1.5","2","2.5","3")) +  xlab("Simulated Error Rate (%)") + scale_y_continuous(labels=yaxis_label_format)+ ylab("")

ggsave("barread_6nt.pdf",width=8,height=4,dpi = 300)





mat = read.table("5barread3_4nt_4r.tsv",header =T,sep="\t")
mat$barcodes = paste(mat$barcodes, "Barcodes")
mat$barcodes  = factor(mat$barcodes, levels=c('8 Barcodes','24 Barcodes','48 Barcodes'))
s = subset(mat,mat$Program != "tagdustallarch")
s$Program  = factor(s$Program, levels=c('tagdust','fastx','btrim'))
m = melt(s, id=c("Program","simerror","barcodes"),measure.vars = c("Recall", "Precision"))
ggplot(m, aes(x = simerror,y = value ,color=Program, fill=factor(Program))) + geom_line()  +facet_grid(scales="free", variable ~ barcodes  ) +   scale_x_continuous(breaks=c(0.01,0.015,0.02,0.025,0.03), labels = expression("1","1.5","2","2.5","3")) +  xlab("Simulated Error Rate (%)") + scale_y_continuous(labels=yaxis_label_format) + ylab("")

ggsave("5barread3_4nt.pdf",width=8,height=4,dpi = 300)


mat = read.table("5barread3_6nt_4r.tsv",header =T,sep="\t")

mat$barcodes = paste(mat$barcodes, "Barcodes")
mat$barcodes  = factor(mat$barcodes, levels=c('8 Barcodes','24 Barcodes','48 Barcodes'))
s = subset(mat,mat$Program != "tagdustallarch")
s$Program  = factor(s$Program, levels=c('tagdust','fastx','btrim'))
m = melt(s, id=c("Program","simerror","barcodes"),measure.vars = c("Recall", "Precision"))
ggplot(m, aes(x = simerror,y = value ,color=Program, fill=factor(Program))) + geom_line()  +facet_grid(scales="free", variable ~ barcodes  ) +   scale_x_continuous(breaks=c(0.01,0.015,0.02,0.025,0.03), labels = expression("1","1.5","2","2.5","3")) +  xlab("Simulated Error Rate (%)") + scale_y_continuous(labels=yaxis_label_format)+ ylab("")

ggsave("5barread3_6nt.pdf",width=8,height=4,dpi = 300)



~~~~




