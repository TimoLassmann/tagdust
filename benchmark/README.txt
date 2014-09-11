Type make to reproduce the paper figures....

Required programs: 

BTRIM: downloaded from: 
http://graphics.med.yale.edu/trim/
 btrim64                 28-Jul-2014 15:56   35K 

FASTX:
version 0.0.13

cutadapt:
./cutadapt --version
1.3

BWA:
Version: 0.7.8-r455

Casava comparison:

1) install bcl2fastq  (bcl2fastq-1.8.4.tar.gz)

2) modify:

/share/bcl2fastq-1.8.4/examples/Validation/110120_P20_0993_A805CKABXX/Data/Intensities/BaseCalls/SampleSheet.csv

so that all samples of lane 1 are de-multiplexed. (default is only PhiX and one of the human samples..)

3) run

configureValidation.pl \
         --source-dir ../share/bcl2fastq-1.8.4/examples/Validation/110120_P20_0993_A805CKABXX \
         --output-dir ./ValidationMultiplexed_casava

4) run basecaller without de-multiplexing.

5) modify configuration makefile:

vi share/bcl2fastq-1.8.4/examples/Validation/110120_P20_0993_A805CKABXX/Unaligned/ValidationConfig.mk

from: 

CONFIGURE_BCL_TO_FASTQ_PARAMS := --use-bases-mask 'y76n*,I6n,y76n*'

to: 

CONFIGURE_BCL_TO_FASTQ_PARAMS := --use-bases-mask 'y76n*,Y6n,y76n*'

(i.e. read 2 is the barcode and is kept in the output - written as compressed fastq file. 

6) modify sample sheet:

/share/bcl2fastq-1.8.4/examples/Validation/110120_P20_0993_A805CKABXX/Data/Intensities/BaseCalls/SampleSheet.csv  

from: 

whatever...

to:

FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator
A805CKABXX,1,not_demultiplexed,,,Test bcl conversion,N,D109LACXX,BCF,testbclconv
#A805CKABXX,1,AR005,human,ACAGTG,Cypress,Y,101+7,CB,Demo
#A805CKABXX,1,AR008,human,ACTTGA,Cypress,Y,101+7,CB,Demo
#A805CKABXX,1,PhiX,phix,TTAGGC,Cypress,Y,101+7,CB,Demo
#A805CKABXX,1,,unknown,Undetermined,Ignored clusters with unmatched barcodes for lane 1,N,101+7,CB,

Note that only lane 1 is included in the test. 


6) run the pipeline  


Now we have demulitplexed data from casava and be-multiplexed data for tagdust. 

##### 

The actual comparison: 

Number of extracted reads cassava: 

zcat  Sample_AR005/AR005_ACAGTG_L001_R1_001.fastq.gz | wc
5890936 /4 =
1472734

zcat  Sample_AR008/AR008_ACTTGA_L001_R1_001.fastq.gz  | wc
6072832 / 4 = 
1518208

zcat Sample_PhiX/PhiX_TTAGGC_L001_R1_001.fastq.gz  | wc
 196416 / 4 = 
49104


Number of extracted reads TagDust:

Run TagDust:

tagdust -arch casava_arch.txt not_demultiplexed_NoIndex_L001_R1_001.fastq.gz not_demultiplexed_NoIndex_L001_R2_001.fastq.gz  not_demultiplexed_NoIndex_L001_R3_001.fastq.gz  -t 80 -o demulti_large

Here is what's in casava_arch.txt:

--------------------
tagdust -1 B:ACAGTG,ACTTGA,TTAGGC
tagdust -1 R:N
--------------------

wc  -l demulti_large_*READ1.fq
6035392 demulti_large_BC_ACAGTG_READ1.fq
6175672 demulti_large_BC_ACTTGA_READ1.fq
205952 demulti_large_BC_TTAGGC_READ1.fq
/4 =

1508848
1543918
51488



Align:

TAgDust:

nice bwa  mem -t 16  /GROUP/REFERENCE/GRCh38/GRCh38.genome.fa demulti_large_BC_ACAGTG_READ1.fq demulti_large_BC_ACAGTG_READ2.fq   | samtools view -Sb - > tagdust_ACAGTG.bam

nice bwa  mem -t 16  /GROUP/REFERENCE/GRCh38/GRCh38.genome.fa demulti_large_BC_ACTTGA_READ1.fq demulti_large_BC_ACTTGA_READ2.fq   | samtools view -Sb - > tagdust_ACTTGA.bam

nice bwa  mem -t 16 phix.fa demulti_large_BC_TTAGGC_READ1.fq demulti_large_BC_TTAGGC_READ2.fq   | samtools view -Sb - > tagdust_TTAGGC.bam

nice bwa  mem -t 16  /GROUP/REFERENCE/GRCh38/GRCh38.genome.fa Sample_AR005/AR005_ACAGTG_L001_R1_001.fastq.gz Sample_AR005/AR005_ACAGTG_L001_R2_001.fastq.gz | samtools view -Sb - > casava_ACAGTG.bam

nice bwa  mem -t 16  /GROUP/REFERENCE/GRCh38/GRCh38.genome.fa Sample_AR008/AR008_ACTTGA_L001_R1_001.fastq.gz Sample_AR008/AR008_ACTTGA_L001_R2_001.fastq.gz | samtools view -Sb - > casava_ACTTGA.bam

nice bwa  mem -t 16 phix.fa Sample_PhiX/PhiX_TTAGGC_L001_R1_001.fastq.gz Sample_PhiX/PhiX_TTAGGC_L001_R2_001.fastq.gz  | samtools view -Sb - > casava_TTAGGC.bam

#count correctly paired reads (excluding secondary and supplementary alignments:)

samtools view -c -f 1 -F 2316 casava_ACAGTG.bam
samtools view -c -f 1 -F 2316  casava_ACTTGA.bam
samtools view -c -f 1 -F 2316 casava_TTAGGC.bam
samtools view -c -f 1 -F 2316 tagdust_ACAGTG.bam
samtools view -c -f 1 -F 2316 tagdust_ACTTGA.bam
samtools view -c -f 1 -F 2316 tagdust_TTAGGC.bam


