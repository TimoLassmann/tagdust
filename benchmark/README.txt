
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



Casava comparison:

1) install bcl2fastq  (bcl2fastq-1.8.4.tar.gz)



2) modify:

/share/bcl2fastq-1.8.4/examples/Validation/110120_P20_0993_A805CKABXX/Data/Intensities/BaseCalls/SampleSheet.csv

so that all samples of lane one are de-multiplexed. (default is only PhiX and one of the human samples..) 

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





