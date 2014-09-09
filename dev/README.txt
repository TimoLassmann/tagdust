Casava regression tests..

Data for this test was obtained from the validation dataset in bcl2fastq-1.8.4.tar:

With several modifications: 

1) install bcl2fastq

2) modify configuration makefile: 

vi /home/san/tlassmann/programs/bcl2fastqgo/share/bcl2fastq-1.8.4/examples/Validation/110120_P20_0993_A805CKABXX/Unaligned/ValidationConfig.mk

from: 

CONFIGURE_BCL_TO_FASTQ_PARAMS := --use-bases-mask 'y76n*,I6n,y76n*'

to: 

CONFIGURE_BCL_TO_FASTQ_PARAMS := --use-bases-mask 'y76n*,Y6n,y76n*'

(i.e. read 2 is the barcode and is kept in the output - written as compressed fastq file. 

3) modify sample sheet:

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


4) run the pipeline  & grab top 100k reads









BARCODE files in this directory were generared using the python script design_edit_metric_tags.py:

`

##############################################################################
#  design_edit_metric_tags.py  - generate sequence tags of arbitrary length  #
#                                                                            #
#  Copyright (c) 2009-2011 Brant C. Faircloth                                #
#  621 Charles E. Young Drive                                                #
#  University of California, Los Angeles, 90095, USA                         #
##############################################################################


Usage: design_edit_metric_tags.py [options]

Options:
  -h, --help            show this help message and exit
  --output=FILE         The path to the file where you want to store the
                        barcodes
  --tag-length=TL       The desired tag length
  --edit-distance=ED    The desired edit distance
  --multiprocessing     Use multiprocessing
  --processors=NPROCS   The number of processing cores to use when using
                        multiprocessing.  Default is # of cores - 2
  --no-polybase         Remove tags with > 2 identical nucleotides in a row
  --gc                  Remove tags with GC content (%) 40 > x > 60
  --comp                Remove tags that are perfect self-complements
  --hamming             Use Hamming distance in place of edit (Levenshtein)
                        distance.
  --min-and-greater     Show tags at all integer values of edit distance >
                        that specified
  --rescan=RESCAN       Rescan a file
  --rescan-length=RESCAN_LENGTH
                        Rescan length

