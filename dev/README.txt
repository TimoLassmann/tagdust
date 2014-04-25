
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

