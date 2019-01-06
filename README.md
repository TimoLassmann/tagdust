# Introduction

Raw sequences produced by next generation sequencing (NGS) machines contain adapter, linker,
barcode and fingerprint (UMI) sequences. TagDust2 is designed to make is as easy as possible to de-multiplex reads from any library preparation method. In addition, TagDust2 can detect which library prep. was used from the raw reads themselves making it possible to automate processing pipelines.


https://user-images.githubusercontent.com/8110320/50732573-c95c4e80-114b-11e9-852f-8ee817af20ba.jpg
![Image of example output](https://user-images.githubusercontent.com/8110320/50732573-c95c4e80-114b-11e9-852f-8ee817af20ba.jpg)

TagDust allows users to specify the expected architecture of a read and converts it into a hidden
Markov model. The latter can assign sequences to a particular barcode (or index) even in the presence
of sequencing errors. Sequences not matching the architecture (primer dimers, contaminants etc.) are
automatically discarded (see Figure 1.1).


# Installation: 

Unpack the tarball:
bash-3.1$ tar -zxvf tagdust-XXX.tar.gz
bash-3.1$ cd tagdust
bash-3.1$ ./autogen.sh
bash-3.1$ ./configure
bash-3.1$ make
bash-3.1$ make check

At this point the TagDust executable appears in the src directory. You can copy it to any directory in your path. To install it system wide type:

bash-3.1$ make install

# Manual

Have a look at the user manual in the doc directory! 

# Please cite:

Lassmann, Timo. "TagDust2: a generic method to extract reads from sequencing data." BMC bioinformatics 16.1 (2015): 24. 
[doi.org/10.1186/s12859-015-0454-y](https://doi.org/10.1186/s12859-015-0454-y)





