% TagDust2 benchmarking
% Timo Lassmann
% 11 Sepember 2014


# Introduction

Analysis of next generation sequencing datasets is complicated and often hard to reproduce. This package contains  a pipeline system to process Fucci single cell data in a _totally_ reproducible manner. Moreover, the package contains all the code, scripts and parameters to run the analysis. Typing __make__ in the root directory will actually start the entire analysis automatically.

For convenience during development all steps are executed within __makefiles__ ensuring that computationally expensive initial steps do not have to be repeated if downstream analysis steps are modified. The makefiles also keep the documentation up to date.

Finally I use the literal programming paradigm to describe and highlight the important steps. This file actually contains the code needed for the analysis. 

# Some details

The package is organized as follows:

1. the configure.ac file contains the instructions to check the versions of external programs used in the pipeline. As a default if the installed version is higher than the required version a warning will be produced. Otherwise the configuration will fail. 
2. once configured, make will compile the C code used to extract code from this document.
3. After all programs are extracted ro the __src__ directory, a makefile will call the __run.mk__ makefile and start the pipelines. 

# Installation and Usage

### Files
Here is the code of the  actual scripts used in the analysis.


#### File: run.mk

The makefile __run.mk__ executes the pipeline __process_run.mk__ on all target libraries. Afterwards additional makefiles collect output files and summarize the results.

~~~~{.Makefile}



.PHONY: mapping 

all: mapping

mapping:
	@echo GAGA;

~~~~



