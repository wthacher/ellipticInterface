# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

makefiles += periodicTest


USE_PETSC=TRUE
MPI=TRUE
PETSC_DIR = ~/Desktop/LBL/petsc

# location of Chombo lib dir

ebase = periodicTest

CHOMBO_HOME = ../../../Chombo/lib

LibNames =  AMRElliptic AMRTools BoxTools BaseTools #EBAMRTools Workshop EBTools

base_dir = .
src_dirs = ../../src

include $(CHOMBO_HOME)/mk/Make.example
