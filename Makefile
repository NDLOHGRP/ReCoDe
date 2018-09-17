# ------------------------------------------------------------------
# This file is part of ReCoDe, a program and library for
# lossless, direct electron microscopy data compression.
#
# This program is released under the terms of the license contained
# in the file LICENSE.
# ------------------------------------------------------------------

SHELL="/bin/sh"

# To assist in cross-compiling
IDIR=include
ODIR=obj
LDIR=lib
SRC=src

CC=gcc

CFLAGS=-I$(IDIR) -O3
LIBS=-lz -lm -fopenmp
LIBS_M=-lz -lbz2 -llzma -lsnappy -llz4 -lm -fopenmp
OBJ=bin/recode

# Where you want it installed when you do 'make install'
PREFIX=/usr/local

recode: $(OBJ)
	$(CC) src/recode.c -o $(OBJ) $(CFLAGS) $(LIBS) -D ENABLE_MULTIPLE_COMPRESSIONS=0
	
recode_m: $(OBJ)
	$(CC) src/recode.c -o $(OBJ) $(CFLAGS) $(LIBS_M) -D ENABLE_MULTIPLE_COMPRESSIONS=1