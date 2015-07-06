#
# Makefile for alignment library
#
# $Id: Makefile,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
#
#
# Copyright (c) 2003 Luke Sheneman.  All rights reserved.
#
# Luke Sheneman
# sheneman@cs.uidaho.edu
#
#

CC	= gcc
CFLAGS = -O3 -g3 -Wall -funroll-loops -march=pentium4 -mcpu=pentium4 -ffast-math -fprefetch-loop-arrays -finline-limit=5000 -fomit-frame-pointer -msse2 -mmmx -fno-strict-aliasing
#CFLAGS = -O3 -Wall -funroll-loops -march=pentium4 -mcpu=pentium4 -fno-unsafe-math-optimizations -fprefetch-loop-arrays -finline-limit=5000 -fomit-frame-pointer -fno-strict-aliasing
#CFLAGS = -O3 -g3 -Wall -funroll-loops -finline-limit=5000 -fomit-frame-pointer -fno-strict-aliasing -fno-unsafe-math-optimizations
#CFLAGS = -pg -Wall -fno-inline -fno-strict-aliasing
#CFLAGS = -g3 -Wall -fno-strict-aliasing
#CLINKER	= /usr/mpi-beowulf/bin/mpicc
INCDIR  = ../prng
LIBS	=	
DEFS 	=

SRCS    = align.c msf.c aln.c fasta.c submat.c syms.c symtab.c cache.c
OBJS	= align.o msf.o aln.o fasta.o submat.o syms.o symtab.o cache.o


all: align.a

clean:
	rm -f *.o a.out *.a libalign.a

.c.o:
	$(CC) -c $(CFLAGS) -I$(INCDIR) $(OPTFLAGS) $(DEFS) $<

align.a: $(OBJS)
	rm -f libalign.a
	ar -rs libalign.a $(OBJS)


