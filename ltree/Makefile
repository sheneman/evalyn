#
# Makefile 
#
# $Id: Makefile,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
#
#
# Copyright (c) 2003 University of Idaho.  All rights reserved.
#
# Luke Sheneman
# 03/08/03
# sheneman@cs.uidaho.edu
#
# VALGRIND ARGS:
# valgrind -v --workaround-gcc296-bugs=yes --leak-check=yes \
#          --show-reachable=yes --logfile-fd=1 ./tree \
#          inputs/n100_s75_l75.seq > foo
#
# --skin=cachegrind
#

CC	= gcc
CFLAGS = -O3 -g3 -Wall -funroll-loops -march=pentium4 -mcpu=pentium4 -fprefetch-loop-arrays -finline-limit=5000 -fomit-frame-pointer -ffast-math -msse2 -mmmx -fno-strict-aliasing
#CFLAGS	= -g3 -Wall -fno-strict-aliasing
#CFLAGS = -pg -g3 -Wall -fno-strict-aliasing
#CFLAGS = -O3 -Wall -funroll-loops -march=i686 -mcpu=i686 -fno-strict-aliasing
#CFLAGS = -O3 -Wall -funroll-loops -march=pentium4 -mcpu=pentium4 -fno-unsafe-math-optimizations -fprefetch-loop-arrays -finline-limit=5000 -fomit-frame-pointer -fno-strict-aliasing


LIBDIRS = -L../lib/align -L../lib/prng -L../nj2
LIBS	= -lalign -lprng -lclearcut -lm
PROF_LIBC = /usr/lib/libc_p.a
INCDIRS = -I /usr/local/include -I ../lib/align -I ../lib/prng -I ../nj2
LINKFLAGS = -g3
#DEFS 	= -D __DEBUG__
DEFS	=


SRCS    = ltree.c dnd.c parsimony.c rnjhooks.c
OBJS	= ltree.o dnd.o parsimony.o rnjhooks.o
LIBOBJS = libalign.a libprng.a libclearcut.a

all: $(LIBOBJS) ltree

clean:
	(cd ../lib/align && make clean)
	(cd ../lib/prng && make clean)
	(cd ../nj2 && make clean)
	rm -f *.o a.out *.a ltree

libalign.a:
	(cd ../lib/align && make)

libprng.a:
	(cd ../lib/prng && make)

libclearcut.a:
	(cd ../nj2 && make)

.c.o:
	$(CC) -c $(INCDIRS) $(CFLAGS) $(DEFS) $<

ltree: $(OBJS) $(LIBOBJS)
	 $(CC) $(LINKFLAGS) -o $@ $(OBJS) $(LIBDIRS) $(LIBS)


