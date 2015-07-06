#
# Makefile for BSD-licensed getopt_long() implementation
#
#

CC	= gcc
CFLAGS = -O3 -Wall -fno-strict-aliasing
LIBDIRS = -L.
LIBS	= -lgetopt_long -lm
INCDIRS =
DEFS	=


SRCS    = getopt_long.c 
OBJS	= getopt_long.o
LIBOBJS = libgetopt_long.a

all: libgetopt_long.a

clean:
	rm -f *.o a.out *.a 

.c.o:
	$(CC) -c $(INCDIRS) $(CFLAGS) $(DEFS) $<

libgetopt_long.a: getopt_long.o
	rm -f libgetopt_long.a
	ar -rs libgetopt_long.a getopt_long.o




