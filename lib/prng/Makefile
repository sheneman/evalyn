#
# Makefile 

CC	= gcc
CFLAGS = -O3 -g3 -Wall -funroll-loops -march=pentium4 -mcpu=pentium4 -ffast-math -fprefetch-loop-arrays -finline-limit=5000 -fomit-frame-pointer -msse2 -mmmx -fno-strict-aliasing
#CFLAGS = -pg -g3 -Wall -fno-strict-aliasing
#CFLAGS	= -g3 -Wall -fno-strict-aliasing
#CFLAGS	= -O3 -Wall -funroll-loops -march=i686 -mcpu=i686 -fno-strict-aliasing
#CFLAGS = -O3 -Wall -funroll-loops -march=pentium4 -mcpu=pentium4 -fno-unsafe-math-optimizations -fprefetch-loop-arrays -finline-limit=5000 -fomit-frame-pointer -fno-strict-aliasing
#CFLAGS = -Wall -fno-strict-aliasing
LIBDIRS = -L.
LIBS	= -lprng -lm
INCDIRS =
DEFS	=


SRCS    = prng.c ptest.c
OBJS	= prng.o ptest.o
LIBOBJS = prng.a

all: prng.a ptest

clean:
	rm -f *.o a.out *.a ptest

.c.o:
	$(CC) -c $(INCDIRS) $(CFLAGS) $(DEFS) $<

ptest: $(OBJS) $(LIBOBJS)
	$(CC) $(LINKFLAGS) -o $@ ptest.o $(LIBDIRS) $(LIBS)

prng.a: prng.o
	rm -f libprng.a
	ar -rs libprng.a prng.o




