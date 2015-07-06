#
# Makefile 
#
# $Id: Makefile,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
#
#
#*****************************************************************************
#
# Copyright (c) 2004,  Luke Sheneman
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions 
# are met:
# 
#  + Redistributions of source code must retain the above copyright 
#    notice, this list of conditions and the following disclaimer. 
#  + Redistributions in binary form must reproduce the above copyright 
#    notice, this list of conditions and the following disclaimer in 
#    the documentation and/or other materials provided with the 
#    distribution. 
#  + The names of its contributors may not be used to endorse or promote 
#    products derived from this software without specific prior 
#    written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.  
#
#*****************************************************************************
#
# AUTHOR:
# 
#   Luke Sheneman
#   sheneman@cs.uidaho.edu
#
#*****************************************************************************


# Specify your compiler here
CC	= gcc
#CC	= cc


# Uncomment the following line to use the GNU getopt_long() instead
#  of the built-in getopt_long()

#OPTFLAG = -DUSE_GNU


############################# THE CFLAGS ####################################
#
# NOTE: The optimal selection of compiler flags is important.
#       Clearcut was written and developed on a Pentium 4 Linux
#       Computer.  As such, compiler optimization flags for 
#       recent versions of GCC on the Pentium 4 platform have been
#       well tested and do make a significant impact on overall 
#       performance.  
#  
#       Several alternative sets of compiler options are listed below.
#       Customize the compiler flags (CFLAGS) to match your target system.
#


######################## ALL ARCHITECTURES ##################################

# DEFAULT GCC OPTIMIZATION CONFIGURATION (ALL ARCHITECTURES)
#CFLAGS = -O3 -Wall -funroll-loops -fomit-frame-pointer -fno-strict-aliasing

# BASIC OPTIMIZATION
#CFLAGS = -O3

######################## INTEL PENTIUM OPTIMIZATIONS ########################

# PENTIUM 4 w/SSE GCC OPTIMIZATIONS.  REQUIRES GCC >= 3.2
CFLAGS = -O3 -Wall -funroll-loops -march=pentium4 -mcpu=pentium4 -mfpmath=sse -ffast-math -momit-leaf-frame-pointer -fno-strict-aliasing

# GENERIC i686 OPTIMIZATIONS GCC OPTIMIZATIONS
#CFLAGS = -O3 -Wall -funroll-loops -march=i686 -mcpu=i686 -fomit-frame-pointer -fno-strict-aliasing

# GENERIC Intel C Compiler (ICC) Optimizations
#CFLAGS = -O3 -axP -xP -unroll100 -IPF_fp_fast -IPF_fltacc+ -ip

#############################################################################


######################### POWER PC OPTIMIZATIONS ############################

# MAC OS X POWER PC OPTIMIZATIONS
#CFLAGS = -O3 -Wall -funroll-loops -fomit-frame-pointer -mcpu=powerpc -mabi=altivec -mtune=powerpc -maltivec

#############################################################################


# DEBUGGING COMPILER FLAGS
#CFLAGS	= -g3 -Wall


INCDIRS = -I /usr/local/include -I ../lib/prng -I ../lib/getopt_long
LIBDIRS = -L ../lib/prng -L ../lib/getopt_long

# *** THE LINK FLAGS ***
LINKFLAGS = -O3
LIBS = -lprng -lgetopt_long -lm

SRCS = nj_main.c
OBJS = nj_main.o
HDRS = nj_main.h

LIBSRCS = nj_clearcut.c nj_dist.c nj_dmat.c nj_fasta.c nj_cmdargs.c nj_bootstrap.c
LIBOBJS	= nj_clearcut.o nj_dist.o nj_dmat.o nj_fasta.o nj_cmdargs.o nj_bootstrap.o
LIBHDRS = nj_clearcut.h nj_dist.h nj_dmat.h nj_fasta.h nj_cmdargs.h nj_bootstrap.h nj_common.h nj_dayhoff.h

LIBCLEARCUT = libclearcut.a


all: libprng.a libgetopt_long.a $(LIBCLEARCUT) clearcut 

clean:
	(cd ../lib/prng && make clean)
	(cd ../lib/getopt_long && make clean)
	rm -f *.o a.out *.a clearcut

.c.o:   $(SRCS) $(LIBSRCS) $(HDRS) $(LIBHDERS)
	$(CC) -c $(INCDIRS) $(CFLAGS) $(OPTFLAG) $<

libprng.a:
	(cd ../lib/prng && make)

libgetopt_long.a:
	(cd ../lib/getopt_long && make)

libclearcut.a: $(LIBOBJS)
	rm -rf libclearcut.a
	ar -rs libclearcut.a $(LIBOBJS)

clearcut: $(OBJS) $(LIBCLEARCUT)
	 $(CC) $(LINKFLAGS) -o $@ $(OBJS) $(LIBCLEARCUT) $(LIBDIRS) $(LIBS)






 
