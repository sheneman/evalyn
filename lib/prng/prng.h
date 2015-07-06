/*
 * $Id: prng.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 * 
 * prng.h -- Mersenne Twisner Pseudo Random Number Generator
 *
 ********************************************************************************
 *
 * 
 * This implementation of the Mersenne Twister random number generator is based
 * on the reference implementation made available by the originators.  The
 * original copyright and license follow:
 *
 * ==============================================================================
 *
 * Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. The names of its contributors may not be used to endorse or promote
 *    products derived from this software without specific prior written
 *    permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <stdint.h>
#include <sys/types.h>


// Period parameters.
#define MT_N 624
#define MT_M 397


#define MT_BUFSIZE (MT_N * sizeof(uint32_t))

#define MT_EXP_CUTOFF 10.0



// Mersenne Twister random number generator state.
typedef struct _struct_MT_state {

  // The array for the state vector.
  uint32_t mt[MT_N];
  int32_t mti;

} MT_state;




// Seed the random number generator.
void
MT_init_genrand(uint32_t aSeed);

// copy state info
void
MT_copy_state(MT_state *dest_state,
	      MT_state *src_state);

void
MT_get_state(MT_state *ret_state);

void
MT_set_state(MT_state *new_state);

// Get a random integer in the range [0,2^63).
int64_t
MT_get_sint64(void);

// Get a random integer in the range [0,aRange).
int64_t
MT_get_range_sint64(int64_t aRange);


// Get a random integer in the range [0,0x7fffffff].
int32_t
MT_get_sint32(void);

// Get a random integer in the range [0,aRange).
int32_t
MT_get_range_sint32(int32_t aRange);

// Get a random integer in the range [0,0xffffffff].
uint32_t
MT_get_uint32(void);

// Get a random integer in the range [0,aRange).
// so if aRange = 5, it will chose values from 0 thru 4
uint32_t
MT_get_range_uint32(uint32_t aRange);


// Get a random real in the range [0,1].
double
MT_get_real1(void);


// Get a random real in the range [0,1).
double
MT_get_real2(void);


// Get a random real in the range (0,1).
double
MT_get_real3(void);

/* get an exponentially distributed number between 0 and 1 */
double
MT_get_bounded_exp(void);


void
MT_print_state(MT_state *state);

int
MT_state_cmp(MT_state *state1,
	     MT_state *state2);

long int
MT_get_sint32_top(long int top);

