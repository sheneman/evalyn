/*
 * $Id: prng.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 * 
 * prng.c -- Mersenne Twister Pseudo Random Number Generator
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

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>

#include "prng.h"



// Constant vector a.
#define MT_matrix_a 0x9908b0dfUL

// Most significant w-r bits.
#define MT_upper_mask 0x80000000UL

// Least significant r bits.
#define MT_lower_mask 0x7fffffffUL


/* the global PRNG state */
MT_state MT_gbl_state;



void
MT_get_state(MT_state *ret_state) {

  MT_copy_state(ret_state, &MT_gbl_state);

  return;

}


void
MT_set_state(MT_state *new_state) {

  MT_copy_state(&MT_gbl_state, new_state);
  
  return;
}


// Initializes mt[MT_N] with a seed.
void
MT_init_genrand(uint32_t aS) {

  MT_state *aMt;
  uint32_t *mt;
  uint32_t mti;
  
  aMt = &MT_gbl_state;
  mt = aMt->mt;
  
  mt[0]= aS & 0xffffffffUL;
  for (mti=1; mti<MT_N; mti++) {
    mt[mti] =
      (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);

    // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
    // In the previous versions, MSBs of the seed affect
    // only MSBs of the array mt[].
    // 2002/01/09 modified by Makoto Matsumoto
    mt[mti] &= 0xffffffffUL;
    // for >32 bit machines
  }

  aMt->mti = mti;
}




// generates a random number on [0,0xffffffff]-interval
uint32_t
MT_get_uint32(void)  {

  MT_state *aMt;
  uint32_t *mt;
  uint32_t mti;
  uint32_t y;
  static uint32_t mag01[2]={0x0UL, MT_matrix_a};

  aMt = &MT_gbl_state;
  mt = aMt->mt;
  mti = aMt->mti;

  if (mti >= MT_N) { // generate MT_N words at one time
    int32_t kk;

    if (mti == MT_N+1) {  // if MT_init_genrand() has not been called,
      MT_init_genrand(5489UL); // a default initial seed is used
    }

    for (kk=0;kk<MT_N-MT_M;kk++) {
      y = (mt[kk]&MT_upper_mask)|(mt[kk+1]&MT_lower_mask);
      mt[kk] = mt[kk+MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (;kk<MT_N-1;kk++) {
      y = (mt[kk]&MT_upper_mask)|(mt[kk+1]&MT_lower_mask);
      mt[kk] = mt[kk+(MT_M-MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mt[MT_N-1]&MT_upper_mask)|(mt[0]&MT_lower_mask);
    mt[MT_N-1] = mt[MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    mti = 0;
  }

  y = mt[mti++];

  // Tempering
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  aMt->mti = mti;
  return y;
}






int64_t
MT_get_sint64(void) {

  int64_t rVal;

  rVal = (int64_t) MT_get_uint32();
  rVal &= 0xfffffffeU; // Clear least significant bit.
  rVal <<= 31;
  rVal |= ((int64_t) MT_get_uint32());

  return rVal;
}




void
MT_copy_state(MT_state *dest_state,
	      MT_state *src_state) {

  memcpy(dest_state, src_state, sizeof(MT_state));

  return;
}





int64_t
MT_get_range_sint64(int64_t aRange) {
		    

  int64_t rVal, above;
  
  above = 0x7fffffffffffffffLL - (0x7fffffffffffffffLL % aRange);
  while (1)
    {
      rVal = MT_get_sint64();
      if (rVal < above)
	{
	  rVal %= aRange;
	  break;
	}
    }

  return rVal;
}


int32_t
MT_get_sint32(void) {
  return (int32_t) (MT_get_uint32() >> 1);
}


int32_t
MT_get_range_sint32(int32_t aRange) {
		    

  int32_t rVal, above;
  
  above = 0x7fffffffL - (0x7fffffffL % aRange);
  while (1)
    {
      rVal = (int32_t) (MT_get_uint32()>>1);
      if (rVal < above)
	{
	  rVal %= aRange;
	  break;
	}
    }

  return rVal;
}



// generates a random number on [0,0xffffffff]-interval
uint32_t
MT_get_range_uint32(uint32_t aRange) {

  uint32_t rVal, above;

  above = 0xffffffffLU - (0xffffffffLU % aRange);
  while (1)
    {
      rVal = MT_get_uint32();
      if (rVal < above)
	{
	  rVal %= aRange;
	  break;
	}
    }

  return rVal;
}



// generates a random number on [0,1]-real-interval
double
MT_get_real1(void) {

  return MT_get_uint32()*(1.0/4294967295.0);
  // divided by 2^32-1
}

// generates a random number on [0,1)-real-interval
double
MT_get_real2(void) {

  return MT_get_uint32()*(1.0/4294967296.0);
  // divided by 2^32
}

// generates a random number on (0,1)-real-interval
double
MT_get_real3(void) {

  return (((double)MT_get_uint32()) + 0.5)*(1.0/4294967296.0);
  // divided by 2^32
}



/*
 * MT_get_bounded_exp() - 
 *
 * Get a number from 0 to upper_bound from an exponential distribution
 *
 */
double
MT_get_bounded_exp() {
  
  double result;
  double p;
  
  
  p = MT_get_real1();
  while(p == 0.0) {
    p = MT_get_real1();
  }

  result = -log(p);

  while(result > MT_EXP_CUTOFF) {

    p = MT_get_real1();
    while(p == 0.0) {
      p = MT_get_real1();
    }
    result = -log(p);
  }
  
  return(result/MT_EXP_CUTOFF);
}



/*
 * MT_print_state() - 
 *
 * Print the state of the random number generatoe
 *
 */
void
MT_print_state(MT_state *state) {
  
  int i;
  
  printf("\n");
  printf("PRNG STATE:\n");
  printf(" mti = %d\n", state->mti);

  for(i=0;i<MT_N;i++) {
    printf("%d ", state->mt[i]);
  }
  printf("\n\n");
  
  return;
}




/*
 * MT_state_cmp() - 
 *
 * Compare two PRNG states, word-by-word
 *
 */
int 
MT_state_cmp(MT_state *state1,
	     MT_state *state2) {

  int i;
  
  if(state1->mti != state2->mti) {
    return(0);  /* not the same */
  }
  
  for(i=0;i<MT_N;i++) {
    if(state1->mt[i] != state2->mt[i]) {
      return(0);  /* not the same */
    }
  }
    
  return(1);  /* the same */
}





/*
 * MT_get_sint32_top() - Returns an int in the range 0..top
 *
 * This function attempts to remove bias in selecting random 
 * integers in a range.
 *
 */
long int
MT_get_sint32_top(long int top) {

  long int overflow;
  long int r;
  long int retval;

  if(top <= 0) {
    return(0); 
  } else {
    overflow = (0x7fffffffUL / top) * top;
  }
  
  while(1) {
    r = MT_get_sint32();
    if(r < overflow) {
      break;
    }
  }
  
  retval = r % top;
  
  return(retval);
}


