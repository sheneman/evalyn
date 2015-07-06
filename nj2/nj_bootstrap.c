/*
 * nj_bootstrap.c
 *
 * $Id: nj_bootstrap.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *****************************************************************************
 *
 * Copyright (c) 2004,  Luke Sheneman
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions 
 * are met:
 * 
 *  + Redistributions of source code must retain the above copyright 
 *    notice, this list of conditions and the following disclaimer. 
 *  + Redistributions in binary form must reproduce the above copyright 
 *    notice, this list of conditions and the following disclaimer in 
 *    the documentation and/or other materials provided with the 
 *    distribution. 
 *  + The names of its contributors may not be used to endorse or promote 
 *    products derived  from this software without specific prior 
 *    written permission. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.  
 *
 *****************************************************************************
 *
 * Functions for bootstrapping
 *
 *****************************************************************************
 *
 * AUTHOR:
 * 
 *   Luke Sheneman
 *   sheneman@cs.uidaho.edu
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include <prng.h>

#include "nj_clearcut.h"
#include "nj_fasta.h"
#include "nj_dmat.h"
#include "nj_dist.h"
#include "nj_bootstrap.h"




/*
 * NJ_bootstrap_alignment() - 
 *
 * Bootstrap a distance matrix by sampling (with replacement) from 
 * columns in a multiple sequence alignment in order to construct a new 
 * bootstrapped distance matrix, which is allocated and returned to the 
 * calling function.
 *
 */
NJ_DMAT *
NJ_bootstrap_alignment(NJ_ARGS *nj_args,
		       NJ_alignment *alignment) {

  NJ_DMAT *dmat;
  NJ_alignment new_alignment;
  unsigned int i, j;
  unsigned int column;
  
  new_alignment.titles = alignment->titles;
  new_alignment.data   = (char *)calloc(alignment->nseq * alignment->length, sizeof(char));
  new_alignment.nseq   = alignment->nseq;
  new_alignment.length = alignment->length;
  
  for(j=0;j<alignment->length;j++) {

    column = MT_get_sint32_top(alignment->length);
    
    for(i=0;i<alignment->nseq;i++) {
      new_alignment.data[i*alignment->length+j] = alignment->data[i*alignment->length+column];
    }
  }
  
  dmat = NJ_compute_dmat(nj_args, &new_alignment);

  free(new_alignment.data);

  return(dmat);
}






/*
 * NJ_build_majority_rule_consensus_tree() - 
 *
 * This function takes an array of trees and constructs a 
 * majority-rule consensus tree.
 *
 */
NJ_TREE *
NJ_build_majority_rule_consensus_tree(NJ_TREE **trees) {

}

