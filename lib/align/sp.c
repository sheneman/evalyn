/*
 * sp.c
 *
 * $Id: sp.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2003 Luke Sheneman.  All rights reserved.
 *
 * Routines for computing sum-of-pairs scores (SPS) given an alignment and 
 * scoring system.
 * 
 *
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */

#include <stdio.h>

#include "align.h"
#include "sp.h"


/*
 * ALIGN_compute_sp() - 
 *
 * Given an alignment and a scoring system, compute the 
 * sum of pairs score
 *
 */
float
ALIGN_compute_sp(ALIGN_scoring_system *ss,
		 ALIGN_alignment *alignment) {

  float score;
  int i;
  
  score = 0.0;
  for(i=0;i<alignment->n;i++) {
    score += 
      ALIGN_score_profile_column(ss,
				 &(alignment->profile[i]),
				 alignment->symtab->nsyms);
  }
  
  return(score);
}


















