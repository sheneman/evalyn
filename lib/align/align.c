/*
 * align.c
 *
 * $Id: align.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2005 Luke Sheneman.  All rights reserved.
 *
 *
 * A general library for performing dynamic programming alignments 
 * using efficient profiles for the following situations:
 * 
 *   SEQUENCE  vs. SEQUENCE
 *   SEQUENCE  vs. ALIGNMENT
 *   ALIGNMENT vs. ALIGNMENT
 *
 * Luke Sheneman
 * 01/28/03
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <values.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#include <prng.h>  /* mersenne twister */

#include "align.h"
#include "cache.h"






/*************************************
 *                                   *
 * BEGIN STATIC FUNCTION PROTOTYPES  *
 *                                   *
 *************************************/


static inline
int
ALIGN_double_ulpcmp(double aA, 
		    double aB, 
		    int64_t aMaxUlps);

static inline
int
ALIGN_double_ulpcmp_eq(double aA, 
		       double aB, 
		       int64_t aMaxUlps);

static inline
void
ALIGN_choose_approximate_match(unsigned int maxulps,
			       int *valid_H,
			       int *valid_V,
			       int *valid_D,
			       double h_diff,
			       double v_diff,
			       double d_diff,
			       int dir);

static inline
void
ALIGN_adjust_affine_gaps(ALIGN_alignment *alignment);


static inline
void
ALIGN_aggregate_profile(ALIGN_scoring_system *ss,
			ALIGN_alignment *dest,
			ALIGN_alignment *source,
			unsigned short int d_col,
			unsigned short int s_col,
			int nsyms);

static inline
void
ALIGN_count_gaptypes(ALIGN_alignment *alignment,
		     int index,
		     int *gopen,
		     int *gext);

static inline
void
ALIGN_HVD_assign(ALIGN_matrix *matrix,
		 int i,
		 int j,
		 double H,
		 double V,
		 double D);

static inline
void
ALIGN_build_profile_column(unsigned short int *prof_c,
			   unsigned short int *a_prof_c,
			   unsigned short int *b_prof_c,
			   int nsyms);

static inline 
double
ALIGN_score_profile_column(ALIGN_scoring_system *ss,
			   unsigned short int *prof_c,
			   int gopen,
			   int gext,
			   int nongaps,
			   int nsyms);

static void
ALIGN_compute_gaps(ALIGN_alignment *a,
		   ALIGN_alignment *b,
		   int matrix_i,
		   int matrix_j,
		   int direction,
		   int prev_direction,
		   int *open_count,
		   int *extend_count,
		   unsigned short int nsyms);

static int
ALIGN_dag_determine_backtrace_move(ALIGN_scoring_system *ss,
				   ALIGN_dag *dag,
				   ALIGN_dagnode *dagnode,
				   ALIGN_alignment *a,
				   ALIGN_alignment *b,
				   int i,
				   int j,
				   int prev,
				   double prev_value);

static int
ALIGN_determine_backtrace_move(ALIGN_scoring_system *ss,
			       ALIGN_matrix *matrix,
			       ALIGN_alignment *a,
			       ALIGN_alignment *b,
			       int i,
			       int j);

static ALIGN_btrace *
ALIGN_backtrace(ALIGN_scoring_system *ss,
		ALIGN_matrix *matrix,
		ALIGN_alignment *a,
		ALIGN_alignment *b);

static ALIGN_btrace *
ALIGN_dag_backtrace(ALIGN_scoring_system *ss,
		    ALIGN_dag *dag,
		    ALIGN_alignment *a,
		    ALIGN_alignment *b);

static void
ALIGN_recurse_dag_backtrace(ALIGN_scoring_system *ss,
			    ALIGN_btrace *btrace,
			    ALIGN_dag *dag,
			    ALIGN_dagnode *dagnode,
			    ALIGN_alignment *a,
			    ALIGN_alignment *b,
			    int i,
			    int j,
			    int prev,
			    double prev_value);
static ALIGN_btrace *
ALIGN_alloc_backtrace(unsigned short int i,
		      unsigned short int j,
		      char action);

static int
ALIGN_break_tie(ALIGN_matrix *matrix,
		ALIGN_dagnode *dagnode,
		int i,
		int j,
		int valid_H,
		int valid_V,
		int valid_D);

static int
ALIGN_valid_check(ALIGN_scoring_system *ss,
		  ALIGN_matrix *matrix,
		  ALIGN_alignment *a,
		  ALIGN_alignment *b,
		  int matrix_i,
		  int matrix_j,
		  int direction,
		  int prev,
		  double *diff);


static int
ALIGN_dag_valid_check(ALIGN_scoring_system *ss,
		      ALIGN_dag *dag,
		      ALIGN_dagnode *dagnode,
		      ALIGN_alignment *a,
		      ALIGN_alignment *b,
		      int direction,
		      int prev,
		      double prev_value,
		      int i,
		      int j,
		      double *diff);


/***********************************
 *                                 *
 * END STATIC FUNCTION PROTOTYPES  *
 *                                 *
 ***********************************/




/***************************************************************/
/***************************************************************/
/***************************************************************/



/**************************
 *                        *
 * BEGIN STATIC FUNCTIONS *
 *                        *
 **************************/




/*
 * ALIGN_double_ulpcmp() - 
 *
 * This function compares two doubles, given a count of 
 * ulps, which is a measure of the overall accumulated error
 * in the contruction of the two doubles.  Evaluating the 
 * equality of floating point numbers only makes sense by
 * accounting for the error in doing such a comparison.
 *
 *
 * More info here:
 *
 * http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
 *
 * Returns:
 *   0 if equal
 *   1 if a > b
 *  -1 if a < b
 *
 */
static inline
int
ALIGN_double_ulpcmp(double aA, 
		    double aB, 
		    int64_t aMaxUlps) {
  int64_t a, b;
  
  // Convert aA and aB to lexicographically ordered ints.
  a = *(int64_t *) &aA;
  if (a < 0) {
    a = 0x8000000000000000LL - a;
  }
  
  b = *(int64_t *) &aB;
  if (b < 0) {
    b = 0x8000000000000000LL - b;
  }
  
  // Check if aA and aB are within aMaxUlps of each other.
  if (llabs(a - b) <= aMaxUlps) {
    return(0);
  } else if (aA > aB) {
    return(1);
  } else {
    return(-1);
  }
}






/*
 * ALIGN_double_ulpcmp_eq() -  
 * 
 * A slightly modified and faster function which checks only for equality
 *
 * Returns: 1 if equal
 *          0 is not equal
 *
 */
static inline
int
ALIGN_double_ulpcmp_eq(double aA, 
		       double aB, 
		       int64_t aMaxUlps) {
  int64_t a, b;

  // Convert aA and aB to lexicographically ordered ints.
  a = *(int64_t *) &aA;
  if (a < 0) {
    a = 0x8000000000000000LL - a;
  }
  
  b = *(int64_t *) &aB;
  if (b < 0) {
    b = 0x8000000000000000LL - b;
  }
  
  // Check if aA and aB are within aMaxUlps of each other.
  if (llabs(a - b) <= aMaxUlps) {
    return(1);
  } else {
    return(0);
  }
}








/*
 * ALIGN_build_profile_column() - 
 *
 * This constructs an aggregate profile column from two
 * profiles.
 *
 */
static inline
void
ALIGN_build_profile_column(unsigned short int *prof_c,
			   unsigned short int *a_prof_c,
			   unsigned short int *b_prof_c,
			   int nsyms) {
  int i;

  /* aggregate the symbols */
  for(i=0;i<nsyms;i++) {
    prof_c[i] =
      a_prof_c[i] +
      b_prof_c[i];
  }

  /* aggregate the gap opens */
  prof_c[nsyms] = 
    a_prof_c[nsyms] + 
    b_prof_c[nsyms];

  /* aggregate the gap extends */
  prof_c[nsyms+1] = 
    a_prof_c[nsyms+1] + 
    b_prof_c[nsyms+1];

  return;
}






/*
 * ALIGN_compute_gaps() - 
 * 
 * This function computes the number of gap opens and
 * gap closes for a cell in the DP matrix for a
 * particular direction in light of the previous directions
 *
 */
static void
ALIGN_compute_gaps(ALIGN_alignment *a,
		   ALIGN_alignment *b,
		   int matrix_i,
		   int matrix_j,
		   int direction,
		   int prev_direction,
		   int *open_count,
		   int *extend_count,
		   unsigned short int nsyms) {
  
  int i, j;
  int gaps, nongaps;
  int gopen, gopen1, gopen2;
  int gext, gext1, gext2;
  
  i = matrix_i - 1;
  j = matrix_j - 1;
  

  switch(direction) {
    
  case ALIGN_H:

    /********************************/
    /*                              */
    /*          HORIZONTAL          */
    /*                              */
    /********************************/
    
    switch(prev_direction) {
    case ALIGN_H:  /* PREV HORIZONTAL */
      /* 
       * Count the number of gap opens and gap extends
       * by looking at adjacent columns in a
       */
      ALIGN_count_gaptypes(a, 
			   j,
			   &gopen,
			   &gext);
      gext += b->k;
      
      break;
      
    case ALIGN_V:  /* PREV_VERTICAL   */

      /* for a, any gaps at this column will be gap extensions */
      gopen = 0;

      //      gext  = a->profile[j].ngap_open + a->profile[j].ngap_ext;
      gext = GAPOPEN(a, nsyms, j) + GAPEXT(a, nsyms, j);
	  
      /* 
       * for b, 1 gap extension for every gap in profile +
       * 1 gap opening for every non-gap in the profile
       */

      //      gaps   = b->profile[i].ngap_open + b->profile[i].ngap_ext;
      gaps = GAPOPEN(b, nsyms, i) + GAPEXT(b, nsyms, i);

      nongaps = (b->k) - gaps;

      gext += gaps;
      gopen = nongaps;

      break;

    case ALIGN_D:  /* PREV_DIAGONAL   */


      ALIGN_count_gaptypes(a, 
			   j,
			   &gopen1,
			   &gext1);
      
      /* 
       * for b, 1 gap extension for every gap in profile +
       * 1 gap opening for every non-gap in the profile
       */
      //      gaps    = b->profile[i].ngap_open + b->profile[i].ngap_ext;
      gaps = GAPOPEN(b, nsyms, i) + GAPEXT(b, nsyms, i);
      nongaps = (b->k) - gaps;
      
      gopen = gopen1 + nongaps;
      gext  = gext1  + gaps;

      break;
      
    default:
      fprintf(stderr, "ALIGN ERROR: invalid prev_direction in ALIGN_compute_gaps()\n");
      exit(-1);
    }

    break;

  case ALIGN_V:

    /********************************/
    /*                              */
    /*          VERTICAL            */
    /*                              */
    /********************************/

    switch(prev_direction) {
    case ALIGN_H:  /* PREV HORIZONTAL */

      /* 
       * for a, 1 gap extension for every gap in profile +
       * 1 gap opening for every non-gap in the profile
       */
      //      gaps    = a->profile[j].ngap_open + a->profile[j].ngap_ext;
      gaps = GAPOPEN(a, nsyms, j) + GAPEXT(a, nsyms, j);
      nongaps = (a->k) - gaps;

      /* for b, any gaps at this column will be gap extensions */
      gopen = 0;
      //      gext  = b->profile[i].ngap_open + b->profile[i].ngap_ext;
      gext = GAPOPEN(b, nsyms, i) + GAPEXT(b, nsyms, i);
	  
      gext  += gaps;
      gopen = nongaps;
      
      break;
      
    case ALIGN_V:  /* PREV_VERTICAL   */

      /* 
       *  Count the number of gap opens and gap extends
       * by looking at adjacent columns in b
       */
      ALIGN_count_gaptypes(b, 
			   i,
			   &gopen,
			   &gext);
      gext += a->k;

      break;

    case ALIGN_D:  /* PREV_DIAGONAL   */

      ALIGN_count_gaptypes(b, 
			   i,
			   &gopen1,
			   &gext1);

      /* 
       * for a, 1 gap extension for every gap in profile +
       * 1 gap opening for every non-gap in the profile
       */
      //     gaps   = a->profile[j].ngap_open + a->profile[j].ngap_ext;
      gaps = GAPOPEN(a, nsyms, j) + GAPEXT(a, nsyms, j);
      nongaps = (a->k) - gaps;

      gopen = gopen1 + nongaps;
      gext  = gext1  + gaps;

      break;
      
    default:
      fprintf(stderr, "ALIGN ERROR: invalid prev_direction in ALIGN_compute_gaps()\n");
      exit(-1);
    }


    break;



  case ALIGN_D:

    /********************************/
    /*                              */
    /*          DIAGONAL            */
    /*                              */
    /********************************/
    
    /* we need to handle the special case where i==1, j==1 */
    if(matrix_i == 1 &&
       matrix_j == 1) {
      
      *open_count =  
	GAPOPEN(a, nsyms, 0) + 
	GAPEXT(a, nsyms, 0)  + 
	GAPOPEN(b, nsyms, 0) + 
	GAPEXT(b, nsyms, 0);

      *extend_count = 0;
      
      return;
    }
    
    /* now handle the general diagonal case */
    switch(prev_direction) {
    case ALIGN_H:  /* PREV HORIZONTAL */
      ALIGN_count_gaptypes(a, 
			   j,
			   &gopen1,
			   &gext1);
	  
      /* for b, any gaps at this column will be gap extensions */
      gopen = 0;
      //      gext = b->profile[i].ngap_open + b->profile[i].ngap_ext;
      gext = GAPOPEN(b, nsyms, i) + GAPEXT(b, nsyms, i);
	  
      gopen = gopen + gopen1;
      gext  = gext  + gext1;
      
      break;
      
    case ALIGN_V:  /* PREV_VERTICAL   */

      /* for a, any gaps at this column will be gap extensions */
      gopen = 0;

      //      gext = a->profile[j].ngap_open + a->profile[j].ngap_ext;
      gext = GAPOPEN(a, nsyms, j) + GAPEXT(a, nsyms, j);
	  
      ALIGN_count_gaptypes(b, 
			   i,
			   &gopen1,
			   &gext1);
	  
      gopen = gopen + gopen1;
      gext  = gext  + gext1;

      break;

    case ALIGN_D:  /* PREV_DIAGONAL   */
      ALIGN_count_gaptypes(a, 
			   j,
			   &gopen1,
			   &gext1);
	  

      ALIGN_count_gaptypes(b, 
			   i,
			   &gopen2,
			   &gext2);
	  
      gopen = gopen1 + gopen2;
      gext  = gext1  + gext2;

      break;
      
    default:
      fprintf(stderr, "ALIGN ERROR: invalid prev_direction in ALIGN_compute_gaps()\n");
      exit(-1);
    }

    break;

  default:
    fprintf(stderr, "ALIGN Error: bogus direction in ALIGN_compute_gaps()\n");
    exit(-1);
  }

  *open_count   = gopen;
  *extend_count = gext;

  return;
}





/* 
 * ALIGN_score_profile_column() - 
 *
 * Given the scoring system and an aggregate profile column,
 * compute the sum-of-pairs score for this cell in the DP 
 * matrix.
 *
 * Now handles affine gap penalties
 *
 */
static inline
double
ALIGN_score_profile_column(ALIGN_scoring_system *ss,
			   unsigned short int *prof_c,
			   int gopen,
			   int gext,
			   int nongaps,
			   int nsyms) {

  double score = 0.0;
  unsigned short int *p, c;
  int i, j, k, increment;
  double *v1;

  v1 = ss->submat->vector;
  increment = nsyms+1;

  /* handle matches */
  p = prof_c;
  for(i=0;i<nsyms;i++) {
    
    c = *p;
    k = 1;
    
    /* all possible matched pairings */
    if(c) {

      /* handle matches */
      score += ((c * (c-1)) >> 1) * (*v1);

      /* handle mismatches */
      for(j=i+1;j<nsyms;j++) {
	if(prof_c[j]) {
	  score += (c * prof_c[j]) * v1[k];
	}
	k++;
      }
    }
    
    p++;
    v1 += increment;
  }

  /* compute the affine gap cost here */
  score += 
    ((double)(nongaps * gopen) * ss->gap_open) +
    ((double)(nongaps * gext)  * ss->gap_extend);
  
  return(score);
}






/*
 * ALIGN_aggregate_profile() - 
 *
 * Increment a profile here
 *
 */
static inline
void
ALIGN_aggregate_profile(ALIGN_scoring_system *ss,
			ALIGN_alignment *dest,
			ALIGN_alignment *source,
			unsigned short int d_col,
			unsigned short int s_col,
			int nsyms) {
  int i;

  for(i=0;i<nsyms;i++) {
    dest->profile[PROFILE(nsyms, d_col, i)] += 
      source->profile[PROFILE(nsyms, s_col, i)];    
  }

  return;
}





/*
 * ALIGN_determine_backtrace_move() - 
 *
 * Given a cell in a matrix, make a decision with respect
 * to which direction to backtrace to.
 *
 * In most cases, there will be only one valid direction.
 * That is, from the current cell in the matrix, and given
 * the value in the current cell, there exists only one 
 * offered path which could produce a cell with such a value.
 * 
 * In the case where there are multiple, viable paths (ties)
 * break the ties in some way. 
 *
 */
static int
ALIGN_determine_backtrace_move(ALIGN_scoring_system *ss,
			       ALIGN_matrix *matrix,
			       ALIGN_alignment *a,
			       ALIGN_alignment *b,
			       int i,
			       int j) {
  int direction;
  static int prev;  /* the previous direction in the backtrace */
  ALIGN_cell *cell;
  int valid_H, valid_V, valid_D;
  double h_diff, v_diff, d_diff;
  
  
  cell = &MATRIX_cell(matrix, i, j);

  /* initialize */
  valid_H = 0;
  valid_V = 0;
  valid_D = 0;

  /* 
   * *** EXTREME LOWER RIGHT ***
   *
   * if this is the lower right cell, any "enabled"
   * directions are equally valid, since we have no 
   * prior backtrace history.  
   *
   * Break remaining ties in an unbiased way.
   *
   */
  if(i == matrix->height - 1 &&
     j == matrix->width  - 1 ) {

    if(GETH(cell->dir)) {
      valid_H = 1;
    }
    if(GETV(cell->dir)) {
      valid_V = 1;
    }
    if(GETD(cell->dir)) {
      valid_D = 1;
    }
    
    direction = 
      ALIGN_break_tie(matrix,
		      NULL,
		      i,
		      j,
		      valid_H,
		      valid_V,
		      valid_D);

    prev = direction;

    return(direction);

  } else {
    
    /* 
     * check to see which of the directions in the cell are valid in the
     * context of the current backtrace
     */

    if(cell->dir == 0) {
      fprintf(stderr, "MATRIX: NO DIRECTIONS FROM THIS CELL (%d,%d)!\n", i, j);
      ALIGN_print_matrix(matrix);
      exit(-1);
    }
    
    /* HORIZONTAL */
    if(GETH(cell->dir)) {
      valid_H =  ALIGN_valid_check(ss, matrix, a, b, i, j, ALIGN_H, prev, &h_diff);
    } 
    
    /* VERTICAL   */
    if(GETV(cell->dir)) {
      valid_V =  ALIGN_valid_check(ss, matrix, a, b, i, j, ALIGN_V, prev, &v_diff);
    }

    /* DIAGONAL   */
    if(GETD(cell->dir)) {
      valid_D =  ALIGN_valid_check(ss, matrix, a, b, i, j, ALIGN_D, prev, &d_diff);
    }

    if(!valid_D &&
       !valid_H &&
       !valid_V) {
      
      /* none of the directions are exact matches so assign at least one which is the closest match */
      ALIGN_choose_approximate_match(matrix->maxulps, &valid_H, &valid_V, &valid_D, h_diff, v_diff, d_diff, cell->dir);
    }

    /* choose a direction here in a fair, unbiased way */
    direction = 
      ALIGN_break_tie(matrix,
		      NULL,
		      i,
		      j,
		      valid_H,
		      valid_V,
		      valid_D);

    /* remember our direction, as it becomes our previous direction */
    prev = direction;
  }

  return(direction);
}





/*
 * ALIGN_dag_determine_backtrace_move() - 
 * 
 * Determine the next backtrace move, using the DAG
 *
 */
static int
ALIGN_dag_determine_backtrace_move(ALIGN_scoring_system *ss,
				   ALIGN_dag *dag,
				   ALIGN_dagnode *dagnode,
				   ALIGN_alignment *a,
				   ALIGN_alignment *b,
				   int i,
				   int j,
				   int prev,
				   double prev_value) {
  
  int direction;
  int valid_H=0, valid_V=0, valid_D=0;
  double h_diff, v_diff, d_diff;
  int adir = 0;
  
  
  /* in the last node, any 
   * "enabled" direction is equally 
   * valid, since we have no prior
   * backtrace history.
   *
   * Break remaining ties in an unbiased way
   *
   */
  if(dagnode == dag->root) {
    
    if(dagnode->horz) {
      valid_H = 1;
    }
    
    if(dagnode->vert) {
      valid_V = 1;
    }
    
    if(dagnode->diag) {
      valid_D = 1;
    }
    
    direction =
      ALIGN_break_tie(NULL,
		      dagnode,
		      0,
		      0,
		      valid_H,
		      valid_V,
		      valid_D);
    
    return(direction);

  } else {
    
    if(!dagnode->horz &&
       !dagnode->vert &&
       !dagnode->diag) {

      fprintf(stderr, "DAG: NO DIRECTIONS FROM THIS NODE\n");
      exit(-1);
    }
    
    if(dagnode->horz) {
      adir += ALIGN_HBIT;
      valid_H = ALIGN_dag_valid_check(ss, dag, dagnode, a, b, ALIGN_H, prev, prev_value, i, j, &h_diff);
    }
    
    if(dagnode->vert) {
      adir += ALIGN_VBIT;
      valid_V = ALIGN_dag_valid_check(ss, dag, dagnode, a, b, ALIGN_V, prev, prev_value, i, j, &v_diff);
    }

    if(dagnode->diag) {
      adir += ALIGN_DBIT;
      valid_D = ALIGN_dag_valid_check(ss, dag, dagnode, a, b, ALIGN_D, prev, prev_value, i, j, &d_diff);
    }

    if(!valid_D &&
       !valid_H &&
       !valid_V) {

      /* none of the directions are exact matches so assign at least one which is the closest match */
      ALIGN_choose_approximate_match(dag->maxulps, &valid_H, &valid_V, &valid_D, h_diff, v_diff, d_diff, adir);
    }

    /* choose a direction here in a fair, unbiased way */
    direction = 
      ALIGN_break_tie(NULL,
		      dagnode,
		      0,
		      0,
		      valid_H,
		      valid_V,
		      valid_D);
    
    /* remember our direction, as it becomes our previous direction */
    prev = direction;
    prev_value = dagnode->value;

  }

  return(direction);
}







/*
 * ALIGN_valid_check() - 
 *
 * As part of the backtrace, we need to check to see if
 * a particular cell is part of the backtrace pattern
 * The cell is in the specified matrix at the
 * matrix_i and matrix_j coordinates.
 *
 * We are checking if the cell at this point can
 * tie together previous steps in the backtrace (prev)
 * to the next step in the backtrace (direction).
 *
 * Basically, we see if we can compute the value
 * in the previous cell from the current cell by 
 * following the specified direction in the backtrace
 * Thus, it takes 3 cells in the matrix to confirm
 * a backtrace move
 *
 */
static int
ALIGN_valid_check(ALIGN_scoring_system *ss,
		  ALIGN_matrix *matrix,
		  ALIGN_alignment *a,
		  ALIGN_alignment *b,
		  int matrix_i,
		  int matrix_j,
		  int direction,
		  int prev,
		  double *diff) {

  ALIGN_cell *c_cell, *p_cell;  
  unsigned short int prof_c[DEFAULT_PROFSIZE];
  int nongaps;

  int result=0, i, j;
  double prev_value, T, score;
  int gopen, gext, nsyms;

  /* the top left corner is always valid from these cells */
  if( (matrix_i == 0 && matrix_j == 1) ||
      (matrix_i == 1 && matrix_j == 0) ) {
    return(1);
  }

  nsyms  = ss->submat->symtab->nsyms;
  
  *diff = 0.0;

  i = matrix_i - 1;
  j = matrix_j - 1;
  
  c_cell  = &(MATRIX_cell(matrix, matrix_i, matrix_j));

  switch(prev) {
  case ALIGN_H:

    p_cell = &(MATRIX_cell(matrix, matrix_i, matrix_j+1));
    prev_value = p_cell->value;

    ALIGN_compute_gaps(a,
		       b,
		       matrix_i,
		       matrix_j+1,
		       ALIGN_H,
		       direction,
		       &gopen,
		       &gext,
		       nsyms);

    nongaps = a->k + b->k - gopen - gext;
	
    score = 
      ALIGN_score_profile_column(ss,
				 &(a->profile[PROFILE(nsyms, matrix_j, 0)]),
				 gopen,
				 gext,
				 nongaps,
				 nsyms);
    
    T = MATRIX_cell(matrix, matrix_i, matrix_j).value + score;
    if(ALIGN_double_ulpcmp_eq(T, prev_value, matrix->maxulps)) {
      result = 1;
    } else {
      *diff = fabs(T - prev_value);
      result = 0;
    }


    break;

  case ALIGN_V:

    /* see if gap open vs. a explains it */

    /* determine the value in the previous cell in the backtrace */
    p_cell = &(MATRIX_cell(matrix, matrix_i+1, matrix_j));
    prev_value = p_cell->value;

    /* compute the gap opening and extensions */
    ALIGN_compute_gaps(a, 
		       b, 
		       matrix_i+1, 
		       matrix_j, 
		       ALIGN_V,   
		       direction,
		       &gopen,
		       &gext,
		       nsyms);

    nongaps = a->k + b->k - gopen - gext;

    score = 
      ALIGN_score_profile_column(ss,
				 &(b->profile[PROFILE(nsyms, matrix_i, 0)]),
				 gopen,
				 gext,
				 nongaps,
				 nsyms);

    T = MATRIX_cell(matrix, matrix_i, matrix_j).value + score;
    if(ALIGN_double_ulpcmp_eq(T, prev_value, matrix->maxulps)) {
      result = 1;
    } else {
      *diff = fabs(T - prev_value);
      result = 0;
    }

    break;

  case ALIGN_D:

    p_cell = &(MATRIX_cell(matrix, matrix_i+1, matrix_j+1));
    prev_value = p_cell->value;
    
    ALIGN_compute_gaps(a,
		       b,
		       matrix_i+1,
		       matrix_j+1,
		       ALIGN_D,
		       direction,
		       &gopen,
		       &gext,
		       nsyms);
    
    ALIGN_build_profile_column(prof_c,
			       &(a->profile[PROFILE(nsyms, matrix_j, 0)]),
			       &(b->profile[PROFILE(nsyms ,matrix_i, 0)]),
			       nsyms);
    nongaps = a->k + b->k - gopen - gext;

    score = 
      ALIGN_score_profile_column(ss,
				 prof_c,
				 gopen,
				 gext,
				 nongaps,
				 nsyms);
    
    T = MATRIX_cell(matrix, matrix_i, matrix_j).value + score;
    if(ALIGN_double_ulpcmp_eq(T, prev_value, matrix->maxulps)) {
      result = 1;
    } else {
      *diff = fabs(T - prev_value);
      result = 0;
    }

    break;

  default:
    fprintf(stderr, "ALIGN ERROR: Incorrect prev value in ALIGN_valid_V()\n");
    exit(-1);
  }

  return(result);
}




/*
 * ALIGN_backtrace() - 
 *
 * Perform backtrace operation in a given DP matrix
 *
 */
static ALIGN_btrace *
ALIGN_backtrace(ALIGN_scoring_system *ss,
		ALIGN_matrix *matrix,
		ALIGN_alignment *a,
		ALIGN_alignment *b) {
  
  ALIGN_btrace *btrace;
  unsigned int i, j;
  int direction;

  /* start at the lower right corner of the matrix */
  i = b->n;
  j = a->n;
  
  /* the end of the list */
  btrace = ALIGN_alloc_backtrace(i, j, ALIGN_ACTION_ALLOC);
  btrace->elms[0].i = i;
  btrace->elms[0].j = j;

  /* 
     here we want to traverse from the lower right corner of
     the matrix up to the upper left corner of the matrix.  At
     each cell, there are three directions (horz, vert, diag).
     if more than one of these pointers is "on", then this 
     indicates that there was more than one path that leads
     to this cell in the DP matrix (a tie).  The tie must be 
     broken.  In most cases, only one direction makes sense 
     in the context of where we came from so far in the
     backtracking, so we choose the direction which makes sense.
     However, it may be possible that multiple tied directions
     are equally viable.  In this case, we must break the tie.
     currently, lets do that by randomly choosing one of the
     viable paths with proportional probability.
  */

  while(i || j) {

    /* find a valid direction to backtrace, breaking ties as necessary */
    direction = 
      ALIGN_determine_backtrace_move(ss,
				     matrix, 
				     a,
				     b,
				     i, 
				     j);

    
    switch(direction) {
    case ALIGN_D:
      i?i--:0;
      j?j--:0;
      break;
    case ALIGN_H:
      j?j--:0;
      break;
    case ALIGN_V:
      i?i--:0;
      break;
    default:
      printf("Should never be here.  ALIGN_backtrace()\n");
      exit(-1);
    }

    btrace->elms[btrace->curlen].dir = direction;
    btrace->curlen++;
    btrace->elms[btrace->curlen].i = i;
    btrace->elms[btrace->curlen].j = j;
    
  }
  
  return(btrace);
}



/*
 * ALIGN_count_gaptypes() - 
 *
 * This function takes the textual portion of the alignment and 
 * counts the number of gap openings and gap extensions if
 * we were to add the column at index to the column at 
 * index-1
 *
 */
static inline
void
ALIGN_count_gaptypes(ALIGN_alignment *alignment,
		     int index,
		     int *gopen,
		     int *gext) {
  
  int i, k;
  int open = 0, extend = 0;
  char *p1, *p2;


  k = alignment->k;
  
  if(index == 0) {
    for(i=0;i<k;i++) {
      if(alignment->text[i] == GAP_CHAR) {
	open++;
      }
    }
  } else {

    p1 = &(alignment->text[ALIGN_index(alignment, index-1, 0)]);
    p2 = &(alignment->text[ALIGN_index(alignment, index,   0)]);
    
    for(i=0;i<k;i++) {

      if(*p1 == GAP_CHAR) {
	if(*p2 == GAP_CHAR) {
	  extend++;
	}
      } else {
	if(*p2 == GAP_CHAR) {
	  open++;
	}
      }

      p1++;
      p2++;

    }
  }
  
  *gopen = open;
  *gext  = extend;

  return;
}




/*
 * ALIGN_adjust_affine_gaps() - 
 *
 * Given an alignment, assign affine gap costs to the 
 * alignment profile
 *
 */
static inline
void
ALIGN_adjust_affine_gaps(ALIGN_alignment *alignment) {

  int i, j, k, n, nsyms;
  char *text, *t1, *t2;
  unsigned short int *profile, *p1, *p2;


  /* some shortcuts to limit dereferencing */
  nsyms   = alignment->symtab->nsyms;
  text    = alignment->text;
  profile = alignment->profile;
  k       = alignment->k;
  n       = alignment->n;

  /* 
   * handle the first column specially, since it
   * starts the alignment, and do it outside
   * the main loop for speed reasons
   */
  profile[GOPENINDEX(nsyms)] = 0;
  profile[GEXTINDEX(nsyms)]  = 0;
  
  t1 = text;
  for(i=0;i<k;i++) {
    if(*t1++ == GAP_CHAR) {
      profile[GOPENINDEX(nsyms)]++;
    }
  }
  
  /* main loop */
  t1 = &(text[k]);
  t2 = text;
  p1 = &(profile[PROFILE(nsyms, 1, GOPENINDEX(nsyms))]);
  p2 = &(profile[PROFILE(nsyms, 1, GEXTINDEX(nsyms))]);
  
  for(j=1;j<n;j++) {

    *p1 = 0;  // gap open
    *p2 = 0;  // gap extend
      
    for(i=0;i<k;i++) {
      if(*t1 == GAP_CHAR) {
	if(*t2 == GAP_CHAR) {
	  (*p2)++;
	} else {
	  (*p1)++;
	}
      }
      
      t1++;
      t2++;
    }
    p1 += nsyms+2;
    p2 += nsyms+2;
  }
  
  return;
}





/*
 * ALIGN_HVD_assign() 
 *
 * Here, we assign directions according to the values H, V, and D
 * We break ties as necessary
 *
 */
static inline
void
ALIGN_HVD_assign(ALIGN_matrix *matrix,
		 int i,
		 int j,
		 double H,
		 double V,
		 double D) {

  double X, Y, p;
  int dir1, dir2;
  int tie = 0;


  /* TOURNAMENT */

  /* determine X */
  if(ALIGN_double_ulpcmp_eq(H, V, matrix->maxulps)) {     // H == V

    tie = 1;
    
    p = MT_get_real1();
    if(p < 0.5) {
      X = H;
    } else {
      X = V;
    }

    dir1 = ALIGN_HV;

  } else {

    if(ALIGN_double_ulpcmp(H, V, matrix->maxulps) > 0) {   // H > V
      X = H;
      dir1 = ALIGN_H;
    } else {
      X = V;
      dir1 = ALIGN_V;
    }

  }

  /* determine Y */
  if(ALIGN_double_ulpcmp_eq(X, D, matrix->maxulps)) {      // X == D

    p = MT_get_real1();

    if(tie) {
      if(p < 0.333) {
	Y = D;
      } else {
	Y = X;
      }
    } else {
      if(p < 0.5) {
	Y = D;
      } else {
	Y = X;
      }
    }

    if(dir1 == ALIGN_H) {
      dir2 = ALIGN_HD;
    } else if(dir1 == ALIGN_V) {
      dir2 = ALIGN_VD;
    } else {
      dir2 = ALIGN_HVD;      
    }

  } else {

    if(ALIGN_double_ulpcmp(X, D, matrix->maxulps) > 0) {
      Y = X;
      dir2 = dir1;
    } else {
      Y = D;
      dir2 = ALIGN_D;
    }

  }

  MATRIX_cell(matrix, i, j).value = Y;

  /* initialize */
  MATRIX_cell(matrix, i, j).dir = 0;
  
  switch(dir2) {

  case ALIGN_H:

    MATRIX_cell(matrix, i, j).dir = ALIGN_HBIT;  // horz
    MATRIX_cell(matrix, i, j).npaths = MATRIX_cell(matrix, i, j-1).npaths;

    break;

  case ALIGN_V:

    MATRIX_cell(matrix, i, j).dir = ALIGN_VBIT;  // vert
    MATRIX_cell(matrix, i, j).npaths = MATRIX_cell(matrix, i-1, j).npaths;

    break;

  case ALIGN_D:

    MATRIX_cell(matrix, i, j).dir = ALIGN_DBIT;  // diag
    MATRIX_cell(matrix, i, j).npaths = MATRIX_cell(matrix, i-1, j-1).npaths;

    break;

  case ALIGN_HV:  

    MATRIX_cell(matrix, i, j).dir = ALIGN_HBIT + ALIGN_VBIT;  // horz + vert

    MATRIX_cell(matrix, i, j).npaths = 
      MATRIX_cell(matrix, i, j-1).npaths +
      MATRIX_cell(matrix, i-1, j).npaths;

    break;

  case ALIGN_HD: 

    MATRIX_cell(matrix, i, j).dir = ALIGN_HBIT + ALIGN_DBIT;  // horz + diag

    MATRIX_cell(matrix, i, j).npaths = 
      MATRIX_cell(matrix, i, j-1).npaths +
      MATRIX_cell(matrix, i-1, j-1).npaths;

    break;

  case ALIGN_VD: 

    MATRIX_cell(matrix, i, j).dir = ALIGN_VBIT + ALIGN_DBIT;  // vert + diag

    MATRIX_cell(matrix, i, j).npaths = 
      MATRIX_cell(matrix, i-1, j).npaths +
      MATRIX_cell(matrix, i-1, j-1).npaths;

    break;

  case ALIGN_HVD:

    MATRIX_cell(matrix, i, j).dir = 
      ALIGN_HBIT + 
      ALIGN_VBIT + 
      ALIGN_DBIT;  // horz + vert + diag

    MATRIX_cell(matrix, i, j).npaths = 
      MATRIX_cell(matrix, i-1, j).npaths +
      MATRIX_cell(matrix, i, j-1).npaths +
      MATRIX_cell(matrix, i-1, j-1).npaths;

    break;

  }

  return;
}
  



/*
 * ALIGN_break_tie() - 
 *
 * Break ties here between valid directions during
 * backtracing
 *
 */
static int
ALIGN_break_tie(ALIGN_matrix *matrix,
		ALIGN_dagnode *dagnode,
		int i,
		int j, 
		int valid_H,
		int valid_V,
		int valid_D) {

  int H_npaths=0, V_npaths=0, D_npaths=0;
  double H_ratio, V_ratio, D_ratio;
  double p;
  
  /* check to see that there are valid paths */
  if(!valid_H && !valid_V && !valid_D) {
    fprintf(stderr, "FATAL ASSERTION: ALIGN ERROR: No valid paths in ALIGN_break_tie()\n");
    exit(-1);
  }

  if(matrix) {
    if(valid_H) {
      H_npaths = MATRIX_cell(matrix, i, j-1).npaths;
    }
    if(valid_V) {
      V_npaths = MATRIX_cell(matrix, i-1, j).npaths;
    }
    if(valid_D) {
      D_npaths = MATRIX_cell(matrix, i-1, j-1).npaths;
    }
  } else if(dagnode) {
    if(valid_H) {
      H_npaths = dagnode->horz->npaths;
    }
    if(valid_V) {
      V_npaths = dagnode->vert->npaths;
    }
    if(valid_D) {
      D_npaths = dagnode->diag->npaths;
    }
  } else {
    fprintf(stderr, "Error: invalid matrix AND dagnode in ALIGN_break_tie().  Exiting.");
    exit(-1);
  }


  if(valid_H && !valid_V && !valid_D) {        /* only H is valid */
    return(ALIGN_H);
  } else if(valid_V && !valid_H && !valid_D) { /* only V is valid */
    return(ALIGN_V);
  } else if(valid_D && !valid_H && !valid_V) { /* only D is valid */
    return(ALIGN_D);
  } else if(valid_H && valid_D && !valid_V) {  /* tie between H and D */
    
    H_ratio = (double)H_npaths/(double)(H_npaths + D_npaths);
    
    p = MT_get_real1();

    // printf("tie-1\n");

    if(p<H_ratio) {
      return(ALIGN_H);
    } else {
      return(ALIGN_D);
    }
  } else if(valid_V && valid_D && !valid_H) {  /* tie between V and D */

    V_ratio = (double)V_npaths/(double)(V_npaths + D_npaths);

    p = MT_get_real1();

    // printf("tie-2\n");

    if(p<V_ratio) {
      return(ALIGN_V);
    } else {
      return(ALIGN_D);
    }
  } else if(valid_H && valid_V && !valid_D) {  /* tie between H and V */

    H_ratio = (double)H_npaths/(double)(H_npaths + V_npaths);

    p = MT_get_real1();

    // printf("tie-3\n");

    if(p<H_ratio) {
      return(ALIGN_H);
    } else {
      return(ALIGN_V);
    }
  } else if(valid_H && valid_V && valid_D) {  /* 3-way tie between H, V, and D */

    H_ratio = (double)H_npaths/(double)(H_npaths + V_npaths + D_npaths);
    V_ratio = (double)V_npaths/(double)(H_npaths + V_npaths + D_npaths);
    D_ratio = (double)D_npaths/(double)(H_npaths + V_npaths + D_npaths);

    p = MT_get_real1();

    // printf("tie-4\n");

    /* total pie = < H_ratio...V_ratio...D_ratio > */

    if(p<H_ratio) {
      return(ALIGN_H);
    } else if(p >= H_ratio && p < (1.0-D_ratio) ) {
      return(ALIGN_V);
    } else {
      return(ALIGN_D);
    }
  } else {
    printf("ALIGN ERROR: Assertion in ALIGN_break_tie()\n");
    printf("valid_D = %d, valid_H = %d, valid_V = %d\n", valid_D, valid_H, valid_V);
    exit(-1);
  }
}





/*
 * ALIGN_alloc_backtrace() - 
 *
 * Allocate a backtrace element here
 *
 */
static ALIGN_btrace *
ALIGN_alloc_backtrace(unsigned short int i, 
		      unsigned short int j,
		      char action) {

  static ALIGN_btrace *retval = NULL;
  unsigned short int nelms = i + j + 1;
  long int size;
  
  switch(action) {

  case ALIGN_ACTION_FREE:
    if(retval) {
      free(retval);
    }
    
    return(NULL);
    
    break;

  case ALIGN_ACTION_ALLOC:

    if (retval != NULL) {
      if (retval->nelms < nelms) {
	
	size = sizeof(ALIGN_btrace) + (nelms - 1) * sizeof(ALIGN_btrace_elm);
	retval = (ALIGN_btrace *)realloc(retval, size);
	
	retval->nelms = nelms;
      }
      retval->curlen = 0;

    } else {

      size = sizeof(ALIGN_btrace) + (nelms - 1) * sizeof(ALIGN_btrace_elm);
      retval = (ALIGN_btrace *)malloc(size);

      retval->nelms = nelms;
      retval->curlen = 0;

    }
    
    return retval;

    break;

  default:
    fprintf(stderr, "ALIGN: Unknown Action in ALIGN_alloc_backtrace\n");
    exit(-1);
    
  }

  fprintf(stderr, "ALIGN: Unknown Action in ALIGN_alloc_backtrace\n");
  exit(-1);
}




/*
 * ALIGN_choose_approximate_match() - 
 *
 * This function is used when no exact match could be found, so 
 * we want to choose the direction leading to the closest possible match
 * in a way which is maximally unbiased
 *
 * If there are multiple minima, all minima can be returned
 *
 */
static void
ALIGN_choose_approximate_match(unsigned int maxulps,
			       int *valid_H, 
			       int *valid_V, 
			       int *valid_D, 
			       double h_diff, 
			       double v_diff, 
			       double d_diff,
			       int dir) {

  double X, p;
  int res;

  *valid_H = 0;
  *valid_V = 0;
  *valid_D = 0;
  
  if(dir == ALIGN_HBIT) {
    *valid_H = 1;
  } else if(dir == ALIGN_VBIT) {
    *valid_V = 1;
  } else if(dir == ALIGN_DBIT) {
    *valid_D = 1;
  } else if(dir == (ALIGN_HBIT + ALIGN_VBIT)) {

    res = ALIGN_double_ulpcmp(h_diff, v_diff, maxulps);
    if(res == 0) {
      *valid_H = 1;
      *valid_V = 1;
    } else if(res < 0) {
      *valid_H = 1;
    } else {
      *valid_V = 1;
    }

  } else if(dir == (ALIGN_HBIT + ALIGN_DBIT)) {

    res = ALIGN_double_ulpcmp(h_diff, d_diff, maxulps);
    if(res == 0) {
      *valid_H = 1;
      *valid_D = 1;
    } else if(res < 0) {
      *valid_H = 1;
    } else {
      *valid_D = 1;
    }

  } else if(dir == (ALIGN_VBIT + ALIGN_DBIT)) {

    res = ALIGN_double_ulpcmp(v_diff, d_diff, maxulps);
    if(res == 0) {
      *valid_V = 1;
      *valid_D = 1;
    } else if(res < 0) {
      *valid_V = 1;
    } else {
      *valid_D = 1;
    }

  } else if(dir == (ALIGN_HBIT + ALIGN_VBIT + ALIGN_DBIT)) {

    /* do a fair tournament to break 3-way */
    /* first, determine X */
    if(GETH(dir) && GETV(dir) && GETD(dir)) {
      if(ALIGN_double_ulpcmp_eq(h_diff, v_diff, maxulps)) {    /* h_diff == v_diff */

	*valid_H = 1;
	*valid_V = 1;

	p = MT_get_real1();
      
	if(p<0.5) {
	  X = h_diff;
	} else {
	  X = v_diff;
	}
      } else {
	if(ALIGN_double_ulpcmp(h_diff, v_diff, maxulps) < 0) {  /* h_diff < v_diff */

	  X = h_diff;
	
	  *valid_H = 1;
      
	} else {

	  X = v_diff;

	  *valid_V = 1;
	}
      }

  
      if(ALIGN_double_ulpcmp_eq(X, d_diff, maxulps)) {          /* X == d_diff */
   
	p = MT_get_real1();

	*valid_D = 1;

      } else {

	if(ALIGN_double_ulpcmp(d_diff, X, maxulps) < 0) {

	  *valid_H = 0;
	  *valid_V = 0;
	  *valid_D = 1;

	}
      }
    }
       
  } else {
    /* assertion, no directions */
    printf("ALIGN_choose_approximate_match() -- no directions from current cell\n");
    exit(-1);
  }

  return;
}





/**************************
 *                        *
 *  END STATIC FUNCTIONS  *
 *                        *
 **************************/



/****************************************************************/
/****************************************************************/
/****************************************************************/



/****************************
 *                          *
 *  BEGIN PUBLIC FUNCTIONS  *
 *                          *
 ****************************/



/*
 * ALIGN_print_alignment() -
 *
 * Print an alignment
 *
 */
void
ALIGN_print_alignment(ALIGN_alignment *alignment,
		      ALIGN_seqset *seqset) {
  
  unsigned short int i, j;

  printf("ALIGNMENT LENGTH: %d\n", alignment->n);
  printf("ALIGNMENT # SEQS: %d\n", alignment->k);
  printf("ALIGNMENT NPATHS: %d\n", alignment->npaths);
  printf("ALIGNMENT TIES:   %d\n", alignment->ties);
  printf("ALIGNMENT ID:     ");
  ALIGN_print_alignment_id(alignment->id);
  printf("\n");
  
  for(i=0;i<alignment->k;i++) {
    if(seqset) {
      printf("%s\t", seqset->seq[alignment->seq_ids[i]].title);
    } else {
      printf("seq_id: %d\t", alignment->seq_ids[i]);
    }
    for(j=0;j<alignment->n;j++) {
      if(alignment->text[ALIGN_index(alignment, j, i)] == GAP_CHAR) {
	printf("%c", GAP_CHAR);
      } else {
	printf("%c", alignment->symtab->syms[(int)(alignment->text[ALIGN_index(alignment, j, i)])]);
      }
    }
    printf("\n");
  }
  
  printf("\n\n");
}






/*
 * ALIGN_free_alignment() - 
 *
 * Free the space allocated for an alignment 
 *
 */
void
ALIGN_free_alignment(ALIGN_alignment *alignment) {
  
  free(alignment->text);
  free(alignment->profile);
  free(alignment->seq_ids);
  free(alignment->id);

  free(alignment);

  return;
}








/*
 * ALIGN_dag_backtrace() - 
 * 
 * Perform a backtrace using the DAG instead of the matrix 
 *
 */
ALIGN_btrace *
ALIGN_dag_backtrace(ALIGN_scoring_system *ss,
		    ALIGN_dag *dag,
		    ALIGN_alignment *a,
		    ALIGN_alignment *b) {

  ALIGN_btrace *btrace;
  unsigned int i, j;
  
  i = b->n;
  j = a->n;
  
  btrace = ALIGN_alloc_backtrace(i, j, ALIGN_ACTION_ALLOC);
  btrace->elms[0].i = i;
  btrace->elms[0].j = j;
  
  ALIGN_recurse_dag_backtrace(ss,
			      btrace,
			      dag,
			      dag->root,
			      a,
			      b,
			      i,
			      j,
			      0,
			      0);
  return(btrace);
}




/*
 * ALIGN_recurse_dag_backtrace() -
 *
 * Recursively backtrace
 *
 */
void
ALIGN_recurse_dag_backtrace(ALIGN_scoring_system *ss,
			    ALIGN_btrace *btrace,
			    ALIGN_dag *dag,
			    ALIGN_dagnode *dagnode,
			    ALIGN_alignment *a,
			    ALIGN_alignment *b,
			    int i,
			    int j,
			    int prev,
			    double prev_value) {
  
  int direction;
  
  /* check to see if we have hit the end of the recursion */
  if( (i == 0) && (j == 0) ) {
    return;
  }

  direction = 
    ALIGN_dag_determine_backtrace_move(ss,
				       dag,
				       dagnode,
				       a,
				       b,
				       i,
				       j,
				       prev,
				       prev_value);

  switch(direction) {
  case ALIGN_H:
    j?j--:0;
    break;
  case ALIGN_V:
    i?i--:0;
    break;
  case ALIGN_D:
    i?i--:0;
    j?j--:0;
    break;
  default:
    printf("Should never be here.  ALIGN_recurse_dag_backtrace()\n");
    exit(-1);
  }

  /* add to the backtrace here */
  btrace->elms[btrace->curlen].dir = direction;
  btrace->curlen++;
  btrace->elms[btrace->curlen].i = i;
  btrace->elms[btrace->curlen].j = j;
  
  /* recurse */
  if(direction == ALIGN_H) {
    ALIGN_recurse_dag_backtrace(ss, btrace, dag, dagnode->horz, a, b, i, j, ALIGN_H, dagnode->value);
  } else if(direction == ALIGN_V) {
    ALIGN_recurse_dag_backtrace(ss, btrace, dag, dagnode->vert, a, b, i, j, ALIGN_V, dagnode->value);
  } else if(direction == ALIGN_D) {
    ALIGN_recurse_dag_backtrace(ss, btrace, dag, dagnode->diag, a, b, i, j, ALIGN_D, dagnode->value);
  } else {
    fprintf(stderr, "Error:  invalid direction in ALIGN_recurse_dag_backtrace()\n");
    exit(-1);
  }

  return;
}




/*
 * ALIGN_print_backtrace() - 
 *
 * Print the contents of a backtrace
 *
 */
void
ALIGN_print_backtrace(ALIGN_btrace *btrace) {

  short int i;

  printf("Backtrace path.  Length = %d: \n", btrace->curlen);

  for(i=0;i<btrace->curlen;i++) {
    
    switch (btrace->elms[i].dir) {
    case ALIGN_H:
      printf("H");
      break;
    case ALIGN_V:
      printf("V");
      break;
    case ALIGN_D:
      printf("D");
      break;
    default:
      printf("XXX\n");
    }
  }

  printf("\nPath through DP matrix:\n");

  for (i = btrace->curlen - 1; i >= 0; i--) {

    switch (btrace->elms[i].dir) {
    case ALIGN_H:
      printf("H");
      break;
    case ALIGN_V:
      printf("V");
      break;
    case ALIGN_D:
      printf("D");
      break;
    default:
      printf("XXX\n");
    }
  }

  printf("\n\n");
 
  return;
}




/*
 * ALIGN_allocate_alignment() -
 *
 * Allocate an alignment here 
 *
 */
ALIGN_alignment *
ALIGN_allocate_alignment(unsigned short int n,
			 unsigned short int k,
			 ALIGN_symtab *symtab) {

  ALIGN_alignment *alignment;

  
  /* allocate the alignment */
  alignment = (ALIGN_alignment *)calloc(1, sizeof(ALIGN_alignment));
  if(!alignment) {
    fprintf(stderr, "In ALIGN_allocate_alignment() - Failed to allocate alignment\n");
    return(NULL);
  }
  
  /* allocate alignment "text" */
  alignment->text = (char *)calloc(n*k, sizeof(char));
  
  /* allocate the profile here */
  //alignment->profile = (ALIGN_profile *)calloc((n+1), sizeof(ALIGN_profile));
  alignment->profile = 
    (unsigned short int *)calloc((symtab->nsyms+2)*(n+1), sizeof(unsigned short int));

  /* allocate seq_ids array here */
  alignment->seq_ids = (unsigned short int *)malloc(k*sizeof(unsigned short int));
  
  /* allocate the alignment ID string here */
  alignment->id = (short int *)malloc(k*3*sizeof(short int));
  memset(alignment->id, ALIGN_ID_EMPTY, k*3*sizeof(short int));
  
  /* associate alignment symtab with provided symtab */
  alignment->symtab = symtab;

  /* set the other parameters here */
  alignment->n = n;  /* length */
  alignment->k = k;  /* number */
  
  alignment->npaths = 0;   /* the number of paths in the associated DP matrix */
  alignment->ties   = 0;   /* whether or not this alignment has *any* ties    */
  
  alignment->score  = ALIGN_SMALL;

  return(alignment);
}




/*
 * ALIGN_allocate_matrix() -
 *
 * Dynamically allocate a dynamic programming matrix of specified size
 *
 */
ALIGN_matrix *
ALIGN_allocate_matrix(unsigned short int width,
		      unsigned short int height) {

  static ALIGN_matrix matrix;
  static unsigned int ncells = 0;
  unsigned int need;
  int i, j;


  /* the matrix has an extra first row and first column */
  matrix.width  = width  + 1;
  matrix.height = height + 1;

  need = matrix.width * matrix.height;

  if(ncells >= need) {

    memset(matrix.cells, 0, need * sizeof(ALIGN_cell));
    return(&matrix);

  } else if(ncells != 0) {

    ncells = need << 1;

    /* reallocate cells. */
    matrix.cells =
      (ALIGN_cell *)realloc(matrix.cells, ncells*sizeof(ALIGN_cell));

  } else {

    ncells = need;

    /* allocate cells */
    matrix.cells = (ALIGN_cell *)malloc(ncells*sizeof(ALIGN_cell));
  }

  /* zero the matrix */
  memset(matrix.cells, 0, ncells*sizeof(ALIGN_cell));
  
  /* initialize the matrix direction signs */
  for(i=0;i<matrix.height;i++) {
    for(j=0;j<matrix.width;j++) {

      MATRIX_cell(&matrix, i, j).dir = 0;

      /* clear npaths, too */
      MATRIX_cell(&matrix, i, j).npaths = 0;
    }
  }
    
  return(&matrix);
}







/*
 * ALIGN_print_matrix() -
 *
 * Pretty-print the specified dynamic programming matrix
 *
 */
void
ALIGN_print_matrix(ALIGN_matrix *matrix) {
  
  unsigned short int i, j;

  printf("\n");

  printf("  ");
  //  for (j = 0; j < matrix->width; j++) {
  for (j = 0; j < 4; j++) {
    printf("%17d ", j);
  }
  printf("\n");

  //  for(i=0;i<matrix->height;i++) {
  for(i=0;i<4;i++) {
    printf("%5d ", i);
    //    for(j=0;j<matrix->width;j++) {
    for(j=0;j<4;j++) {
      ALIGN_print_matrix_cell(matrix, i, j);
    }
    printf("\n");
  }
  
  printf("\n");
  
  return;
}









/*
 * ALIGN_build_alignment() -
 *
 * Given a backtrace, construct an alignment 
 *
 */
ALIGN_alignment *
ALIGN_build_alignment(ALIGN_scoring_system *ss,
		      ALIGN_btrace *btrace,
		      ALIGN_alignment *a,
		      ALIGN_alignment *b) {
  
  ALIGN_alignment *result;
  unsigned int alignment_length;
  int j, h, pos1, pos2, nsyms;
  

  nsyms = ss->submat->symtab->nsyms;

  /* construct an alignment from backtrace */
  alignment_length = ALIGN_alignment_length(btrace);
  result = 
    ALIGN_allocate_alignment(alignment_length, 
			     (a->k + b->k),
			     ss->submat->symtab);

  /* copy seq ids here */
  memcpy(result->seq_ids, a->seq_ids, (a->k)*sizeof(unsigned short int));
  memcpy(&(result->seq_ids[a->k]), b->seq_ids, (b->k)*sizeof(unsigned short int));

  /* traverse backtrace and construct alignment and profile */
  j = pos1 = pos2 = 0;
  for (h=btrace->curlen-1;h>=0;h--) {
      
    switch(btrace->elms[h].dir) {
    case ALIGN_V:
      memset( &(result->text[ALIGN_index(result, j, 0)]), 
	      GAP_CHAR, 
	      a->k );

      memcpy( &(result->text[ALIGN_index(result, j, a->k)]),
	      &(b->text[ALIGN_index(b, pos2, 0)]),
	      b->k );

      ALIGN_aggregate_profile(ss, result, b, j, pos2, nsyms); 

      pos2++;

      break;

    case ALIGN_H:

      memcpy( &(result->text[ALIGN_index(result, j, 0)]),
	      &(a->text[ALIGN_index(a, pos1, 0)]),
	      a->k );

      memset( &(result->text[ALIGN_index(result, j, a->k)]),
	      GAP_CHAR,
	      b->k );
      
      ALIGN_aggregate_profile(ss, result, a, j, pos1, nsyms); 
      pos1++;

      break;

    case ALIGN_D:

      memcpy( &(result->text[ALIGN_index(result, j, 0)]),
	      &(a->text[ALIGN_index(a, pos1, 0)]),
	      a->k );
	      
      memcpy( &(result->text[ALIGN_index(result, j, a->k)]),
	      &(b->text[ALIGN_index(b, pos2, 0)]),
	      b->k );
	

      ALIGN_aggregate_profile(ss, result, a, j, pos1, nsyms); 
      ALIGN_aggregate_profile(ss, result, b, j, pos2, nsyms); 
      pos1++;
      pos2++;

      break;

    default:
      fprintf(stderr, "UNKNOWN DIR IN BACKTRACE\n");
      exit(-1);
    }

    j++;
  }
  
  /* 
   * At this point, the alignment has been constructed, and 
   * the profile has been updated with the exception of 
   * affine gap foo.  We need to do one pass through the
   * "textual" alignment to properly assign gaps to the 
   * the profile.
   */
  ALIGN_adjust_affine_gaps(result);

  /* compute and assign the alignment id here */
  if(result->id) {
    free(result->id);
  }
  result->id = ALIGN_compute_alignment_id(a, b);
  
  return(result);
}




/*
 * ALIGN_print_profile() -
 *
 * Print a profile for the given alignment
 */
void
ALIGN_print_profile(ALIGN_alignment *alignment) {
  
  unsigned short int i, j;
  int nsyms;
  
  nsyms = alignment->symtab->nsyms;
  
  printf("In ALIGN_print_profile()\n");
  printf("Number of symbols in profile: %d\n", alignment->symtab->nsyms);

  /* iterate through all symbols in profile */
  for(i=0;i<alignment->symtab->nsyms;i++) {

    printf("[%c]\t", alignment->symtab->syms[i]);

    /* iterate through all columns in the alignment */
    for(j=0;j<alignment->n;j++) {
      printf("%d ", alignment->profile[PROFILE(nsyms, j, i)]);
    }

    printf("\n");

  }
  
  /* print the gap openings */
  printf("[-]\t");
  for(j=0;j<alignment->n;j++) {
    printf("%d ", GAPOPEN(alignment, nsyms, j));
  }
  printf("\n");

  /* print the gap extensions */
  printf("[.]\t");
  for(j=0;j<alignment->n;j++) {
    printf("%d ", GAPEXT(alignment, nsyms, j));
  }
  printf("\n\n");
  
  return;
}






/*
 * ALIGN_init_alignment() -
 *
 * Initialize an alignment here
 *
 */
ALIGN_alignment *
ALIGN_init_alignment(char *sequence,
		     ALIGN_symtab *symtab,
		     unsigned int length,
		     unsigned short int seq_id) {

  ALIGN_alignment *alignment;
  unsigned short int i, nsyms;


  /* some shortcuts to limit dereferencing */
  nsyms = symtab->nsyms;
 
  /* allocate alignment and initialize content and seq_ids */
  alignment = ALIGN_allocate_alignment(length, 1, symtab);
  alignment->seq_ids[0] = seq_id;
  alignment->npaths = 1;
  alignment->ties = 0;
  alignment->score = ALIGN_SMALL;

  alignment->id[0] = seq_id;

  /* init alignment and profile */
  for(i=0;i<length;i++) {
    alignment->text[ALIGN_index(alignment, i, 0)] = sequence[i];
    alignment->profile[PROFILE(nsyms, i, sequence[i])]++;
  }

  return(alignment);
}





/*
 * ALIGN_print_seqset() - 
 *
 * Print a set of sequences
 *
 */
void 
ALIGN_print_seqset(ALIGN_seqset *seqset) {

  unsigned int i, j;

  printf("\n\n");
  
  for(i=0;i<seqset->num;i++) {
    printf("[%s]\n", seqset->seq[i].title);
    for(j=0;j<seqset->seq[i].length;j++) {
      /* printf("%d ", seqset->seq[i].data[j]); */  /* the index */
      printf("%c", seqset->symtab->syms[(int)(seqset->seq[i].data[j])]);
    }
    printf("\n");
  }
  
  return;
}





/*
 * ALIGN_seqset_type() - 
 * 
 * Example seqset and automatically determine if these sequences are
 * protein or DNA
 *
 * DNA input sequences can include ambiguity codes, so be sure to
 * handle those
 *
 * returns:
 *   ALIGN_DNA 
 *   ALIGN_PROTEIN
 *   -1 (if things are ambiguous)
 *
 */
int
ALIGN_seqset_type(ALIGN_seqset *seqset) {
  
  int i, flag;


  /* check to see if this is DNA */
  flag = 1;
  for(i=0;i<seqset->num;i++) {
    if(!(ALIGN_SEQ_is_DNA(seqset->seq[i].data))) {
      flag = 0;
      break;
    }
  }
  
  if(flag) {
    return(ALIGN_DNA);
  }



  /* check to see if this is dna with ambiguity symbols */
  flag = 1;
  for(i=0;i<seqset->num;i++) {
    if(!(ALIGN_SEQ_is_DNA_AMBIGUITY(seqset->seq[i].data))) {
      flag = 0;
      break;
    }
  }
  
  if(flag) {
    /*    return(ALIGN_DNA_AMBIGUITY);*/
    printf("ALIGN_DNA_AMBIGUITY\n");
    return(-1);    /* for now */
  }


  /* check to see if this is rna  */
  flag = 1;
  for(i=0;i<seqset->num;i++) {
    if(!(ALIGN_SEQ_is_RNA(seqset->seq[i].data))) {
      flag = 0;
      break;
    }
  }
  
  if(flag) {
    /*    return(ALIGN_RNA); */
    printf("ALIGN_RNA\n");
    return(-1);   /* for now */
  }


  /* finally, check to see if it the sequences are protein */
  flag = 1;
  for(i=0;i<seqset->num;i++) {
    if(!(ALIGN_SEQ_is_PROTEIN(seqset->seq[i].data))) {
      flag = 0;
      printf("flag = 0, i=%d\n", i);
      break;
    }
  }
  
  if(flag) {
    return(ALIGN_PROTEIN);
  }
  
  /* WOAH.  Bad.  This isn't anything we can handle */
  printf("BADNESS\n");
  return(-1);
}




/*
 * ALIGN_free_matrix() -
 * 
 * Free the DP matrix and backtrace data structure
 *
 */
void
ALIGN_free_matrix(void) {
  
  ALIGN_matrix *matrix;
  
  /* free the matrix here */
  matrix = ALIGN_allocate_matrix(0,0);
  if(matrix) {
    if(matrix->cells) {
      free(matrix->cells);
    }
  }
  
  /* free the backtrace */
  ALIGN_alloc_backtrace(0, 0, ALIGN_ACTION_FREE);

  return;
}




/*
 * ALIGN_free_seqset() - 
 *
 * Free the seqset here
 *
 */
void 
ALIGN_free_seqset(ALIGN_seqset *seqset) {

  unsigned short int i;
  
  for(i=0;i<seqset->num;i++) {
    if(seqset->seq && seqset->seq[i].data) {
      free(seqset->seq[i].data);
    }
  }

  if(seqset->seq) {
    free(seqset->seq);
  }
  
  if(seqset) {
    free(seqset);
  }
  
  return;
}




/*
 * ALIGN_find_seq() - 
 *
 * Given an alignment and a seq_id, return the conceptual row
 * of the given alignment which corresponds to the seq_id
 *
 * RETURNS:
 * --------
 *
 *  -1 if specified seq_id is not in the alignment
 *   otherwise, return index in alignment to this sequence
 *
 */
int
ALIGN_find_seq(ALIGN_alignment *alignment,
	       int seq_id) {

  int i;
  
  for(i=0;i<alignment->k;i++) {
    if(alignment->seq_ids[i] == seq_id) {
      return(i);
    }
  }
  
  return(-1);
}




/*
 * ALIGN_alloc_seqset() - 
 *
 * Dynamically allocate a seqset and initialize it using
 * the symbol table used in the scoring system submat
 *
 */
ALIGN_seqset *
ALIGN_alloc_seqset(ALIGN_symtab *symtab) {
  
  ALIGN_seqset *seqset; 

  seqset = (ALIGN_seqset *)calloc(1, sizeof(ALIGN_seqset));
  seqset->symtab = symtab;

  return(seqset);
}





/*
 * ALIGN_seqset_lookup_sym()
 *
 * Lookup the index of a particular symbol
 *
 * RETURNS:
 * --------
 *
 *  -1 if the symbol is not preent
 *  otherwise, the index to the symbol
 *
 */
int
ALIGN_seqset_lookup_sym(ALIGN_seqset *seqset,
			char c) {
  
  int i;
  
  for(i=0;i<seqset->symtab->nsyms;i++) {
    if(seqset->symtab->syms[i] == toupper(c)) {
      return(i);
    }
  }
  
  return(-1);
}









/*
 * ALIGN_print_profile_column() - 
 *
 * Print a profile column
 *
 */
void
ALIGN_print_profile_column(unsigned short int *prof,
			   int nsyms) {
  
  int i;
  
  for(i=0;i<nsyms;i++) {
    printf("%d] %d\n", i, prof[i]);
  }
  
  printf("-] %d\n", prof[GOPENINDEX(nsyms)]);
  printf(".] %d\n", prof[GEXTINDEX(nsyms)]);
  
  return;
}
















/*
 * ALIGN_dag_valid_check() - 
 *
 * As part of the backtrace, we need to check to see if
 * a particular dagnode is part of the backtrace pattern

 * We are checking if the node at this point can
 * tie together previous steps in the backtrace (prev)
 * to the next step in the backtrace (direction).
 *
 * Basically, we see if we can compute the value
 * in the previous node from the current node by 
 * following the specified direction in the backtrace
 * Thus, it takes 3 nodes in the DAG to confirm
 * a backtrace move
 *
 */
int
ALIGN_dag_valid_check(ALIGN_scoring_system *ss,
		      ALIGN_dag *dag,
		      ALIGN_dagnode *dagnode,
		      ALIGN_alignment *a,
		      ALIGN_alignment *b,
		      int direction,
		      int prev,
		      double prev_value,
		      int i,
		      int j,
		      double *diff) {
  int result = 0;
  int nsyms, nongaps;
  int gopen, gext;
  double T, score;
  unsigned short int prof_c[DEFAULT_PROFSIZE];
  
  nsyms = ss->submat->symtab->nsyms;
  
  *diff = 0.0;
  
  switch(prev) {
  case ALIGN_H:

    ALIGN_compute_gaps(a,
		       b,
		       i,
		       j+1,
		       ALIGN_H,
		       direction,
		       &gopen,
		       &gext,
		       nsyms);

    nongaps = a->k + b->k - gopen - gext;
		
    score = 
      ALIGN_score_profile_column(ss,
				 &(a->profile[PROFILE(nsyms, j, 0)]),
				 gopen,
				 gext,
				 nongaps,
				 nsyms);
    T = dagnode->value + score;

    if(ALIGN_double_ulpcmp_eq(T, prev_value, dag->maxulps)) {
      result = 1;
    } else {
      *diff = fabs(T - prev_value);
      result = 0;
    }

    break;

  case ALIGN_V:

    ALIGN_compute_gaps(a,
		       b,
		       i+1,
		       j,
		       ALIGN_V,
		       direction,
		       &gopen,
		       &gext,
		       nsyms);
    
    nongaps = a->k + b->k - gopen - gext;
    
    score = 
      ALIGN_score_profile_column(ss,
				 &(b->profile[PROFILE(nsyms, i, 0)]),
				 gopen,
				 gext,
				 nongaps,
				 nsyms);
    T = dagnode->value + score;

    if(ALIGN_double_ulpcmp_eq(T, prev_value, dag->maxulps)) {
      result = 1;
    } else {
      *diff = fabs(T - prev_value);
      result = 0;
    }

    break;

  case ALIGN_D:

    ALIGN_compute_gaps(a,
		       b,
		       i+1,
		       j+1,
		       ALIGN_D,
		       direction,
		       &gopen,
		       &gext,
		       nsyms);
    
    ALIGN_build_profile_column(prof_c,
			       &(a->profile[PROFILE(nsyms, j, 0)]),
			       &(b->profile[PROFILE(nsyms, i, 0)]),
			       nsyms);
    nongaps = a->k + b->k - gopen - gext;
    
    score = 
      ALIGN_score_profile_column(ss,
				 prof_c,
				 gopen,
				 gext,
				 nongaps,
				 nsyms);

    T = dagnode->value + score;

    if(ALIGN_double_ulpcmp_eq(T, prev_value, dag->maxulps)) {
      result = 1;
    } else {
      *diff = fabs(T - prev_value);
      result = 0;
    }
    
    break;

  default:
    printf("fuckored\n");
    exit(-1);
    break;
   
  }
  
  return(result);
}






/* 
 * ALIGN_dup_alignment() - 
 *
 * Allocate a new alignment and copy the specified alignment into it 
 *
 */
ALIGN_alignment *
ALIGN_dup_alignment(ALIGN_alignment *source) {
  
  ALIGN_alignment *dest;
  int n, k, nsyms;
  
  nsyms = source->symtab->nsyms;
  n = source->n;
  k = source->k;
  
  dest = ALIGN_allocate_alignment(n, k, source->symtab);

  memcpy(dest->text,    source->text,    n*k*sizeof(char));
  memcpy(dest->profile, source->profile, (nsyms+2)*(n+1)*sizeof(unsigned short int));
  memcpy(dest->seq_ids, source->seq_ids, k*sizeof(unsigned short int));
  memcpy(dest->id,      source->id,      (k*3)*sizeof(short int));
  
  dest->npaths = source->npaths;
  dest->ties   = source->ties;
  dest->score  = source->score;
  
  return(dest);
}





/*
 * ALIGN_compute_sp() - 
 *
 * Given an alignment and a scoring system, compute the 
 * sum of pairs score
 *
 */
double
ALIGN_compute_sp(ALIGN_scoring_system *ss,
		 ALIGN_alignment *alignment) {

  double score;
  int i;
  int gopen, gext, nongaps;
  int nsyms;
  
  nsyms = alignment->symtab->nsyms;

  score = 0.0;
  for(i=0;i<alignment->n;i++) {

    ALIGN_count_gaptypes(alignment,
			 i,
			 &gopen,
			 &gext);

    nongaps = alignment->k - gopen - gext;

    score += 
      ALIGN_score_profile_column(ss,
				 &(alignment->profile[PROFILE(nsyms,i, 0)]),
				 gopen,
				 gext,
				 nongaps,
				 nsyms);
  }
  
  return(score);
}






/*
 * ALIGN_print_matrix_neighborhood() - 
 *
 * Print all cells neighboring the cell at i,j
 * Used for debugging purposes.
 * 
 */
void
ALIGN_print_matrix_neighborhood(ALIGN_matrix *matrix,
				int i,
				int j) {
  
  printf("ALIGN_print_matrix_neighborhood()\n");
  printf("  Printing neighbors of i=%d, j=%d\n", i, j);
  
  ALIGN_print_matrix_cell(matrix, i-1, j-1);
  ALIGN_print_matrix_cell(matrix, i-1, j);
  ALIGN_print_matrix_cell(matrix, i-1, j+1);
  printf("\n");

  ALIGN_print_matrix_cell(matrix, i,   j-1);
  ALIGN_print_matrix_cell(matrix, i,   j);
  ALIGN_print_matrix_cell(matrix, i,   j+1);
  printf("\n");

  ALIGN_print_matrix_cell(matrix, i+1, j-1);
  ALIGN_print_matrix_cell(matrix, i+1, j);
  ALIGN_print_matrix_cell(matrix, i+1, j+1);
  printf("\n");
  
  return;
}



/*
 * ALIGN_print_matrix_cell() - 
 *
 * Print the contents of the matrix at position i, j
 *
 */
void 
ALIGN_print_matrix_cell(ALIGN_matrix *matrix,
			int i,
			int j) {
  
  /* do some simple bounds checking here */
  if(i < 0 || j < 0 || i >= matrix->height || j >= matrix->width) {
    printf("              <END>");
    return;
  }

  printf("%12.20lf", MATRIX_cell(matrix, i, j).value);
  //  printf("%18Lx", MATRIX_cell(matrix, i, j).value);
      
  printf("|");

  if(GETH(MATRIX_cell(matrix, i, j).dir)) {
    printf("H");
  } else {
    printf(" ");
  }

  if(GETV(MATRIX_cell(matrix, i, j).dir)) {
    printf("V");
  } else {
    printf(" ");
  }

  if(GETD(MATRIX_cell(matrix, i, j).dir)) {
    printf("D");
  } else {
    printf(" ");
  }
  printf("|");

  printf("%d", MATRIX_cell(matrix, i, j).npaths);

  return;
}







/* 
 * ALIGN_alignment_alignment() - 
 *
 * Align an alignment to an alignment
 *
 */
ALIGN_alignment *
ALIGN_alignment_alignment(ALIGN_scoring_system *ss,
			  ALIGN_alignment *a,
			  ALIGN_alignment *b) {
  
  ALIGN_alignment *result;
  ALIGN_alignment *c;
  ALIGN_matrix *matrix;
  ALIGN_btrace *btrace;
  ALIGN_submat *submat;
  ALIGN_dag *dag;

  char pdir;
  unsigned short int i, j, prev, nsyms;
  int gopen, gext, nongaps, newsize;
  double score, H, V, D, T;

  /* scratch profile column */
  unsigned short int prof_c[DEFAULT_PROFSIZE];

  /* the alignment cache */
  static ALIGN_alignment_cache *cache = NULL;
  ALIGN_cache_hashnode *hashnode;
  short int *id;
  ALIGN_alignment *acopy;
  
  /* we need to align things in a canonical order to pass all regression tests */
  if(ALIGN_determine_id_order(a->id, b->id) > 0) {
    c = a;
    a = b;
    b = c;
  }

  /* 
   * the new alignment will have the number of sequences equal to the sum of the 
   * number of sequences of its two sub-alignments
   */
  newsize = a->k + b->k;

  /*
  printf("ALIGNING THESE:\n");
  ALIGN_print_alignment(a, NULL);
  ALIGN_print_alignment(b, NULL);
  fflush(stdout);
  */
  
  /* if this is the first time through, allocate the alignment cache */
  if(cache == NULL) {
    cache = ALIGN_init_alignment_cache();
  }
  
  //  ALIGN_print_alignment_cache(cache);
  
  /* check to see if this alignment already exists in the alignment cache */

  id = ALIGN_compute_alignment_id(a, b);
  hashnode = ALIGN_find_hashnode(cache, id, (a->k)+(b->k));
  free(id);

  if(hashnode) {
    
    ALIGN_cache_promote(cache, hashnode, hashnode->listnode);
    
    /*
    if(hashnode->alignment && hashnode->dag) {
      printf("DOOM1\n");
      exit(-1);
    }
    
    if(!hashnode->alignment && !hashnode->dag) {
      printf("DOOM2\n");
      exit(-1);
    }
    */

    if(hashnode->alignment) {

      acopy = ALIGN_dup_alignment(hashnode->alignment);

      return(acopy);

    } else if(hashnode->dag) {
      
      btrace = ALIGN_dag_backtrace(ss, hashnode->dag, a, b);
      result = ALIGN_build_alignment(ss, btrace, a,  b);

      result->score  = hashnode->dag->score;
      result->ties   = hashnode->dag->ties;
      result->npaths = hashnode->dag->npaths;
      
      return(result);

    } else {
      printf("FATAL ASSERTION: HASHNODE WITH NEITHER AN ALIGNMENT NOR A DAG in ALIGN_alignment_alignment()\n");
      exit(-1);
    }
  } 


  /* some "aliases" for speed (less dereferencing) and code readability */
  nsyms  = ss->submat->symtab->nsyms;
  submat = ss->submat;

  /* allocate dp matrix */
  matrix = ALIGN_allocate_matrix(a->n, b->n);
  //matrix->maxulps = 36*nsyms*nsyms+10;  /* the maxulps used in float comparison */
  matrix->maxulps = 4*nsyms*nsyms+10;
  
  // ALIGN_print_matrix(matrix);
  
  /* initialize first (upper left) cell */
  MATRIX_cell(matrix, 0, 0).value  = 0.0;
  MATRIX_cell(matrix, 0, 0).npaths = 1;

  MATRIX_cell(matrix, 0, 0).dir = ALIGN_HBIT + ALIGN_VBIT + ALIGN_DBIT;  /* all directions */

  
  /***********************************
    
       MATRIX FILL CODE
       
  ***********************************/

  /* FILL IN THE TOP ROW HERE */
  for(j=1;j<matrix->width;j++) {
    if(j==1) {
      gopen = (b->k) + GAPOPEN(a, nsyms, j-1) + GAPEXT(a, nsyms, j-1);
      gext  = 0;
    } else {
      gopen = GAPOPEN(a, nsyms, j-1);
      gext  = GAPEXT(a,  nsyms, j-1) + (b->k);
    }

    nongaps = newsize - gopen - gext;

    MATRIX_cell(matrix, 0, j).dir = ALIGN_HBIT;   /* set to horz */
    MATRIX_cell(matrix, 0, j).npaths = 1;

    H = ALIGN_score_profile_column(ss,
				   &(a->profile[PROFILE(nsyms, j-1, 0)]),
				   gopen,
				   gext,
				   nongaps,
				   nsyms);
    
    MATRIX_cell(matrix, 0, j).value = 
      MATRIX_cell(matrix, 0, j-1).value + H;
  }
  

  /* FILL IN THE LEFTMOST COLUMN */
  for(i=1;i<matrix->height;i++) {
    if(i==1) {
      gopen = (a->k) + GAPOPEN(b, nsyms, i-1) + GAPEXT(b, nsyms, i-1);
      gext  = 0;

    } else {
      gopen = GAPOPEN(b, nsyms, i-1);
      gext  = GAPEXT(b, nsyms, i-1) + (a->k);

    }

    nongaps = newsize - gopen - gext;

    MATRIX_cell(matrix, i, 0).dir = ALIGN_VBIT;  /* set to vert */
    MATRIX_cell(matrix, i, 0).npaths = 1;
    
    V = ALIGN_score_profile_column(ss,
				   &(b->profile[PROFILE(nsyms, i-1, 0)]),
				   gopen,
				   gext,
				   nongaps,
				   nsyms);
    
    MATRIX_cell(matrix, i, 0).value = 
      MATRIX_cell(matrix, i-1, 0).value + V;
  }




  /* FILL IN THE REMAINDER OF THE MATRIX */
  for(i=1;i<matrix->height;i++) {
    for(j=1;j<matrix->width;j++) {

      /* HORIZONTAL */
      /* Here, we need to compute what the value of a cell is if we
	 go horizontal.  However, the horizontal (left) cell might have been 
	 derived through ties.  If so, we need to determine the value 
	 of the current cell in light of all of the ties which might have
	 been used to derive the value in the previous (horizontal/left) cell.  
	 We will use the value which ultimately gives us the best (largest)
	 value for the current cell.
      */
      H = ALIGN_SMALL;
      pdir = MATRIX_cell(matrix, i, j-1).dir;
      for(prev=ALIGN_H;prev<=ALIGN_D;prev++) {

	if( (prev == ALIGN_H && GETH(pdir)) ||
	    (prev == ALIGN_V && GETV(pdir)) ||
	    (prev == ALIGN_D && GETD(pdir)) ) {

	  ALIGN_compute_gaps(a, b, i, j, ALIGN_H, prev, &gopen, &gext, nsyms);
	  nongaps = newsize - gopen - gext;

	  score = 
	    ALIGN_score_profile_column(ss,
				       &(a->profile[PROFILE(nsyms, j-1, 0)]),
				       gopen,
				       gext,
				       nongaps,
				       nsyms);

	  T = MATRIX_cell(matrix, i, j-1).value + score;
	  if(ALIGN_double_ulpcmp(T, H, matrix->maxulps) > 0) {
	    H = T;
	  }
	}
      }


      /* VERTICAL */
      /* Here, we need to compute what the value of a cell is if we
	 go vertical.  However, the vertical (up) cell might have been 
	 derived through ties.  If so, we need to determine the value 
	 of the current cell in light of all of the ties which might have
	 been used to derive the value in the previous (vertical/up) cell.  
	 We will use the value which ultimately gives us the best (largest)
	 value for the current cell.
      */
      V = ALIGN_SMALL;
      pdir = MATRIX_cell(matrix, i-1, j).dir;
      for(prev=ALIGN_H;prev<=ALIGN_D;prev++) {

	if( (prev == ALIGN_H && GETH(pdir)) ||
	    (prev == ALIGN_V && GETV(pdir)) ||
	    (prev == ALIGN_D && GETD(pdir)) ) {


	  ALIGN_compute_gaps(a, b, i, j, ALIGN_V, prev, &gopen, &gext, nsyms);
	  nongaps = newsize - gopen - gext;

	  score = 
	    ALIGN_score_profile_column(ss,
				       &(b->profile[PROFILE(nsyms, i-1, 0)]),
				       gopen,
				       gext,
				       nongaps,
				       nsyms);

	  T = MATRIX_cell(matrix, i-1, j).value + score;
	  if(ALIGN_double_ulpcmp(T, V, matrix->maxulps) > 0) {
	    V = T;
	  }
	}
      }

      /* DIAGONAL */
      /* Here, we need to compute what the value of a cell is if we
	 go diagonal.  However, the diagonal cell might have been 
	 derived through ties.  If so, we need to determine the value 
	 of the current cell in light of all of the ties which might have
	 been used to derive the value in the previous (diagonal) cell.  
	 We will use the value which ultimately gives us the best (largest)
	 value for the current cell.
      */
      D = ALIGN_SMALL;
      pdir = MATRIX_cell(matrix, i-1, j-1).dir;
      for(prev=ALIGN_H;prev<=ALIGN_D;prev++) {

	if( (prev == ALIGN_H && GETH(pdir)) ||
	    (prev == ALIGN_V && GETV(pdir)) ||
	    (prev == ALIGN_D && GETD(pdir)) ) {

	  ALIGN_compute_gaps(a, b, i, j, ALIGN_D, prev, &gopen, &gext, nsyms);
	
	  ALIGN_build_profile_column(prof_c,
				     &(a->profile[PROFILE(nsyms, j-1, 0)]),
				     &(b->profile[PROFILE(nsyms, i-1, 0)]),
				     nsyms);
	  nongaps = newsize - gopen - gext;

	  score = 
	    ALIGN_score_profile_column(ss,
				       prof_c,
				       gopen,
				       gext,
				       nongaps,
				       nsyms);

	  T = MATRIX_cell(matrix, i-1, j-1).value + score;
	  if(ALIGN_double_ulpcmp(T, D, matrix->maxulps) > 0) {
	    D = T;
	  }
	}
      }

      /* assign directions and break ties as necessary */
      ALIGN_HVD_assign(matrix,
		       i,
		       j,
		       H,
		       V,
		       D);

    }
  }
  
  /* perform the backtrace */
  btrace = ALIGN_backtrace(ss, matrix, a, b);
  
  //  printf("BACKTRACE:\n");
  //  ALIGN_print_backtrace(btrace);

  /* actually build the alignment here */
  result = 
    ALIGN_build_alignment(ss,
			  btrace,
			  a,
			  b);

  /* record the value in the lower-right corner of the matrix */
  result->score = 
    MATRIX_cell(matrix, 
		matrix->height-1, 
		matrix->width-1).value;

  result->npaths = 
    MATRIX_cell(matrix,
		matrix->height-1,
		matrix->width-1).npaths;
  
  if(result->npaths > 1) {
    result->ties = 1;
  }
  

  /* if an alignment is built of ties, it is itself a tied alignment */
  if( a->ties      == 0 &&
      b->ties      == 0 &&
      result->ties == 0 ) {
    result->ties = 0;
  } else {
    result->ties = 1;
  }

  /* construct a DAG if we need to */
  if(a->ties == 0 &&
     b->ties == 0 &&
     result->ties == 1) {

    dag = ALIGN_build_dag(matrix, a, b);

    dag->ties   = result->ties;
    dag->npaths = result->npaths;
    dag->score  = result->score;
  } else {
    dag = NULL;
  }

  /* Add both the alignment and the dag to the cache
   *
   * NOTE that the alignment will only really be added to the cache if
   * the alignment has no ties.  The DAG is only really added if
   * this DP matrix has ties, but previous alignments did not.
   */
  if(result->ties == 0) {
    ALIGN_add_cache_node(cache, result, NULL, result->n, result->k, result->id);
  } else if(((result->ties == 1) && (a->ties == 0) && (b->ties == 0))) {
    ALIGN_add_cache_node(cache, NULL, dag, result->n, result->k, result->id);
  }

  /*
  ALIGN_print_alignment(result, NULL);
  printf("LR SCORE: %lf\n", result->score);
  printf("SP SCORE: %lf\n", ALIGN_compute_sp(ss, result));
  ALIGN_print_profile(result);
  ALIGN_print_matrix(matrix);
  */

  return(result);
}




/***************************************************************/
/***************************************************************/
/***************************************************************/



/************************
 *                      *
 * END PUBLIC FUNCTIONS *
 *                      *
 ************************/

