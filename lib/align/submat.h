/*
 * submat.h
 *
 * $Id: submat.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2003 Luke Sheneman.  All rights reserved.
 *
 * A set of routines for loading and handling substitution matrices
 * Here, we represent substitution matrices in optimized data structures
 *
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */


#ifndef _INC_SUBMAT_H_
#define _INC_SUBMAT_H_ 1

#include "symtab.h"
#include "align.h"


/* compute array index manually */
#define INDEX(n,x,y)  ((n)*(x)+(y))                              /* a = nsyms,  x&y are indices */
#define VECTOR(s,x,y) (s->vector[INDEX(s->symtab->nsyms, x, y)]) /* s = submat, x&y are indices */
#define VECTOR1(v,n,x,y) (v[INDEX(n,x,y)])



typedef struct ALIGN_submat {
  
  char matrix_type;   /* ALIGN_DNA, ALIGN_DNA_AMBIGUITY, 
			 ALIGN_RNA, ALIGN_PROTEIN */

  /* symbols in the substitution matrix */
  ALIGN_symtab *symtab;
  
  /* the substitution matrix, as a vector */
  double *vector;
  
  /* score categories */
  short int *scores;
  short int nscores;
  
} ALIGN_submat;



/* some function prototypes */
void
ALIGN_print_submat(ALIGN_submat *submat);


int
ALIGN_count_submat_syms(FILE *fp);

int
ALIGN_add_score(int *scorebuf, int score);

short int *
ALIGN_make_scores(int *scorebug, int scorebuf_sz);

void
ALIGN_print_scorebuf(int *scorebuf, int scorebuf_sz);

unsigned int
ALIGN_lookup_score(ALIGN_submat *submat,
		   short int score);

void
ALIGN_print_submat_syms(ALIGN_submat *submat);


void
ALIGN_read_submat_pass1(FILE *fp, ALIGN_submat *submat);

void
ALIGN_read_submat_pass2(FILE *fp, ALIGN_submat *submat);

ALIGN_submat *
ALIGN_read_submat(char *matrix_filename);

void
ALIGN_free_submat(ALIGN_submat *submat);

#endif /* _INC_SUBMAT_H_ */










