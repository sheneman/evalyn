/*
 * align.h
 *
 * $Id: align.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2005 Luke Sheneman.  All rights reserved.
 *
 *
 * A general library for performing dynamic programming alignments for
 * the following situations:
 * 
 *   SEQUENCE  vs. SEQUENCE
 *   SEQUENCE  vs. ALIGNMENT
 *   ALIGNMENT vs. ALIGNMENT
 *
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 * 01/28/03
 *
 */




#ifndef _INC_ALIGN_H_
#define _INC_ALIGN_H_ 1


#include <stdint.h>

#include "syms.h"
#include "symtab.h"
#include "submat.h"


/* some definitions */
#define GAP_CHAR  '-'  /* dec 45 */
#define WILD_CHAR '*'  /* dec 42 */

#define ALIGN_DNA           100
#define ALIGN_DNA_AMBIGUITY 101
#define ALIGN_RNA           102
#define ALIGN_PROTEIN       103

#define ALIGN_SMALL  (-10e300)

#define ALIGN_H     0
#define ALIGN_V     1
#define ALIGN_D     2

#define ALIGN_HBIT  4
#define ALIGN_VBIT  2
#define ALIGN_DBIT  1

#define ALIGN_HV    100
#define ALIGN_HD    101
#define ALIGN_VD    102
#define ALIGN_HVD   103

#define ALIGN_ACTION_ALLOC 100
#define ALIGN_ACTION_FREE  101


/* CACHE DEFINES */
#define ALIGN_HASHSIZE  500
#define ALIGN_CACHESIZE 10000

#define ALIGN_L_PAREN  -100
#define ALIGN_R_PAREN  -101

#define ALIGN_ID_EMPTY -1



/* some limits */
#define MAX_SEQ_TITLE_LENGTH 128    /* remove this limit */
#define MAX_SEQ_LENGTH       10000  /* remove this limit */

#define DEFAULT_PROFSIZE MAX_SYMS+2

/* some macros */
#define ISGAP(a,b) (a->nsyms-1 == b?1:0)  /* a == symtab, b == index */
#define ALIGN_index(a, b, c) ((b) * (a->k) + (c))

#define ALIGN_alignment_length(b) (b->curlen)
#define	MATRIX_cell(a_matrix, a_i, a_j) \
  ((a_matrix)->cells[(a_i) * (a_matrix)->width + (a_j)])

/* macros for indexing alignment profiles */
#define PROFILE(n, c, s) ( ((c)*(n+2)) + (s) )
#define GAPOPEN(a, n, c) ( (a->profile[((c)*(n+2)) + (n)]) )
#define GAPEXT(a, n, c)  ( (a->profile[((c)*(n+2)) + (n+1)]) )
#define GOPENINDEX(n)    ( (n) )
#define GEXTINDEX(n)     ( (n+1) )

#define GETH(c)          ( ( ((c)&4)?1:0) )  /* check for horz */
#define GETV(c)          ( ( ((c)&2)?1:0) )  /* check for vert */
#define GETD(c)          ( ( ((c)&1)?1:0) )  /* check for diag */


/* some datatypes */
typedef struct ALIGN_STRUCT_seq {

  unsigned int length;
  char *data;
  char title[MAX_SEQ_TITLE_LENGTH];

} ALIGN_seq;



typedef struct ALIGN_STRUCT_seqset {

  /* an array of sequences */
  unsigned int num;
  ALIGN_seq *seq;
  
  /* our symbols */
  ALIGN_symtab *symtab;

} ALIGN_seqset;



typedef struct ALIGN_STRUCT_dagnode {

  struct ALIGN_STRUCT_dagnode *horz;
  struct ALIGN_STRUCT_dagnode *vert;
  struct ALIGN_STRUCT_dagnode *diag;
  
  double value;

  unsigned short int npaths;
  
  unsigned int node_id;
  char visited;
  char refcount;

} ALIGN_dagnode;



typedef struct ALIGN_STRUCT_dag {
  
  ALIGN_dagnode *root;

  unsigned int maxulps;

  double score;
  unsigned short int npaths;
  char ties;
  
  unsigned int k;
  short int *id;  
  
} ALIGN_dag;



typedef struct ALIGN_STRUCT_cell {
  
  double value;               /* the DP value                                 */
  char dir;                   /* directions bit-packed as 00000HVD            */
  unsigned short int npaths;  /* number of paths leading to cell              */

  ALIGN_dagnode *dagnode;     /* a dagnode pointer, used in constructing dags */

} ALIGN_cell;


typedef struct ALIGN_STRUCT_matrix {
  
  unsigned short int width;
  unsigned short int height;
  
  unsigned int maxulps;

  ALIGN_cell *cells;

} ALIGN_matrix;


typedef struct ALIGN_STRUCT_alignment {

  unsigned short int n;     /* length */
  unsigned short int k;     /* num of sequences */

  unsigned short int *seq_ids;  
  ALIGN_symtab *symtab;

  char *text;                  /* actual sequences in alignment */
  unsigned short int *profile; /* the alignment profile         */
  
  unsigned short int npaths;   /* the number of paths through the last DP matrix */
  char ties;                   /* flag (0|1) indicating *any* ties               */
  
  short int *id;               /* an array of sequence numbers indicating order of alignment */
  double score;

} ALIGN_alignment;




typedef struct ALIGN_STRUCT_btrace_elm {

  unsigned short int i;
  unsigned short int j;
  
  char dir;  /* ALIGN_D, ALIGN_H, ALIGN_V */

} ALIGN_btrace_elm;



typedef struct ALIGN_STRUCT_btrace {

  unsigned int nelms;
  unsigned short int curlen;

  ALIGN_btrace_elm elms[1];

} ALIGN_btrace;



typedef struct ALIGN_STRUCT_scoring_system {
  
  ALIGN_submat *submat;

  double gap_open;
  double gap_extend;

} ALIGN_scoring_system;


/* need to forward-declare some datatypes */
struct ALIGN_STRUCT_cache_listnode;

typedef struct ALIGN_STRUCT_cache_hashnode {
  
  struct ALIGN_STRUCT_cache_hashnode *next;
  struct ALIGN_STRUCT_cache_hashnode *prev;
  
  struct ALIGN_STRUCT_cache_listnode *listnode;

  ALIGN_alignment *alignment;
  ALIGN_dag *dag;
  
  short int *id;
  unsigned short int n;
  unsigned short int k;
  
  
} ALIGN_cache_hashnode;



typedef struct ALIGN_STRUCT_cache_listnode {
  
  struct ALIGN_STRUCT_cache_listnode *next;
  struct ALIGN_STRUCT_cache_listnode *prev;

  struct  ALIGN_STRUCT_cache_hashnode *hashnode;

} ALIGN_cache_listnode;



typedef struct ALIGN_STRUCT_alignment_cache {
  
  ALIGN_cache_listnode *list;
  ALIGN_cache_hashnode **hash;
  
  ALIGN_cache_listnode *listend;

  int size;

} ALIGN_alignment_cache;


/* llabs is sometimes not defined */
extern long long int llabs(long long int __x);


/* some function prototypes */
ALIGN_alignment *
ALIGN_init_alignment(char *sequence,
		     ALIGN_symtab *symtab,
		     unsigned int length,
		     unsigned short int seq_id);

ALIGN_alignment *
ALIGN_dup_alignment(ALIGN_alignment *source);


ALIGN_matrix *
ALIGN_allocate_matrix(unsigned short int width,
		      unsigned short int height);

void
ALIGN_print_matrix(ALIGN_matrix *matrix);


void
ALIGN_print_backtrace(ALIGN_btrace *btrace);

void
ALIGN_free_alignment(ALIGN_alignment *alignment);

void
ALIGN_free_matrix(void);

ALIGN_alignment *
ALIGN_allocate_alignment(unsigned short int n,
			 unsigned short int k,
			 ALIGN_symtab *symtab);

void
ALIGN_print_alignment(ALIGN_alignment *alignment,
		      ALIGN_seqset *seqset);

void 
ALIGN_print_profile(ALIGN_alignment *alignment);

void
ALIGN_print_seqset(ALIGN_seqset *seqset);

void
ALIGN_free_seqset(ALIGN_seqset *seqset);

int
ALIGN_seqset_type(ALIGN_seqset *seqset);

int 
ALIGN_find_seq(ALIGN_alignment *alignment,
	       int seq_id);

ALIGN_seqset *
ALIGN_alloc_seqset(ALIGN_symtab *symtab);

int
ALIGN_seqset_lookup_sym(ALIGN_seqset *seqset,
			char c);

void
ALIGN_print_profile_column(unsigned short int *prof,
			   int nsyms);


double
ALIGN_compute_sp(ALIGN_scoring_system *ss,
		 ALIGN_alignment *alignment);



ALIGN_alignment *
ALIGN_build_alignment(ALIGN_scoring_system *ss,
		      ALIGN_btrace *btrace,
		      ALIGN_alignment *a,
		      ALIGN_alignment *b);

void
ALIGN_print_matrix_neighborhood(ALIGN_matrix *matrix,
				int i,
				int j);

void
ALIGN_print_matrix_cell(ALIGN_matrix *matrix,
			int i,
			int j);


ALIGN_alignment *
ALIGN_alignment_alignment(ALIGN_scoring_system *ss,
			  ALIGN_alignment *alignment_a, 
			  ALIGN_alignment *alignment_b);

		
#endif /* _INC_ALIGN_H_ */






