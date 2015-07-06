/*
 * fasta.h 
 *
 * $Id: fasta.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2003 Luke Sheneman.  All rights reserved.
 *
 * A set of routines for handling the FASTA file format
 *
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */

#ifndef _INC_FASTA_H_
#define _INC_FASTA_H_ 1

#include "align.h"

/* the states of the parser state machine */
#define ALIGN_FASTA_UNKNOWN_STATE  0
#define ALIGN_FASTA_BRACKET_STATE  1
#define ALIGN_FASTA_TITLE_STATE    2
#define ALIGN_FASTA_SEQUENCE_STATE 3
#define ALIGN_FASTA_WS_STATE       4
#define ALIGN_FASTA_EOF_STATE      5


#define ALIGN_FASTA_INITIAL_BUFSIZE   2
#define ALIGN_FASTA_INITIAL_SEQNUM    2


/* some data structures */
typedef struct _ALIGN_FASTA_STRUCT {
  
  char *buf;
  long int bufsize;
  int type;

} ALIGN_FASTA_TOKEN;



/* some function prototypes */
/*
static inline
int
ALIGN_FASTA_is_alpha(char c);

static inline
int
ALIGN_FASTA_is_whitespace(char c);

static inline
int
ALIGN_FASTA_is_bracket(char c);

static inline
int
ALIGN_FASTA_get_token(FILE *fp,
		      ALIGN_FASTA_TOKEN *token);

static inline
void
ALIGN_remove_leading_ws(char *buf);
*/

int
ALIGN_read_FASTA(ALIGN_seqset *seqset,
		 char *filename);

char *
ALIGN_FASTA_add_seqset(ALIGN_seqset *seqset,
		       int seq,
		       char *buf);

#endif /* _INC_FASTA_H_ */





