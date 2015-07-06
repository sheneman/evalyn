/*
 * fasta.c
 *
 * $Id: fasta.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "align.h"
#include "fasta.h"




/*
 *
 * ALIGN_FASTA_is_whitespace() - determine if character is a whitespace character
 *
 * INPUT:
 * ------
 *  c -- character to test
 *
 * RETURN:
 * -------
 *   int -- 1 if character is whitespace (space, tab, CR, LF)
 *          0 if character is NOT whitespace
 *
 */
static inline
int
ALIGN_FASTA_is_whitespace(char c) {

  if( c == ' '  ||   /* space           */
      c == '\n' ||   /* newline         */
      c == '\r' ||   /* carriage-return */
      c == '\v' ||   /* vertical tab    */
      c == '\f' ||   /* form feed       */
      c == '\t' ) {  /* horizontal tab  */
    return(1);
  } else {
    return(0);
  }

}




/*
 *
 * ALIGN_FASTA_is_bracket() - determine if character is an open bracket
 *
 * INPUT:
 * ------
 *  c -- character to test
 *
 * RETURN:
 * -------
 *   int -- 1 if character is a bracket
 *          0 if character is NOT a bracket
 *
 */
static inline
int
ALIGN_FASTA_is_bracket(char c) {
  if(c == '>') {
    return(1);
  } else {
    return(0);
  }
}




/*
 * ALIGN_remove_leading_ws() - 
 * 
 * A simple function for removing the leading whitespace from
 * sequence titles 
 *
 */
static inline
void
ALIGN_remove_leading_ws(char *buf) {

  int length, i;
  int leading_flag;
  char *tmp;
  
  length = strlen(buf);
  
  tmp = (char *)calloc(length+1, sizeof(char));
  
  leading_flag = 1;
  for(i=0;i<length;i++) {

    if(leading_flag) {
      if(buf[i] != ' ' && buf[i] != '\t' ) {
	leading_flag = 0;
	strcpy(tmp, &buf[i]);
      }
    }
  }
  
  strcpy(buf, tmp);
  free(tmp);
  
  return;
}






/*
 * ALIGN_read_FASTA() - 
 *
 * Read a FASTA file using a simple finite state automaton parser
 *
 */
int
ALIGN_read_FASTA(ALIGN_seqset *seqset,
		 char *filename) {

  FILE *fp = NULL;
  char *buf = NULL;
  int bufsize;
  int index, seq;
  int first_sequence_flag;
  int c, state;
  char *uniq;



  /* allocate our initial buffer */
  bufsize = ALIGN_FASTA_INITIAL_BUFSIZE;
  buf = (char *)calloc(bufsize, sizeof(char));
  
  /* allocate out seqset foo here */
  seqset->num = ALIGN_FASTA_INITIAL_SEQNUM;
  seqset->seq = (ALIGN_seq *)calloc(seqset->num, sizeof(ALIGN_seq));

  /* open the file for reading */
  fp = fopen(filename, "r");
  if(!fp) {
    fprintf(stderr, "ALIGN:  Failed to open input FASTA file: %s\n", filename);
    goto XIT_BAD;
  }

  first_sequence_flag = 1;

  index      = 0;  /* this is the pointer into the buffer */
  seq        = 0;
  
  state = ALIGN_FASTA_UNKNOWN_STATE;  /* the initial state */

  while(1) {
    
    /* get a character */
    c = fgetc(fp);
    if(feof(fp)) {
      if(state == ALIGN_FASTA_SEQUENCE_STATE ||
	 state == ALIGN_FASTA_WS_STATE) {
	buf[index] = '\0';
	
	if(seq >= seqset->num) {
	  seqset->num *= 2;
	  seqset->seq = (ALIGN_seq *)realloc(seqset->seq, seqset->num*sizeof(ALIGN_seq));
	}

	uniq = 
	  ALIGN_FASTA_add_seqset(seqset,
				 seq++,
				 buf);
	if(uniq == NULL) {
	  fprintf(stderr, "ALIGN: Error parsing fasta file.\n");
	  goto XIT_BAD;
	}

      } else {
	fprintf(stderr, "ALIGN: Premature EOF in fasta sequence file.\n");
	goto XIT_BAD;
      }
      break;
    }

    /* expand our buf by 2x if necessary */
    if(index >= bufsize) {
      bufsize *= 2;
      buf = (char *)realloc(buf, bufsize);
    }
    
    switch(state) {
    case ALIGN_FASTA_UNKNOWN_STATE:
      if(ALIGN_FASTA_is_whitespace(c)) {
	/* do nothing */
      } else if(ALIGN_FASTA_is_bracket(c)) {
	state = ALIGN_FASTA_TITLE_STATE;
      } else {
	fprintf(stderr, "ALIGN: Unexpected characters in FASTA sequence file.\n");
	goto XIT_BAD;
      }

      break;

    case ALIGN_FASTA_TITLE_STATE:

      if(c == '\n') {
	state = ALIGN_FASTA_SEQUENCE_STATE;
	buf[index] = '\0';

	if(seq >= seqset->num) {
	  seqset->num *= 2;
	  seqset->seq = (ALIGN_seq *)realloc(seqset->seq, (seqset->num*sizeof(ALIGN_seq)) );
	}

	ALIGN_remove_leading_ws(buf);
	strcpy(seqset->seq[seq].title, buf);
	index = 0;
      } else {
	if(index < MAX_SEQ_TITLE_LENGTH-1) {
	  buf[index++] = c;
	} else {
	  buf[MAX_SEQ_TITLE_LENGTH-1] = '\0';
	}
      }
      break;
      
    case ALIGN_FASTA_WS_STATE:

      if(ALIGN_FASTA_is_bracket(c)) {

	state = ALIGN_FASTA_TITLE_STATE;
	buf[index] = '\0';

	uniq = 
	  ALIGN_FASTA_add_seqset(seqset,
				 seq++,
				 buf);
	if(uniq == NULL) {
	  fprintf(stderr, "ALIGN: Error parsing fasta file.\n");
	  goto XIT_BAD;
	}
	index = 0;

      } else {

	if(!ALIGN_FASTA_is_whitespace(c)) {
	  buf[index++] = c; 
	}
      }

      break;

    case ALIGN_FASTA_SEQUENCE_STATE:

      if(ALIGN_FASTA_is_whitespace(c)) {
	state = ALIGN_FASTA_WS_STATE;
      } else {
	buf[index++] = c;
      }
      break;
      
    default:
      fprintf(stderr, "ALIGN: Unknown state in ALIGN_read_FASTA()\n");
      goto XIT_BAD;
    }

  }

  if(buf) {
    free(buf);
  }

  /* close the file */
  if(fp) {
    fclose(fp);
  }
  
  seqset->num = seq;
  
  /* determine type of sequences based on types of symbols */
  if(strlen(uniq) > ALIGN_NUM_DNA_SYMS) {
    return(ALIGN_PROTEIN);
  } else {
    return(ALIGN_DNA);
  }


 XIT_BAD:

  /* free allocated memory */
  if(buf) {
    free(buf);
  }

  /* close the file */
  if(fp) {
    fclose(fp);
  }
  
  return(-1);
}











/*
 * ALIGN_FASTA_add_seqset()
 *
 * Add a sequence to the seqset
 * 
 */
char *
ALIGN_FASTA_add_seqset(ALIGN_seqset *seqset,
		       int seq,
		       char *buf) {
  int found;
  int j, k, index;
  static char uniq[MAX_SYMS], tmp[2];

  seqset->seq[seq].data = (char *)calloc(strlen(buf)+1, sizeof(char));
  found = 0;
  for(j=0;j<strlen(buf);j++) {

    /* try to add every character to the string of uniq symbols */
    /* if it is not there, add it                               */
    /* we do this to count the number of unique symbols         */
    found = 0;
    for(k=0;k<strlen(uniq) && !found;k++) {
      if(uniq[k] == toupper(buf[j])) {
	found = 1;
	break;
      }
    }
    
    sprintf(tmp, "%c", toupper(buf[j]));
    if(!found) {
      strcat(uniq, tmp);
    }

    /* lookup the symbol here */
    index = ALIGN_seqset_lookup_sym(seqset, buf[j]);
    if(index == -1) {
      fprintf(stderr, "ALIGN: Seq. file contains symbol \"%c\" that is not in substitution matrix.\n", buf[j]);
      return(NULL);
    } else {
      seqset->seq[seq].data[j] = index;
    }
  }

  seqset->seq[seq].length = j;
  
  return(uniq);
}





