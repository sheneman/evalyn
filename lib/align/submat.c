/*
 * submat.c
 *
 * $Id: submat.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
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


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "align.h"
#include "syms.h"
#include "submat.h"





/*
 * ALIGN_print_scorebuf() 
 * 
 * Print a scorebuf
 *
 */
void
ALIGN_print_scorebuf(int *scorebuf,
		     int scorebuf_sz) {
  
  int i;
  

  printf("SCOREBUF:\n");
  
  for(i=0;i<scorebuf_sz;i++) {
    printf("  %d> %d\n", i, scorebuf[i]);
  }
  
  printf("\n");
  
  return;
}





/*
 * ALIGN_add_score() - 
 *
 * If given score doesn't exist in scorebuf, then add it.
 * increment the scorebuf size and return scorebuf size
 * 
 */
int
ALIGN_add_score(int *scorebuf,
		int score) { 

  static int scorebuf_sz = 0;
  int i;
  
  for(i=0;i<scorebuf_sz;i++) {
    if(scorebuf[i] == score) {
      return(scorebuf_sz);
    }
  }
  
  scorebuf[scorebuf_sz++] = score;

  return(scorebuf_sz);
}




/*
 * ALIGN_read_submat_pass1() - 
 *
 * Here we build the unique list of scores in the submat
 * and we fill in the symbols.
 *
 * These steps are a prerequisite to the second pass
 *
 */
void
ALIGN_read_submat_pass1(FILE *fp,
			ALIGN_submat *submat) {

  char buf[128], *string;
  char firstline_flag = 0;
  int i, j, *scorebuf=NULL, scorebuf_sz=0;
  double score;


  /* make sure fp is valid */
  if(!fp) {
    fprintf(stderr, "ALIGN: NULL File Pointer in ALIGN_read_submat_pass1()\n");
    exit(-1);
  }
  
  /* rewind our file */
  rewind(fp);
  
  /* a temporary buffer for tracking unique scores */
  scorebuf = (int *)calloc(submat->symtab->nsyms*submat->symtab->nsyms, sizeof(int));

  /* read the whole file here, ignoring comment lines */
  /* this is a 2-pass thing, the first pass is needed to build the scores vector */
  firstline_flag = 1;
  j = 0;
  while(!feof(fp)) {
    
    i=0;
    fgets(buf, 128, fp);
    
    if(buf[0]!='#') {

      /* if this is the first line in the matrix file */
      if(firstline_flag) {
	string = (char *)strtok(buf, " \t\n");
	if(string && string[0] >= WILD_CHAR) {
	  submat->symtab->syms[i++] = (toupper(string[0]));
	}

	while((string = (char *)strtok(NULL, " \t\n"))) {
	  if(string && string[0] >= WILD_CHAR) {
	    submat->symtab->syms[i++]=(toupper(string[0]));
	  }
	}

	firstline_flag = 0;

      } else {  /* not the first line in the matrix file */
	
	string = (char *)strtok(buf, " \t\n");

	if(j < submat->symtab->nsyms) {
	  if(submat->symtab->syms[j] != toupper(string[0])) {
	    fprintf(stderr, "ALIGN: Format error in substitution matrix.  Symbols along both axes must match.\n");
	  }
	  
	  string = (char *)strtok(NULL, " \t");
	  if(string) {
	    sscanf(string, "%lf", &score);

	    scorebuf_sz = ALIGN_add_score(scorebuf, score);

	    i++;

	    while((string = (char *)strtok(NULL, " \t"))) {
	      if(string && i < (submat->symtab->nsyms) ) {
		sscanf(string, "%lf", &score);

		scorebuf_sz = ALIGN_add_score(scorebuf, score);

		i++;
	      }
	    }
	  }
	}
	j++;
      }
    }
  }

  /* set the substitution matrix type here */
  if(submat->symtab->nsyms > 20) {
    submat->matrix_type = ALIGN_PROTEIN;
  } else {
    submat->matrix_type = ALIGN_DNA;
  }

  /* ALIGN_print_scorebuf(scorebuf, scorebuf_sz); */
  submat->scores  = ALIGN_make_scores(scorebuf, scorebuf_sz);
  submat->nscores = scorebuf_sz;

  /* free the temporary score buffer here */
  free(scorebuf);

  rewind(fp);

  return;
}







/*
 * ALIGN_read_submat_pass2() - 
 *
 * Here, all we do is determine the score category and fill in the 
 * substitution matrix
 *
 */
void
ALIGN_read_submat_pass2(FILE *fp,
			ALIGN_submat *submat) {

  char buf[128], *string;
  char firstline_flag = 0;
  int i, j;
  double score;


  /* make sure fp is valid */
  if(!fp) {
    fprintf(stderr, "ALIGN: NULL File Pointer in ALIGN_read_submat_pass2()\n");
    exit(-1);
  }
  
  /* rewind our file */
  rewind(fp);
  
  /* read the whole file here, ignoring comment lines */
  /* this is a 2-pass thing, the first pass is needed to build the scores vector */
  firstline_flag = 1;
  j = 0;
  while(!feof(fp)) {
    
    i=0;
    fgets(buf, 128, fp);

    if(buf[0]!='#') {

      /* if this is the first line in the matrix file */
      if(firstline_flag) {
	firstline_flag = 0;

      } else {  /* not the first line in the matrix file */
	
	string = (char *)strtok(buf, " \t");

	if(j < submat->symtab->nsyms) {
	  
	  string = (char *)strtok(NULL, " \t");
	  if(string) {
	    sscanf(string, "%lf", &score);

	    VECTOR(submat, i, j) = score;

	    i++;

	    while((string = (char *)strtok(NULL, " \t"))) {
	      if(string && i < (submat->symtab->nsyms) ) {
		sscanf(string, "%lf", &score);

		VECTOR(submat, i, j) = score; 

		i++;
	      }
	    }
	  }
	}
	j++;
      }
    }
  }

  rewind(fp);

  return;
}





/*
 * ALIGN_determine_submat_syms() - 
 *
 * Determine the number of symbols in the substitution matrix
 *
 */
int
ALIGN_count_submat_syms(FILE *fp) {

  int firstline_flag=1, nsyms=0;
  char buf[128], *string;

  /* make sure fp is valid */
  if(!fp) {
    fprintf(stderr, "ALIGN: NULL File Pointer in ALIGN_count_submat_syms()\n");
    exit(-1);
  }
  
  /* rewind our file */
  rewind(fp);

  /* read the whole file here, ignoring comment lines */
  while(firstline_flag && !feof(fp)) {
    
    fgets(buf, 128, fp);

    if(buf[0] != '#') {

      /* if this is the first line in the matrix file */
      if(firstline_flag) {
	string = (char *)strtok(buf, " \t");
	if(string && string[0] >= WILD_CHAR) {
	  nsyms++;
	}

	while((string = (char *)strtok(NULL, " \t"))) {
	  if(string && string[0] >= WILD_CHAR) {
	    nsyms++;
	  }
	}
	firstline_flag = 0;
      }
    }
  }

  rewind(fp);

  return(nsyms);
}




/*
 * ALIGN_print_submat() - 
 *
 * Pretty-print the substitution matrix *
 */
void
ALIGN_print_submat(ALIGN_submat *submat) {

  int i, j;


  printf("In ALIGN_print_submat()\n");  

  /* bail if we have some NULL pointers */
  if(!submat || !submat->symtab || 
     !submat->symtab->syms || !submat->vector) {

    fprintf(stderr, "ALIGN: Submat is corrupt.\n");
    exit(-1);
  }

  /* print the symbol table first */
  printf("SYMBOL TABLE:\n\n");
  printf("  nsyms = %d\n", submat->symtab->nsyms);
  printf("  ");
  for(i=0;i<submat->symtab->nsyms;i++) {
    printf("%c ", submat->symtab->syms[i]);
  }
  printf("\n\n");
  
  
  printf("SUBSTITUTION MATRIX:\n\n");

  /* print the symbols across the top here */
  for(i=0;i<submat->symtab->nsyms;i++) {
    printf("%6c", submat->symtab->syms[i]);
  }
  printf("\n");
  
  for(i=0;i<submat->symtab->nsyms;i++) {

    for(j=0;j<submat->symtab->nsyms;j++) {
      
      /* print the symbols across the left */
      if(j==0) {
	printf("%c ", submat->symtab->syms[i]);
      }
      printf("%6.2lf", VECTOR(submat, i, j)); 
    }
    printf("\n");
  }

  printf("\n\n");
  
  return;
}






/*
 * ALIGN_read_submat() - 
 *
 * Read a substitution matrix from disk and return it
 *
 */
ALIGN_submat *
ALIGN_read_submat(char *matrix_filename) {

  FILE *fp;
  ALIGN_submat *submat;

  /* some debugging output */
  printf("In ALIGN_read_submat()\n");
  printf("  matrix_filename = %s\n", matrix_filename);

  /* open matrix file here */
  fp = fopen(matrix_filename, "r");
  if(!fp) {
    fprintf(stderr, "ALIGN: Could not open substitution matrix: %s\n", matrix_filename);
    exit(-1);
  }

  /* dynamically allocate the submat, the symbol array, and the vector */
  submat = (ALIGN_submat *)calloc(1, sizeof(ALIGN_submat));

  /* 
   * dynamically determine the number of symbols in the substitution matrix
   * and allocate a symbol table of the appropriate size
   */
  submat->symtab = 
    ALIGN_alloc_symtab(ALIGN_count_submat_syms(fp));
  
  /* this is the substitution matrix, which is represented as a 1D array
     for cache performance reasons
  */
  printf("in ALIGN_read_submat() -- %d\n", submat->symtab->nsyms);
  
  submat->vector = 
    (double *)calloc( ( (submat->symtab->nsyms) * (submat->symtab->nsyms) ), sizeof(double) );

  /* this takes two passes */
  ALIGN_read_submat_pass1(fp, submat); 
  ALIGN_read_submat_pass2(fp, submat);
  
  fclose(fp);
  
  /* automatically determine type of residues (protein or dna) in substitution matrix */
  if( submat->symtab->nsyms > ALIGN_NUM_DNA_SYMS) {
    submat->matrix_type = ALIGN_PROTEIN;
  } else {
    submat->matrix_type = ALIGN_DNA;
  }
  
  return(submat);
}








/*
 * ALIGN_free_submat() -
 *
 * Free a substitution matrix
 */
void
ALIGN_free_submat(ALIGN_submat *submat) {
  
  /* free the dynamically allocated symbols from the matrix */
  ALIGN_free_symtab(submat->symtab);
  
  /* free the dynamically allocated vector */
  if(submat && submat->vector) {
    free(submat->vector);
  }
  
  /* free the scores array here */
  if(submat && submat->scores) {
    free(submat->scores);
  }
  
  /* free the submat itself */
  if(submat) {
    free(submat);
  }
  
  return;
}













/*
 * ALIGN_make_scores()
 *
 * Given a scorebuf, make a proper score vector 
 *
 */
short int *
ALIGN_make_scores(int *scorebuf,
		  int scorebuf_sz) {
  
  int i;
  short int *scores;
  
  scores = (short int *)calloc(scorebuf_sz, sizeof(short int));

  for(i=0;i<scorebuf_sz;i++) {
    scores[i] = scorebuf[i];
  }
  
  return(scores);
}





/*
 * ALIGN_lookup_score()
 *
 * Linearly search for a score and return its index
 *
 */
unsigned int
ALIGN_lookup_score(ALIGN_submat *submat,
		   short int score) {

  int i;
  
  /* just a simple linear search */
  for(i=0;i<submat->nscores;i++) {
    if(submat->scores[i] == score) {
      return(i);
    }
  }

  fprintf(stderr, "ALIGN: Looking up a score which does not exist.  Something very wrong.\n");
  exit(-1);
}




/*
 * ALIGN_print_submat_syms()
 *
 * Print the symbols in the submat
 *
 */
void
ALIGN_print_submat_syms(ALIGN_submat *submat) {

  int i;
  
  printf("Number of Symbols in Submat: %d\n", submat->symtab->nsyms);
  for(i=0;i<submat->symtab->nsyms;i++) {
    printf("%d> [%c]\n", i, submat->symtab->syms[i]);
  }
  
  return;
}



