/*
 * syms.c 
 *
 * $Id: syms.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2003 Luke Sheneman.  All rights reserved.
 *
 * A collection of lookup routines that look for specific symbols
 * in one of the classified symbol tables
 *
 * These functions look highly inefficient, and they are in
 * worst-case scenarios.  However, in practice, this will not 
 * be a problem.
 *
 *
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */


#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "syms.h"





/*******************************
 *                             *
 *  BEGIN STATIC DECLARATIONS  *
 *                             *
 *******************************/


/* Four DNA symbols */
static char ALIGN_dna_syms[ALIGN_NUM_DNA_SYMS] = { 'A', 'G', 'C', 'T' };

/* Four DNA symbols */
static char ALIGN_rna_syms[ALIGN_NUM_RNA_SYMS] = { 'A', 'G', 'C', 'U' };

/* All 12 IUPAC-IUB/GCG DNA ambiguity codes */
static char ALIGN_dna_ambiguity_syms[ALIGN_NUM_DNA_AMBIGUITY_SYMS] =
{
  'M', 'R', 'W', 'S', 'Y', 'K',
  'V', 'H', 'D', 'B', 'X', 'N',
  'A', 'C', 'G', 'T'
};
 
/* All 23 protein residue symbols */
static char ALIGN_protein_syms[ALIGN_NUM_PROTEIN_SYMS] = 
{
  'A', 'R', 'N', 'B', 'D', 
  'C', 'Q', 'Z', 'E', 'G', 
  'H', 'I', 'L', 'K', 'M', 
  'F', 'P', 'S', 'T', 'W', 
  'Y', 'V', 'X'  
};


/* 
 *  DNA substitution matrices might include ambiguity codes.
 *  We track the 8 symbols which indicate amino acids but are
 *  not valid DNA nucleotides.  This is useful when automatically
 *  classifying a set of sequences or an inputted substitution 
 *  matrix as protein or DNA
 */   

#ifdef _NOT_USED_
static char ALIGN_nondna_syms[ALIGN_NUM_NONDNA_SYMS] = 
{ 'Q', 'Z', 'E', 'I', 'L', 'F', 'P' };
#endif /* _NOT_USED_ */


/* 
 * The following is all valid symbols from the gap to the 
 * proteins (with ambiguity) and DNA/RNA (with ambiguity).
 * In short, if the symbol isn't here, its not DNA, RNA or
 * protein and we can't handle it
 * 
 * (may need to add support for '*' later) 
 *
 */
static char ALIGN_all_syms[ALIGN_NUM_ALL_SYMS] = 
{ '-', 'A', 'R', 'N', 'B', 'D', 
  'C', 'Q', 'Z', 'E', 'G', 'H', 
  'I', 'L', 'K', 'M', 'F', 'P', 
  'S', 'T', 'W', 'Y', 'V', 'X',
  'U'
 };


/*******************************
 *                             *
 *   END STATIC DECLARATIONS   *
 *                             *
 *******************************/







/*
 * ALIGN_SYM_is_DNA() - 
 *
 * Answers the question "is this symbol a valid DNA symbol?"
 *
 * 1 = YES
 * 0 = NO
 *
 */
int
ALIGN_SYM_is_DNA(char c) {
  
  int i;
  
  for(i=0;i<ALIGN_NUM_DNA_SYMS;i++) {
    if(ALIGN_dna_syms[i] == toupper(c)) {
      return(1);
    }
  }

  return(0);
}





/*
 * ALIGN_SYM_is_DNA_AMBIGUITY() - 
 *
 * Answers the question "is this symbol a valid DNA symbol (with ambiguity codes)?"
 *
 * 1 = YES
 * 0 = NO
 *
 */
int 
ALIGN_SYM_is_DNA_AMBIGUITY(char c) {

  int i;
  
  for(i=0;i<ALIGN_NUM_DNA_AMBIGUITY_SYMS;i++) {
    if(ALIGN_dna_ambiguity_syms[i] == toupper(c)) {
      return(1);
    }
  }

  return(0);
}




/*
 * ALIGN_SYM_is_RNA() - 
 *
 * Answers the question "is this symbol a valid RNA symbol?"
 *
 * 1 = YES
 * 0 = NO
 *
 */
int 
ALIGN_SYM_is_RNA(char c) {

  int i;
  
  for(i=0;i<ALIGN_NUM_RNA_SYMS;i++) {
    if(ALIGN_rna_syms[i] == toupper(c)) {
      return(1);
    }
  }

  return(0);
}



/*
 * ALIGN_SYM_is_PROTIEN() - 
 *
 * Answers the question "is this symbol a valid RNA symbol?"
 *
 * 1 = YES
 * 0 = NO
 *
 */
int 
ALIGN_SYM_is_PROTEIN(char c) {

  int i;
  
  for(i=0;i<ALIGN_NUM_PROTEIN_SYMS;i++) {
    if(ALIGN_protein_syms[i] == toupper(c)) {
      return(1);
    }
  }

  return(0);
}





/*
 * ALIGN_SEQ_is_DNA() - 
 *
 * Answers the question "Is this a DNA sequence?"
 *
 * 1 = YES
 * 0 = NO
 *
 */
int
ALIGN_SEQ_is_DNA(char *sequence) {
  
  int i;
  
  for(i=0;i<strlen(sequence);i++) {
    if(!ALIGN_SYM_is_DNA(sequence[i])) {
      return(0);
    }
  }

  return(1);
}




/*
 * ALIGN_SEQ_is_DNA_AMBIGUITY() - 
 *
 * Answers the question "Is this a DNA sequence (with ambiguity codes)?"
 *
 * 1 = YES
 * 0 = NO
 *
 */
int
ALIGN_SEQ_is_DNA_AMBIGUITY(char *sequence) {
  
  int i;
  
  for(i=0;i<strlen(sequence);i++) {
    if(!ALIGN_SYM_is_DNA_AMBIGUITY(sequence[i])) {
      return(0);
    }
  }

  return(1);
}



/*
 * ALIGN_SEQ_is_RNA() - 
 *
 * Answers the question "Is this an RNA sequence?"
 *
 * 1 = YES
 * 0 = NO
 *
 */
int
ALIGN_SEQ_is_RNA(char *sequence) {
  
  int i;
  
  for(i=0;i<strlen(sequence);i++) {
    if(!ALIGN_SYM_is_RNA(sequence[i])) {
      return(0);
    }
  }

  return(1);
}





/*
 * ALIGN_SEQ_is_PROTEIN() - 
 *
 * Answers the question "Is this a protein sequence?"
 *
 * 1 = YES
 * 0 = NO
 *
 */
int
ALIGN_SEQ_is_PROTEIN(char *sequence) {
  
  int i;
  
  for(i=0;i<strlen(sequence);i++) {
    if(!ALIGN_SYM_is_PROTEIN(sequence[i])) {
      return(0);
    }
  }

  return(1);
}




/*
 * ALIGN_lookup_index()
 *
 * Given a symbol, lookup its index
 * This can be slow (linear), because it should only be used at startup
 * when loading the substitution matrices
 *
 * RETURNS:
 * --------
 *
 *  -1 if specified character is not in table
 *  index into table, otherwise
 * 
 */
int
ALIGN_lookup_index(char c) {
  
  int i;
  
  for(i=0;i<ALIGN_NUM_ALL_SYMS;i++) {
    if(c == ALIGN_all_syms[i]) {
      return(i);
    }
  }
  
  return(-1);
}














