/*
 * syms.h
 *
 * $Id: syms.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 * Copyright (c) 2003 Luke Sheneman.  All rights reserved.
 *
 * Some tables listing known symbols for biological sequence characters
 *
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */


#ifndef _INC_SYMS_H_ 
#define _INC_SYMS_H_ 1


/* some defines */
#define ALIGN_NUM_DNA_SYMS 4
#define ALIGN_NUM_RNA_SYMS 4
#define ALIGN_NUM_DNA_AMBIGUITY_SYMS 16
#define ALIGN_NUM_PROTEIN_SYMS 23
#define ALIGN_NUM_NONDNA_SYMS 8
#define ALIGN_NUM_ALL_SYMS 25



/* some function prototypes */
int
ALIGN_SYM_is_DNA(char c);

int
ALIGN_SYM_is_DNA_AMBIGUITY(char c);

int
ALIGN_SYM_is_RNA(char c);

int
ALIGN_SYM_is_PROTEIN(char c);

int
ALIGN_SEQ_is_DNA(char *sequence);

int
ALIGN_SEQ_is_DNA_AMBIGUITY(char *sequence);

int
ALIGN_SEQ_is_RNA(char *sequence);

int
ALIGN_SEQ_is_PROTEIN(char *sequence);

int
ALIGN_lookup_index(char c);
 
#endif  /* _INC_SYMS_H_ */







