/*
 * symtab.h
 *
 * $Id: symtab.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2003 Luke Sheneman.  All rights reserved.
 *
 * A set of routines for handling symbol tables, a very important
 * data structure here
 *
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */

#ifndef _INC_SYMTAB_H_
#define _INC_SYMTAB_H_ 1

#include "syms.h"


/* 
   for speed, it is better to statically specify a maximum number
   of symbols here
*/
#define MAX_SYMS ALIGN_NUM_ALL_SYMS


/* some type definitions */
typedef struct ALIGN_STRUCT_symtab {

  char *syms;
  unsigned short int nsyms;
  
} ALIGN_symtab;


/* some function prototypes */
ALIGN_symtab *
ALIGN_alloc_symtab(unsigned short int nsyms);

void
ALIGN_free_symtab(ALIGN_symtab *symtab);

int
ALIGN_symtab_lup(ALIGN_symtab *symtab,
		 char c);



#endif  /* _INC_SYMTAB_H_ */





