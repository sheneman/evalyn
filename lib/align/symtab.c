/*
 * symtab.c
 *
 * $Id: symtab.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
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


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "syms.h"
#include "symtab.h"
#include "align.h"



/*
 * ALIGN_alloc_symtab()
 *
 * Allocate the symbol table here
 *
 */
ALIGN_symtab *
ALIGN_alloc_symtab(unsigned short int nsyms) {

  ALIGN_symtab *symtab;

  symtab = (ALIGN_symtab *)calloc(1, sizeof(ALIGN_symtab));
  symtab->nsyms = nsyms;
  symtab->syms  = (char *)calloc(symtab->nsyms, sizeof(char));

  return(symtab);
}



/*
 * ALIGN_free_symtab()
 *
 * Free the symbol table here
 *
 */
void
ALIGN_free_symtab(ALIGN_symtab *symtab) {

  if(symtab) {
    if(symtab->syms) {
      free(symtab->syms);
    }
  }

  free(symtab);
  
  return;
}





/*
 * ALIGN_symtab_lup()
 *
 * Lookup the index of a particular symbol
 *
 * RETURNS:
 * --------
 *
 *  -1 if the symbol is not present
 *  otherwise, the index to the symbol
 *
 */
int
ALIGN_symtab_lup(ALIGN_symtab *symtab,
		 char c) {
  int i;
  
  for(i=0;i<symtab->nsyms;i++) {
    if(symtab->syms[i] == toupper(c)) {
      return(i);
    }
  }
  
  return(-1);
}




