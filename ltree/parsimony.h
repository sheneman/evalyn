/*
 * parsimony.h
 *
 * $Id: parsimony.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2003 University of Idaho.  All rights reserved.
 *
 *
 * Routines for performing parsimony scoring
 * 
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */


#ifndef _INC_PARSIMONY_H_
#define _INC_PARSIMONY_H_

#include "ltree.h"



/* some function prototypes */
int 
LTREE_parsimony_score(LTREE_tree *tree);

void
LTREE_test_parsimony(void);


#endif /* _INC_PARSIMONY_H_ */












