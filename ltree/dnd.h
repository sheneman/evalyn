/*
 * dnd.h
 *
 * $Id: dnd.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2003 University of Idaho.  All rights reserved.
 *
 *
 * Routines for writing phylogenetic (guide) trees in Newick format
 *
 * The specification for the Newick tree format can be found in
 * "The Newick tree format", J. Felsenstein (1999), or at 
 * 
 *   http://evolution.genetics.washington.edu/phylip/newicktree.html
 * 
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */

#ifndef _INC_DND_H_
#define _INC_DND_H_ 1

#include "ltree.h"

/* some function prototypes */
void 
LTREE_write_tree_DND(LTREE_tree *ltree,
		     char *filename);

#endif /* _INC_DND_H_ */



