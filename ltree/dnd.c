/*
 * dnd.c
 *
 * $Id: dnd.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2003 University of Idaho.  All rights reserved.
 *
 *
 * Routines for printing phylogenetic trees (guide trees) in Newick format
 * 
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */

#include <stdio.h>
#include <string.h>

#include <align.h>

#include "ltree.h"
#include "dnd.h"


/*
 * LTREE_write_tree_DND() - 
 *
 * Write the guide tree to a file in the newick (dnd) format
 *
 */
void
LTREE_write_tree_DND(LTREE_tree *tree,
		     char *filename) {
  
  FILE *fp;
  
  fp = fopen(filename, "w");
  if(!fp) {
    fprintf(stderr, "Error:  Could not open DND file %s for writing.\n", filename);
    return;
  }
  
  LTREE_write_tree(tree, fp);

  fprintf(fp, " ;\n");

  fclose(fp);
  
  return;
}




