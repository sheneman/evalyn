/*
 * parsimony.c
 *
 * $Id: parsimony.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 * 
 *
 * Copyright (c) 2003 University of Idaho.  All rights reserved.
 *
 *
 * Some routines for computing Fitch parsimony scores.
 *
 * The implementation is as described in Joseph Felsenstein's
 * "Inferring Phylogenies" (Sinauer 2004)
 *
 * Optimizations made using Gladstein's Parsimony Calculation optimizations,
 * from:
 *
 *   David S. Gladstein
 *   "Efficient Incremental Character Optimization"
 *   The Willi Hennig Society, 1997
 *
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <align.h>

#include "ltree.h"
#include "parsimony.h"



/*************************************
 *                                   *
 * BEGIN STATIC FUNCTION PROTOTYPES  *
 *                                   *
 *************************************/

static int
LTREE_parsimony_informative(ALIGN_alignment *alignment,
			    int column);
static int
LTREE_fitch_intersection(LTREE_treenode *node);

static int
LTREE_fitch_union(LTREE_treenode *node);

static int
LTREE_fitch_find(LTREE_treenode *node, char c);

static 
void LTREE_print_pbuf(char *pbuf,
		      int pbuf_p);
static int
LTREE_parsimony_recurse(LTREE_treenode *node,
			ALIGN_alignment *alignment,
			int column);

/***********************************
 *                                 *
 * END STATIC FUNCTION PROTOTYPES  *
 *                                 *
 ***********************************/




/*
 * LTREE_parsimony_score() -
 *
 * Return the parsimony score of the given tree
 *
 * Note that I treat gaps as wildcard states (any state), such that 
 * the intersection of an A and a GAP is an A with no state change.
 *
 */
int
LTREE_parsimony_score(LTREE_tree *tree) {

  int i, score = 0;
  
  /* across all columns */
  for(i=0;i<tree->root->alignment->n;i++) {
    
    /* for each column in the alignment, compute the parsimony scores */
    switch(LTREE_parsimony_informative(tree->root->alignment, i)) {
    case 1:

      score += 
	LTREE_parsimony_recurse(tree->root,
				tree->root->alignment,
				i);
      break;

    default: 

      break;

    }

  }

  return(score);
}





/*
 * LTREE_parsimony_recurse() - 
 *
 * The recursive component of the parsimony scorer
 *
 */
static int
LTREE_parsimony_recurse(LTREE_treenode *node,
			ALIGN_alignment *alignment,
			int column) {

  
  int score = 0;
  int index;


  /* if this is a leaf node, then setup initial pbuf */
  if(!node->left &&
     !node->right) {
    
    index = ALIGN_find_seq(alignment, node->seq_id);
    if(index < 0) {
      fprintf(stderr, "EVALYN: seq_id invalid in LTREE_parsimony_recurse\n");
      exit(0);
    }

    //    node->pbuf[0] = alignment->alignment[column][index];
    node->pbuf[0] = alignment->text[ALIGN_index(alignment, column, index)];
    node->pbuf_p = 1;
  }

  /* recurse left */
  if(node->left) {

    score += 
      LTREE_parsimony_recurse(node->left,
			      alignment,
			      column);
  }

  /* recurse right */
  if(node->right) {

    score += 
      LTREE_parsimony_recurse(node->right,
			      alignment,
			      column);
  }
  
  /* Fitch Parsimony */
  if(node->left && node->right) {
    if(!LTREE_fitch_intersection(node)) {
      LTREE_fitch_union(node);
      score++;
    }
  }
  
  return score;
}



 


/*
 *
 * LTREE_parsimony_informative() -
 *
 * Determine if the site (alignment column) is informative.
 * Use alignment profile to do this.
 *
 * An informative column must have at least 2 unique symbols
 *
 * RETURNS:
 * --------
 * 
 *   0 = uninformative
 *   1 = informative
 *
 * If 0, the column can be ignored, as it will not add
 * anything to the parsimony score.  If 1, the column
 * is informative in the general sense.  
 *
 */
static int
LTREE_parsimony_informative(ALIGN_alignment *alignment,
			    int column) {
  
  int i;
  int nsyms;
  
  nsyms = alignment->symtab->nsyms;
  
  /* iterate across all symbols in the profile */
  /* this loop will need to change in order to be fast */
  for(i=0;i<nsyms;i++) {

    if(alignment->profile[PROFILE(nsyms, column, i)]) {
      if(alignment->profile[PROFILE(nsyms, column, i)]) {
	return(0);
      } else {
	return(1);
      }
    }
  }
  
  exit(-1);
}







/*
 * LTREE_fitch_intersection() - 
 *
 * Compute the intersection of two nodes for use in Fitch parsimony
 *
 * RETURNS:
 * --------
 *
 *  Number of intersections
 *
 */
static int
LTREE_fitch_intersection(LTREE_treenode *node) {
  
  int i, j, k = 0;
  int intersection_flag = 0;

  /*
  printf("In LTREE_fitch_intersection() - \n");
  printf("LEFT PBUF:\n\t");
  LTREE_print_pbuf(tree->left->pbuf, tree->left->pbuf_p);
  printf("RIGHT PBUF:\n\t");
  LTREE_print_pbuf(tree->right->pbuf, tree->right->pbuf_p);
  */
  
  /* for every symbol in left pbuf */
  for(i=0;i<node->left->pbuf_p;i++) {

    /* right pbuf */
    for(j=0;j<node->right->pbuf_p;j++) {

      if(ISGAP(node->alignment->symtab,  node->left->pbuf[i]) &&
	 !ISGAP(node->alignment->symtab, node->right->pbuf[j])) {

	node->pbuf[k++] = node->right->pbuf[j];
	intersection_flag = 1;
      }

      if(!ISGAP(node->alignment->symtab, node->left->pbuf[i]) &&
	 ISGAP(node->alignment->symtab, node->right->pbuf[j])) {

	node->pbuf[k++] = node->left->pbuf[i];
	intersection_flag = 1;
      }

      if(ISGAP(node->alignment->symtab, node->left->pbuf[i]) &&
	 ISGAP(node->alignment->symtab, node->right->pbuf[j])) {

	node->pbuf[k++] = node->alignment->symtab->nsyms-1;
	intersection_flag = 1;

      } else {

	if(node->left->pbuf[i] == node->right->pbuf[j]) {
	  node->pbuf[k++] = node->left->pbuf[i];

	  intersection_flag = 1;
	}
      }
    }
  }
  
  node->pbuf_p = k;

  return(intersection_flag);
}






/*
 *
 * LTREE_fitch_union() - 
 *
 * Compute the union of two nodes for use in Fitch parsimony computation
 *
 * RETURNS:
 * --------
 *
 *   Number of symbols in the union
 * 
 */
static int
LTREE_fitch_union(LTREE_treenode *node) {
  
  int i;

  
  node->pbuf_p = 0;

  /* iterate through left */
  for(i=0;i<node->left->pbuf_p;i++) {
    if(!LTREE_fitch_find(node, node->left->pbuf[i])) {
      node->pbuf[(int)(node->pbuf_p++)] = node->left->pbuf[i];
    }
  }

  /* iterate through right */
  for(i=0;i<node->right->pbuf_p;i++) {
    if(!LTREE_fitch_find(node, node->right->pbuf[i])) {
      node->pbuf[(int)(node->pbuf_p++)] = node->right->pbuf[i];
    }
  }
  
  return(node->pbuf_p);
}






/*
 * LTREE_fitch_find() - 
 *
 * Find a given symbol in a node for use in Fitch parsimony
 *
 * RETURNS:
 * --------
 *
 *   0 if symbol absent
 *   1 if symbol present
 *
 */
static int
LTREE_fitch_find(LTREE_treenode *node,
		 char c) {

  int i;
		
  /* just do a linear search here */
  for(i=0;i<node->pbuf_p;i++) {
    if(node->pbuf[i] == c) {
      return(1);
    }
  }
  
  return(0);
}





/*
 * LTREE_test_parsimony() - 
 *
 * A function which manually builds a simple tree and tests parsimony
 *
 * Here, we use the alignment
 *
 *   1> CATTTGAT-A
 *   2> AA-TTCAT-A
 *   3> CAATTTCA-T
 *   4> AA-ACAC--T
 *   5> GA-TCACTAT
 *
 * Which is constructed on a tree which looks like:
 *
 *  ((1,2),(3,(4,5)))
 *
 * That tree topology is constructed manually.
 *
 * The proper parsimony score for that alignment + tree is 12
 *
 */
void
LTREE_test_parsimony(void) {
  
  int i, j, pscore;
  LTREE_treenode *nodes[9];
  LTREE_tree *tree;
  char buf[2], c;
  char cs[5][10] = { "CATTTGAT-A",
		     "AA-TTCAT-A",
		     "CAATTTCA-T",
		     "AA-ACAC--T",
		     "GA-TCACTAT"  };


  ALIGN_scoring_system ss;

  


  /* load out substitution matrix here */
  ss.submat = ALIGN_read_submat(DEFAULT_DNA_MATRIX);
  if(!ss.submat) {
    fprintf(stderr, "EVALYN: Fatal error parsing substitution matrix file\n");
    exit(-1);
  }

  /* convert the character sequences above to this */
  for(i=0;i<5;i++) {
    for(j=0;j<10;j++) {
      c = cs[i][j];
      cs[i][j] = ALIGN_symtab_lup(ss.submat->symtab, cs[i][j]);
      if(cs[i][j] < 0) {
	printf("looked up [%c] and got -1\n", c);
      }
    }
  }
  
  /* print the converted foo here */
  printf("The converted matrix\n");
  for(i=0;i<5;i++) {
    for(j=0;j<10;j++) {
      printf("%2d ", cs[i][j]);
    }
    printf("\n");
  }
  printf("\n\n");
  
  ss.gap_open   = -999.0;
  ss.gap_extend = -999.0;


  printf("EVALYN: In LTREE_test_parsimony() - \n");

  /* allocate all of the nodes */
  for(i=0;i<9;i++) {
    nodes[i] = (LTREE_treenode *)calloc(1, sizeof(LTREE_treenode));
    nodes[i]->seq_id = -1;
  }
  
  /* construct the leaves */
  buf[1] = '\0';
  for(i=0;i<5;i++) {

    nodes[i]->seq_id = i;
    nodes[i]->remove_flag = 0;
    nodes[i]->left   = NULL;
    nodes[i]->right  = NULL;
    nodes[i]->parent = NULL;

    nodes[i]->alignment = ALIGN_init_alignment(cs[i], ss.submat->symtab, 10, i);
  }
  
  /* Node that connects C and A */
  nodes[5]->left  = nodes[0];
  nodes[5]->right = nodes[1];
  nodes[0]->parent = nodes[1]->parent = nodes[5];

  /* Node that connects A and G */
  nodes[6]->left = nodes[3];
  nodes[6]->right = nodes[4];
  nodes[3]->parent = nodes[4]->parent = nodes[6];
  
  /* Node that connects C to (AG) */
  nodes[7]->left= nodes[2];
  nodes[7]->right = nodes[6];
  nodes[2]->parent = nodes[6]->parent = nodes[7];

  /* Root node, joining (CA) and (C,(AG)) */
  nodes[8]->left = nodes[5];
  nodes[8]->right = nodes[7];
  nodes[5]->parent = nodes[7]->parent = nodes[8];

  tree = (LTREE_tree *)calloc(1, sizeof(LTREE_tree));
  tree->root = nodes[8];

  /* construct the alignment */  
  LTREE_align_tree(tree, &ss);
  
  ALIGN_print_alignment(tree->root->alignment, NULL);
  
  pscore = LTREE_parsimony_score(tree);
  printf("Parsimony Score = %d\n", pscore);
  
  return;
}




/*
 * LTREE_print_pbuf() - 
 *
 * A debugging function used for printing the contents of the 
 * pbuf, which is used in parsimony scoring.  
 *
 */
static
void LTREE_print_pbuf(char *pbuf, 
		      int pbuf_p) {
  int i;
  
  printf("PBUF: ");
  for(i=0;i<pbuf_p;i++) {
    printf("[%d]", pbuf[i]);
  }
  printf("\n");
}



