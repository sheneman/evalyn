/*
 * rnjhooks.c
 *
 * $Id: rnjhooks.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2005, Luke Sheneman.  All rights reserved.
 *
 *
 * Routines for integrating 
 * 
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */

#include <stdio.h>
#include <string.h>

#include <align.h>
#include <nj_fasta.h>
#include <nj_dist.h>

#include "ltree.h"
#include "rnjhooks.h"




/*
 * LTREE_create_rnj_tree()
 *
 * Construct a Relaxed Neighbor Joining tree
 * and stick it into the population at some 
 * randomly chosen point 
 *
 */
LTREE_tree *
LTREE_create_rnj_tree(ALIGN_scoring_system *ss,
		      LTREE_population *pop,
		      ALIGN_seqset *seqset) {
  

  NJ_DMAT *dmat, *pw_dmat;
  NJ_ARGS rnj_args;
  NJ_TREE *rnj_tree;

  ALIGN_alignment *a, *b, *ltree_alignment;
  NJ_alignment *rnj_alignment;
  LTREE_tree *ltree;
  int i, j;


  printf("In LTREE_create_rnj_tree() - \n");
  

  /*************************************************************
   *
   * ALLOCATE AND INITIALIZE DISTANCE MATRIX
   *
   *************************************************************/
  dmat        = (NJ_DMAT *)calloc(1, sizeof(NJ_DMAT));
  dmat->ntaxa = seqset->num;
  dmat->size  = dmat->ntaxa;

  dmat->taxaname = (char **)calloc(dmat->ntaxa, sizeof(char *));
  
  /* copy taxanames from seqset into dmat */
  for(i=0;i<dmat->ntaxa;i++) {
    dmat->taxaname[i] = (char *)calloc(strlen(seqset->seq[i].title)+1, sizeof(char));
    strncpy(dmat->taxaname[i], seqset->seq[i].title, strlen(seqset->seq[i].title)+1);
  }

  /* allocate the floats in the matrix here */
  dmat->val = (float *)calloc(NJ_NCELLS(dmat->ntaxa), sizeof(float));
  dmat->valhandle = dmat->val;
  
  dmat->rhandle  = dmat->r  = (float *)calloc(dmat->ntaxa, sizeof(float));
  dmat->r2handle = dmat->r2 = (float *)calloc(dmat->ntaxa, sizeof(float));
  
  dmat->maxulps = 4*(dmat->ntaxa+4) + 100;

  /*************************************************************
   *
   *  END OF ALLOCATE AND INITIALIZE OF DISTANCE MATRIX 
   *
   *************************************************************/
  

  /* initialize rnj_args, basically just to choose correction model */
  rnj_args.kimura_flag      = 1;
  rnj_args.jukes_flag       = 0;
  rnj_args.correction_model = NJ_MODEL_KIMURA;
  
  /* lame hack, but this will have to do for now */
  if(seqset->symtab->nsyms < 20) {
    rnj_args.dna_flag     = 1;
    rnj_args.protein_flag = 0;
  } else {
    rnj_args.dna_flag     = 0;
    rnj_args.protein_flag = 1;
  }


  /* 
   * Construct a distance matrix for the input sequences
   * by pairwise aligning each pair of sequences in the 
   * seqset in the context of the specified scoring system.
   * 
   * For now, use k2p distance correction.
   *
   */
  a = b = ltree_alignment = NULL;
  for(i=0;i<seqset->num;i++) {
    
    /* construct an alignment for a, which is sequence i in the seqset */
    if(a) {
      ALIGN_free_alignment(a);
    }

    a = 
      ALIGN_init_alignment(seqset->seq[i].data,
			   seqset->symtab,
			   seqset->seq[i].length,
			   i);

    for(j=i+1;j<seqset->num;j++) {

      /* construct an alignment for b, which is sequence j in the seqset */
      if(b) { 
	ALIGN_free_alignment(b); 
      }

      b = 
	ALIGN_init_alignment(seqset->seq[j].data,
			     seqset->symtab,
			     seqset->seq[j].length,
			     j);

      /* we have both one-sequence alignments, so now do pairwise alignment */
      if(ltree_alignment) { free(ltree_alignment); }
      ltree_alignment = ALIGN_alignment_alignment(ss, a, b);
      
      /* convert alignment formats */
      rnj_alignment = LTREE_ev_alignment_to_rnj(seqset, ltree_alignment);
      
      /* compute the pairwise distances given a pairwise alignment */
      pw_dmat = NJ_compute_dmat(&rnj_args, rnj_alignment);
      
      NJ_free_alignment(rnj_alignment);
      
      /* fill cells in "big" distance matrix (dmat) */
      dmat->val[NJ_MAP(i, j, dmat->ntaxa)] = pw_dmat->val[NJ_MAP(0, 1, 2)];
      
      NJ_free_dmat(pw_dmat);
    }
  }
  
  /* build a RNJ tree based on the distance matrix */
  rnj_tree = NJ_relaxed_nj(&rnj_args, dmat);
  //  rnj_tree = NJ_neighbor_joining(&rnj_args, dmat);
  
  NJ_output_tree(&rnj_args,
		 rnj_tree,
		 dmat,
		 1);
  
  /* convert the RNJ tree to an LTREE */
  ltree = 
    LTREE_rnjtree_to_ltree(rnj_tree,
			   dmat,
			   seqset);
  
  LTREE_print_tree(ltree);
  printf(" ;\n\n");
  
  /* free everything */
  NJ_free_dmat(dmat);                /* free the distance matrix */
  ALIGN_free_alignment(a);           /* free evalyn alignment a  */
  ALIGN_free_alignment(b);           /* free evalyn alignment b  */

  return(ltree);
}




/*
 * LTREE_rnjtree_to_ltree() - 
 *
 * Converts a tree data structure used by Clearcut (RNJ Implementation) 
 * into a usable tree datastructure usable by Evalyn.  
 *
 */
LTREE_tree *
LTREE_rnjtree_to_ltree(NJ_TREE *rnj_tree,
		       NJ_DMAT *dmat,
		       ALIGN_seqset *seqset) {

  LTREE_tree *tree;
  
  tree = (LTREE_tree *)calloc(1, sizeof(LTREE_tree));
  
  tree->root = 
    LTREE_recurse_rnjtree_to_ltree(rnj_tree,
				   dmat,
				   seqset);
  
  tree->seqset = seqset;
  
  return(tree);
}






/*
 * LTREE_recurse_rnjtree_to_ltree() - 
 *
 * Converts a tree data structure used by Clearcut (RNJ Implementation) 
 * into a usable tree datastructure usable by Evalyn.  
 *
 */
LTREE_treenode *
LTREE_recurse_rnjtree_to_ltree(NJ_TREE *rnj_tree,
			       NJ_DMAT *dmat,
			       ALIGN_seqset *seqset) {
  
  LTREE_treenode *ev_root;
  
  if(!rnj_tree) {
    printf("RETURNING NULL IN LTREE_rnjtree_to_ltree()\n");
    return(NULL);
  }

  if(rnj_tree) {
    
    /* allocate the ev_root */
    ev_root = 
      (LTREE_treenode *)calloc(1, sizeof(LTREE_treenode));

    if(rnj_tree->taxa_index == NJ_INTERNAL_NODE) {
      ev_root->seq_id = -1;
    } else {
      ev_root->seq_id = rnj_tree->taxa_index;
    }

  }

  /* recurse left */
  if(rnj_tree->left) {

    ev_root->left =  
      LTREE_recurse_rnjtree_to_ltree(rnj_tree->left,
				     dmat,
				     seqset);

  } else {
    ev_root->left = NULL;
  }
  
  /* recurse right */
  if(rnj_tree->right) {

    ev_root->right = 
      LTREE_recurse_rnjtree_to_ltree(rnj_tree->right,
				     dmat,
				     seqset);
  } else {
    ev_root->right = NULL;
  }

  return(ev_root);

}






/*
 * LTREE_ev_alignment_to_rnj() - 
 *
 * Convert an Evalyn alignment to a clearcut alignment
 *
 */
NJ_alignment *
LTREE_ev_alignment_to_rnj(ALIGN_seqset *seqset,
			  ALIGN_alignment *ev_alignment) {
  
  NJ_alignment *rnj_alignment;
  ALIGN_symtab *symtab;
  int i, j;
  char c;
  
  symtab = seqset->symtab;

  /* allocate our clearcut alignment structure */
  rnj_alignment         = (NJ_alignment *)calloc(1, sizeof(NJ_alignment));
  rnj_alignment->titles = (char **)calloc(ev_alignment->k, sizeof(char *));
  rnj_alignment->data   = (char *)calloc(ev_alignment->k * (ev_alignment->n+1), sizeof(char));

  /* set our dimensions */
  rnj_alignment->nseq   = ev_alignment->k;
  rnj_alignment->length = ev_alignment->n;
  
  for(i=0;i<ev_alignment->k;i++) {

    /* allocate and copy titles */
    rnj_alignment->titles[i] = (char *)calloc(strlen(seqset->seq[i].title), sizeof(char));
    strncpy(rnj_alignment->titles[i], seqset->seq[i].title, strlen(seqset->seq[i].title));

    /* copy alignment text */
    for(j=0;j<ev_alignment->n;j++) {

      c = 
	symtab->syms[(int)(ev_alignment->text[ALIGN_index(ev_alignment, j, i)])];

      if( ( (c >= 'A') && (c <= 'Z') ) ||
	  ( (c >= 'a') && (c <= 'z') ) ) {
	rnj_alignment->data[i*rnj_alignment->length+j] = c;
      } else {
	rnj_alignment->data[i*rnj_alignment->length+j] = NJ_AMBIGUITY_CHAR;
      }
	
    }
  }
  
  return(rnj_alignment);
}
