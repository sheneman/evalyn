/*
 * cache.h
 *
 * $Id: cache.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2005 Luke Sheneman.  All rights reserved.
 *
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 * 04/01/05
 *
 */



#ifndef _INC_CACHE_H_
#define _INC_CACHE_H_ 1


#include "align.h"


/* some function prototypes */
ALIGN_alignment_cache *
ALIGN_init_alignment_cache(void);

void
ALIGN_free_alignment_cache(void);


void
ALIGN_merge_alignment_ids(ALIGN_alignment *dest,
			  ALIGN_alignment *source1,
			  ALIGN_alignment *source2);

int
ALIGN_determine_id_order(short int *a,
			 short int *b);


short int *
ALIGN_compute_alignment_id(ALIGN_alignment *a,
			   ALIGN_alignment *b);

void
ALIGN_print_alignment_id(short int *id);


ALIGN_cache_hashnode *
ALIGN_find_hashnode(ALIGN_alignment_cache *cache,
		    short int *id,
		    int k);

int
ALIGN_compare_alignments_by_id(ALIGN_alignment *a,
			       ALIGN_alignment *b);



void
ALIGN_print_alignment_cache(ALIGN_alignment_cache *cache);

int
ALIGN_cmp_id_string(short int *a,
		    short int *b);

unsigned int
ALIGN_hash(short int *id);

void
ALIGN_add_hashnode(ALIGN_alignment_cache *cache,
		   ALIGN_cache_hashnode *node,
		   short int *id);

void
ALIGN_free_hashnode(ALIGN_cache_hashnode *hashnode);

void
ALIGN_del_hashnode(ALIGN_alignment_cache *cache,
		   ALIGN_cache_hashnode *node);

void
ALIGN_add_listnode(ALIGN_alignment_cache *cache,
		   ALIGN_cache_listnode *node);

void
ALIGN_del_listnode(ALIGN_alignment_cache *cache,
		   ALIGN_cache_listnode *node);

void 
ALIGN_add_cache_node(ALIGN_alignment_cache *cache,
		     ALIGN_alignment *alignment,
		     ALIGN_dag *dag,
		     unsigned int n,
		     unsigned int k,
		     short int *id);

void
ALIGN_cache_promote(ALIGN_alignment_cache *cache,
		    ALIGN_cache_hashnode *hashnode,
		    ALIGN_cache_listnode *listnode);

ALIGN_dag *
ALIGN_build_dag(ALIGN_matrix *matrix,
		ALIGN_alignment *a,
		ALIGN_alignment *b);

ALIGN_dagnode *
ALIGN_recurse_build_dag(ALIGN_matrix *matrix,
			unsigned int i,
			unsigned int j);

void
ALIGN_print_dag(ALIGN_dag *dag);

void
ALIGN_recurse_print_dag(ALIGN_dagnode *node);

void
ALIGN_free_dag(ALIGN_dag *dag);

void
ALIGN_recurse_free_dag(ALIGN_dagnode *dagnode);

void
ALIGN_clear_visited(ALIGN_dag *dag);

void
ALIGN_recurse_clear_visited(ALIGN_dagnode *dagnode);

		
#endif /* _INC_CACHE_H_ */







