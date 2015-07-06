/*
 * cache.c
 *
 * $Id: cache.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2005 Luke Sheneman.  All rights reserved.
 *
 *
 * Luke Sheneman
 * 04/01/05
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "align.h"
#include "cache.h"


int global_node_id;
ALIGN_alignment_cache *global_cache_ptr = NULL;



/*
 * ALIGN_init_alignment_cache() -- 
 *
 * Allocate our LRU alignment cache
 *
 */
ALIGN_alignment_cache *
ALIGN_init_alignment_cache(void) {

  ALIGN_alignment_cache *cache;
  int i;

  cache = (ALIGN_alignment_cache *)calloc(1, sizeof(ALIGN_alignment_cache));

  cache->hash    = (ALIGN_cache_hashnode **)calloc(ALIGN_HASHSIZE, sizeof(ALIGN_cache_hashnode *));
  for(i=0;i<ALIGN_HASHSIZE;i++) {
    cache->hash[i] = NULL;
  }

  cache->list    = NULL;
  cache->listend = NULL;

  cache->size = 0;
  
  global_cache_ptr = cache;

  return(cache);
}







/*
 * ALIGN_merge_alignment_ids() - 
 *
 * Merge alignment IDs in a canonical way 
 *
 */
void
ALIGN_merge_alignment_ids(ALIGN_alignment *dest,
			  ALIGN_alignment *source1,
			  ALIGN_alignment *source2) {

  if(dest->id) {
    free(dest->id);
  }

  dest->id = ALIGN_compute_alignment_id(source1, source2);
  
  return;
}






/*
 * ALIGN_compute_alignment_id() - 
 *
 * Allocates a buffer big enough and stuffs it 
 * full with the canonical alignment representation.
 *
 */
short int *
ALIGN_compute_alignment_id(ALIGN_alignment *a,
			   ALIGN_alignment *b) {

  short int *result;
  int i, j, max;

  max = (a->k+b->k)*3;

  /* allocate the space for the id and initialize it */
  result = (short int *)malloc(max*sizeof(short int));
  memset(result, ALIGN_ID_EMPTY, max*sizeof(short int));

  result[0] = ALIGN_L_PAREN;

  j = 1;
  if(ALIGN_determine_id_order(a->id, b->id) < 0) {  // a < b

    i = 0;
    while(a->id[i] != ALIGN_ID_EMPTY) {
      result[j++] = a->id[i++];
    }

    i = 0;
    while(b->id[i] != ALIGN_ID_EMPTY) {
      result[j++] = b->id[i++];
    }
    
    
  } else {                                  // b < a

    i = 0;
    while(b->id[i] != ALIGN_ID_EMPTY) {
      result[j++] = b->id[i++];
    }

    i = 0;
    while(a->id[i] != ALIGN_ID_EMPTY) {
      result[j++] = a->id[i++];
    }

  }
  
  result[j] = ALIGN_R_PAREN;

  return(result);
}








/*
 * ALIGN_determine_id_order() -
 *
 *  1 if a > b
 * -1 if a < b
 *
 */
int
ALIGN_determine_id_order(short int *a,
			 short int *b) {

  int i;
  short int aval, bval;
  
  i = 0;
  while(1) {
    if(a[i] != ALIGN_L_PAREN &&
       a[i] != ALIGN_R_PAREN) {
      aval = a[i];
      break;
    }
    i++;
  }

  i = 0;
  while(1) {
    if(b[i] != ALIGN_L_PAREN &&
       b[i] != ALIGN_R_PAREN) {
      bval = b[i];
      break;
    }
    i++;
  }
  
  if(aval > bval) {
    return(1);
  } else if (aval < bval) {
    return(-1);
  } else {
    fprintf(stderr, "ASSERTION IN ALIGN_determine_id_order()\n");
    exit(-1);
  }
}






/*
 * ALIGN_print_alignment_id() -
 *
 * Print the alignment ID string here
 *
 */
void
ALIGN_print_alignment_id(short int *id) {

  int i = 0;
  
  while(id[i] != ALIGN_ID_EMPTY) {
    if(id[i] == ALIGN_L_PAREN) {
      printf("(");
    } else if(id[i] == ALIGN_R_PAREN) {
      printf(")");
    } else {
      printf("%d", id[i]);
    }
    
    i++;
  }

  return;
}





/*
 * ALIGN_find_hashnode() - 
 *
 * Find a hashnode in the cache, given only an alignment ID
 *
 */
ALIGN_cache_hashnode *
ALIGN_find_hashnode(ALIGN_alignment_cache *cache,
		    short int *id,
		    int k) {

  ALIGN_cache_hashnode *node;
  unsigned int hashindex;
  
  hashindex = ALIGN_hash(id);
  
  node = cache->hash[hashindex];
  
  while(node) {
    
    if(node->k == k) {
      if(ALIGN_cmp_id_string(node->id, id)) {
	return(node);
      }
    }

    node = node->next;
  }

  return(NULL);
}




/*
 * ALIGN_compare_alignments_by_id() -
 *
 * Compare two alignments by first comparing their
 * dimensions and then comparing their canonized id strings
 *
 */
int
ALIGN_compare_alignments_by_id(ALIGN_alignment *a,
			       ALIGN_alignment *b) {
  if( (a->n != b->n) ) {
    return(0);
  }
  
  if( (a->k != b->k) ) {
    return(0);
  }
  
  if(ALIGN_cmp_id_string(a->id, b->id)) {
    return(1);
  } else {
    return(0);
  }
}



/*
 * ALIGN_compare_hashnodes_by_id() -
 *
 * Compare two hashnodes by first comparing their
 * dimensions and then comparing their canonized id strings
 *
 */
int
ALIGN_compare_hashnodes_by_id(ALIGN_cache_hashnode *a,
			      ALIGN_cache_hashnode *b) {
  if( (a->n != b->n) ) {
    return(0);
  }
  
  if( (a->k != b->k) ) {
    return(0);
  }
  
  if(ALIGN_cmp_id_string(a->id, b->id)) {
    return(1);
  } else {
    return(0);
  }
}





/*
 * ALIGN_free_alignment_cache() --
 *
 * Free the cache nodes
 *
 */
void
ALIGN_free_alignment_cache(void) {

  ALIGN_alignment_cache *cache;
  ALIGN_cache_listnode *node, *t;
  
  cache = global_cache_ptr;

  if(cache) {
    
    /* free the hash */
    free(cache->hash);

    /* free the cache list and all of the hash pointers */
    node = cache->list;
    while(node) {

      /* free the hashnode and everything in it */
      if(node->hashnode) {
	ALIGN_free_hashnode(node->hashnode);
      }

      t = node->next;
      free(node);
      node = t;
    }

    free(cache);
  }
  
  return;
}


			      






/*
 * ALIGN_print_alignment_cache() - 
 *
 * Print details of the cache
 * 
 */
void
ALIGN_print_alignment_cache(ALIGN_alignment_cache *cache) {
  
  int i;
  ALIGN_cache_listnode *listnode;
  ALIGN_cache_hashnode *hashnode;

  
  printf("\n\n");
  printf("**************************************\n");
  printf("ALIGNMENT CACHE: maxsize = %d, current_size = %d\n", ALIGN_CACHESIZE, cache->size);
  
  /* print the list */
  printf("CACHE LIST:  START = %p, END = %p, size = %d\n", cache->list, cache->listend, cache->size);

  i = 0;
  listnode = cache->list;
  while(listnode) {
    printf("NODE #%d: (%p)\n", i, listnode);
    printf("  PREV:   (%p)\n", listnode->prev);
    printf("  NEXT:   (%p)\n", listnode->next);
    printf("  HASH:   (%p)\n", listnode->hashnode);  
    
    listnode = listnode->next;
    i++;
  }
  
  printf("\n---\n");

  printf("HASH:  Hash size = %d\n", ALIGN_HASHSIZE);
  for(i=0;i<ALIGN_HASHSIZE;i++) {
    printf("SLOT #%d:\n", i);
    hashnode = cache->hash[i];
    if(!hashnode) {
      printf(" [NULL]\n");
    }
    while(hashnode) {
      printf("  HASHNODE: (%p)\n", hashnode);
      printf("    PREV: (%p)\n", hashnode->prev);
      printf("    NEXT: (%p)\n", hashnode->next);
      printf("    LIST: (%p)\n", hashnode->listnode);
      printf("    ALGN: (%p)\n", hashnode->alignment);
      printf("     DAG: (%p)\n", hashnode->dag);

      hashnode = hashnode->next;
    }
  }

  printf("**************************************\n");
  
  printf("\n\n\n");

  return;
}










/*
 * ALIGN_cmp_id_string() - 
 *
 * Compare two raw id strings of the same length 
 *
 */
int
ALIGN_cmp_id_string(short int *a,
		    short int *b) {

  int i;

  i = 0;
  while(a[i] != ALIGN_ID_EMPTY &&
	b[i] != ALIGN_ID_EMPTY) {

    if(a[i] != b[i]) {
      return(0);
    }

    i++;
  }
  
  /* make sure one string is not continuing while the other is empty */
  if(a[i] != b[i]) {
    return(0);
  }
  
  return(1);  /* return same */
}







/*
 * ALIGN_hash() -- 
 *
 * A hash function for hashing canonical alignment ids 
 *
 */
unsigned int
ALIGN_hash(short int *id) {
  
  unsigned int hash = 0;
  int i = 0;
  
  while(id[i] != ALIGN_ID_EMPTY) {
    hash = id[i] + (hash << 6) + (hash << 16) - hash;
    i++;
  }
  
  hash = hash % ALIGN_HASHSIZE;

  return(hash);
}





/*
 * ALIGN_add_hashnode() - 
 *
 * Add a hashnode to the hash in the cache 
 *
 */
void
ALIGN_add_hashnode(ALIGN_alignment_cache *cache,
		   ALIGN_cache_hashnode *hashnode,
		   short int *id) {

  unsigned int hashkey;

  hashkey = ALIGN_hash(id);

  if(cache->hash[hashkey]) {
    cache->hash[hashkey]->prev = hashnode;
  }

  hashnode->prev = NULL;
  hashnode->next = cache->hash[hashkey];
  cache->hash[hashkey] = hashnode;
  
  return;
}




/*
 * ALIGN_free_hashnode() - 
 *
 * Free the hashnode and internal items
 *
 */
void
ALIGN_free_hashnode(ALIGN_cache_hashnode *hashnode) {

  if(hashnode) {

    if(hashnode->alignment) {
      ALIGN_free_alignment(hashnode->alignment);
    }

    if(hashnode->dag) {
      ALIGN_free_dag(hashnode->dag);
    }
      
    if(hashnode->id) {
      free(hashnode->id);
    }

    free(hashnode);  /* free the associated hashnode */
  }

  return;
}





/*
 * ALIGN_del_hashnode() - 
 *
 * Delete the hashnode from the hash in the cache
 *
 */
void
ALIGN_del_hashnode(ALIGN_alignment_cache *cache,
		   ALIGN_cache_hashnode *node) {

  unsigned int hashkey;

  
  hashkey = ALIGN_hash(node->id);

  if(node->next) {
    node->next->prev = node->prev;
  }
  
  if(node->prev) {
    node->prev->next = node->next;
  }

  if(cache->hash[hashkey] == node) {
    cache->hash[hashkey] = NULL;
  }
  
  ALIGN_free_hashnode(node);

  return;
}






/*
 * ALIGN_add_listnode() - 
 *
 * Add a listnode to the list in the cache
 *
 */
void
ALIGN_add_listnode(ALIGN_alignment_cache *cache,
		   ALIGN_cache_listnode *node) {
  
  node->prev = NULL;
  node->next = cache->list;
  if(node->next) {
    node->next->prev = node;
  }

  if(cache->list == NULL) {
    cache->listend = node;
  }

  cache->list = node;
  
  /* if the cache is filled up */
  if(cache->size == ALIGN_CACHESIZE) {
    ALIGN_del_listnode(cache, cache->listend);
  } else {
    cache->size++;
    if(cache->size == ALIGN_CACHESIZE) {
      printf("alignment cache is full\n");
    }
  }

  return;
}





/*
 * ALIGN_del_listnode() -
 *
 * Delete a node from the list structure in the cache
 *
 */
void
ALIGN_del_listnode(ALIGN_alignment_cache *cache,
		   ALIGN_cache_listnode *node) {

  if(node->prev) {
    node->prev->next = node->next;
    if(cache->listend == node) {
      cache->listend = node->prev;
    }
  }
  
  if(node->next) {
    node->next->prev = node->prev;
  }
  
  if(cache->list == node) {
    cache->list    = NULL;
    cache->listend = NULL;
  }
  
  if(node->hashnode) {
    ALIGN_del_hashnode(cache, node->hashnode);
  }
  
  free(node);
  
  return;
}




/*
 * ALIGN_add_cache_node() -
 *
 * Add a node w/alignment to the cache
 *
 */
void
ALIGN_add_cache_node(ALIGN_alignment_cache *cache,
		     ALIGN_alignment *alignment,
		     ALIGN_dag *dag,
		     unsigned int n,
		     unsigned int k,
		     short int *id) {

  ALIGN_cache_hashnode *hashnode;
  ALIGN_cache_listnode *listnode;

  /* do some debugging checks.  one or other, but not both should be set */
  if(!alignment && !dag) {
    printf("In ALIGN_add_cache_node() and alignment and DAG are both NULL\n");
    exit(-1);
  } else if(alignment && dag) {
    printf("In ALIGN_add_cache_node() and alignment and DAG are both set\n");
    exit(-1);
  }
  
  /* allocate our nodes here */
  hashnode = (ALIGN_cache_hashnode *)calloc(1, sizeof(ALIGN_cache_hashnode));
  listnode = (ALIGN_cache_listnode *)calloc(1, sizeof(ALIGN_cache_listnode));

  /* cross reference hashnodes and listnodes */
  hashnode->listnode  = listnode;
  listnode->hashnode  = hashnode;

  /* duplicate our alignment */
  if(alignment) {
    hashnode->alignment = ALIGN_dup_alignment(alignment);
  }

  /* assign the dag to the hashnode */
  hashnode->dag = dag;
  
  /* assign some specific values to the hashnode */
  hashnode->n = n;
  hashnode->k = k;
  hashnode->id = (short int *)calloc(k*3, sizeof(short int));
  memcpy(hashnode->id, id,((k*3)*sizeof(short int)));
  
  /* handle listnode by prepending it to the list */
  ALIGN_add_listnode(cache, listnode);

  /* handle hashnode by hashing alignnment id and prepending hashnode to cache's hash */
  ALIGN_add_hashnode(cache, hashnode, id);
  
  return;
}






/*
 * ALIGN_cache_promote() - 
 *
 * Promote the specified hashnode to the beginning of the hashlist
 * Promote the specified listnode to the beginning of the list 
 *
 */
void
ALIGN_cache_promote(ALIGN_alignment_cache *cache,
		    ALIGN_cache_hashnode *hashnode,
		    ALIGN_cache_listnode *listnode) {

  short int *id;
  unsigned int hashkey;
  ALIGN_cache_hashnode *hash_head;
  ALIGN_cache_listnode *list_head;
  
  /* take care of the hashnode */  
  id = hashnode->id;
  hashkey = ALIGN_hash(id);
  hash_head = cache->hash[hashkey];
  if(hashnode != hash_head) {
    
    /* remove hashnode from the list */
    if(hashnode->prev) {
      hashnode->prev->next = hashnode->next;
    }
    if(hashnode->next) {
      hashnode->next->prev = hashnode->prev;
    }

    /* put hashnode at the beginning of the list */
    hash_head->prev = hashnode;
    hashnode->next = hash_head;
    hashnode->prev = NULL;
    
    cache->hash[hashkey] = hashnode;
  }
  
  
  
  /* take care of the list */
  list_head = cache->list;
  if(listnode != list_head) {

    /* remove listnode from the list */
    if(listnode->prev) {
      listnode->prev->next = listnode->next;
    }
    if(listnode->next) {
      listnode->next->prev = listnode->prev;
    }
    
    /* put the listnode at the beginning of the list */
    list_head->prev = listnode;
    listnode->next = list_head;
    listnode->prev = NULL;

    cache->list = listnode;
  }

  return;
}



/*** The DAG Functions ***/


/*
 * ALIGN_build_dag() - 
 *
 * Given a DP matrix, construct a DAG which encodes all 
 * optimal paths through the matrix
 *
 */
ALIGN_dag *
ALIGN_build_dag(ALIGN_matrix *matrix,
		ALIGN_alignment *a,
		ALIGN_alignment *b) {
  
  ALIGN_dag *dag;
  
  if(!matrix) {
    fprintf(stderr, "Matrix pointer was NULL in ALIGN_build_dag()\n");
    return(NULL);
  }
  
  global_node_id = 0;

  /* allocate the DAG here */
  dag = (ALIGN_dag *)calloc(1, sizeof(ALIGN_dag));
  
  dag->root = 
    ALIGN_recurse_build_dag(matrix, 
			    matrix->height-1, 
			    matrix->width-1);
  
  dag->maxulps = matrix->maxulps;
  
  dag->id = ALIGN_compute_alignment_id(a, b);
  dag->k = a->k + b->k;

  dag->root->refcount = 1;
  
  return(dag);
}




/*
 * ALIGN_recurse_build_dag() - 
 *
 * Recursively construct a directed acyclic graph (DAG) from 
 * a DP/alignment matrix
 *
 */
ALIGN_dagnode *
ALIGN_recurse_build_dag(ALIGN_matrix *matrix,
			unsigned int i,
			unsigned int j) {
  
  ALIGN_dagnode *retnode = NULL;
  ALIGN_cell *cell;
  char dir;
  
  cell = &(MATRIX_cell(matrix, i, j));
  dir = cell->dir;

  /*  
  printf("dir = %d\n", dir);
  printf("val = %lf\n", MATRIX_cell(matrix, i, j).value);
  printf("i=%d, j=%d\n", i, j);
  */

  /* allocate our retnode here */
  retnode = (ALIGN_dagnode *)calloc(1, sizeof(ALIGN_dagnode));
  retnode->npaths   = cell->npaths;
  retnode->value    = cell->value;
  retnode->visited  = 0;
  retnode->refcount = 0;

  if(i==0 && j==0) {

    retnode->horz = NULL;
    retnode->vert = NULL;
    retnode->diag = NULL;
    
    retnode->node_id  = global_node_id++;

    cell->dagnode = retnode;

    return(retnode);

  }

  retnode->node_id = global_node_id++;
  
  /* check horizontal direction */
  if(GETH(dir)) {
    if(MATRIX_cell(matrix, i, j-1).dagnode) {
      retnode->horz = MATRIX_cell(matrix, i, j-1).dagnode;
    } else {
      retnode->horz = ALIGN_recurse_build_dag(matrix, i, j-1);
      MATRIX_cell(matrix, i, j-1).dagnode = retnode->horz;
    }
    
    retnode->horz->refcount++;
  }
  
  if(GETV(dir)) {
    if(MATRIX_cell(matrix, i-1, j).dagnode) {
      retnode->vert = MATRIX_cell(matrix, i-1, j).dagnode;
    } else {
      retnode->vert = ALIGN_recurse_build_dag(matrix, i-1, j);
      MATRIX_cell(matrix, i-1, j).dagnode = retnode->vert;
    }

    retnode->vert->refcount++;
  }
  
  if(GETD(dir)) {
    if(MATRIX_cell(matrix, i-1, j-1).dagnode) {
      retnode->diag = MATRIX_cell(matrix, i-1, j-1).dagnode;
    } else {
      retnode->diag = ALIGN_recurse_build_dag(matrix, i-1, j-1);
      MATRIX_cell(matrix, i-1, j-1).dagnode = retnode->diag;
    }

    retnode->diag->refcount++;
  }

  return(retnode);
}




/*
 * ALIGN_print_dag() - 
 *
 * Print an alignment DAG in a reasonable way
 * by outputting the graph in DOT format
 * for viewing with graph visualizer software
 * such as graphviz.
 *
 * http://www.graphviz.org
 *
 */
void
ALIGN_print_dag(ALIGN_dag *dag) {
  
  ALIGN_clear_visited(dag);
  
  printf("/*  Evalyn -- DAG output for DP matrix.  */\n");

  printf("\n");
  printf("digraph evalyn_dag {\n");
  
  if(dag && dag->root) {
    ALIGN_recurse_print_dag(dag->root);
  } else {
    printf("NO DAG\n");
  }

  printf("}\n\n");

  ALIGN_clear_visited(dag);
  
  return;
}




/*
 * ALIGN_recurse_print_dag() - 
 *
 * Recursively traverse the graph and print node ids
 *
 */
void
ALIGN_recurse_print_dag(ALIGN_dagnode *node) {
  
  if(node->visited == 0) {

    if(node->horz) {
      ALIGN_recurse_print_dag(node->horz);
      printf("%d -> %d ;\n", node->node_id, node->horz->node_id);
    }

    if(node->vert) {
      ALIGN_recurse_print_dag(node->vert);
      printf("%d -> %d ;\n", node->node_id, node->vert->node_id);
    }

    if(node->diag) {
      ALIGN_recurse_print_dag(node->diag);
      printf("%d -> %d ;\n", node->node_id, node->diag->node_id);
    }
  
    node->visited = 1;

  }

  return;
}




/*
 * ALIGN_clear_visited() - 
 *
 * Clear a DAG's visited flags
 *
 */
void
ALIGN_clear_visited(ALIGN_dag *dag) {

  ALIGN_recurse_clear_visited(dag->root);
  
  return;
}



/*
 * ALIGN_recurse_clear_visited() - 
 *
 * Recurse through a DAG clearing the visited flag  in
 * every node
 *
 */
void
ALIGN_recurse_clear_visited(ALIGN_dagnode *dagnode) {
  
  if(dagnode->visited == 0) {
    return;
  }
  
  if(dagnode->horz) {
    ALIGN_recurse_clear_visited(dagnode->horz);
  }

  if(dagnode->vert) {
    ALIGN_recurse_clear_visited(dagnode->vert);
  }

  if(dagnode->diag) {
    ALIGN_recurse_clear_visited(dagnode->diag);
  }

  dagnode->visited = 0;
  
  return;
}



/*
 * ALIGN_free_dag() - 
 *
 * Frees the given DAG
 *
 */
void
ALIGN_free_dag(ALIGN_dag *dag) {

  if(dag) {
    
    if(dag->root) {
      ALIGN_recurse_free_dag(dag->root);
    }
    
    if(dag->id) {
      free(dag->id);
    }

    free(dag);
  }

  return;
}
 


/*
 * ALIGN_recurse_free_diag() - 
 *
 * Recursively traverse the DAG, freeing nodes 
 *
 */
void
ALIGN_recurse_free_dag(ALIGN_dagnode *dagnode) {
  
  if(dagnode->horz) {
    ALIGN_recurse_free_dag(dagnode->horz);
    dagnode->horz = NULL;
  }
  
  if(dagnode->vert) {
    ALIGN_recurse_free_dag(dagnode->vert);
    dagnode->vert = NULL;
  }

  if(dagnode->diag) {
    ALIGN_recurse_free_dag(dagnode->diag);
    dagnode->diag = NULL;
  }
  
  dagnode->refcount--;
  if(dagnode->refcount == 0) {
    free(dagnode);
  }
  
  return;
}
