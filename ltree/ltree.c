/*
 * ltree.c
 *
 * $Id: ltree.c,v 1.3 2006/05/09 21:10:51 sheneman Exp $
 *
 *
 * Copyright (c) 2005 Luke Sheneman.  All rights reserved.
 *
 *
 * A GA for evolving trees using crossover and mutation with constraints, 
 * using progressive alignment scores based on trees as the fitness function 
 *
 *
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/param.h>
#include <unistd.h>
#include <values.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#include <getopt.h>

#include <nj_fasta.h>
#include <nj_clearcut.h>

#include <align.h>
#include <cache.h>
#include <msf.h>
#include <fasta.h>
#include <aln.h>


#include "ltree.h"
#include "dnd.h"
#include "parsimony.h"
#include "rnjhooks.h"








/*************************************
 *                                   *
 * BEGIN STATIC FUNCTION PROTOTYPES  *
 *                                   *
 *************************************/

static LTREE_treenode *
LTREE_find_cp(LTREE_treenode *node);

static void
LTREE_recursive_print_tree(LTREE_treenode *node,
			   ALIGN_seqset *seqset);

static void
LTREE_print_population(LTREE_population *population);

static void
LTREE_free_subtree(LTREE_treenode *node);

static double
LTREE_recursive_align_tree(LTREE_treenode *node,
			   ALIGN_scoring_system *ss);

static int
LTREE_bsearch(LTREE_population *population,
	      double fitness);

static int
LTREE_bsearch_compare(const void *a,
		      const void *b);

static void
LTREE_recursive_write_tree(LTREE_treenode *node,
			   ALIGN_seqset *seqset,
			   FILE *fp);

static void
LTREE_recursive_compress_tree(LTREE_treenode *node);

static unsigned int
LTREE_nodes_at_depth(LTREE_treenode *node,
		     unsigned int target_depth,
		     unsigned int current_depth);

static void
LTREE_get_vector(LTREE_treenode *node,
		 unsigned int target_depth,
		 unsigned int current_depth,
		 LTREE_treenode **vector);

static void
LTREE_print_fitnesses(LTREE_population *pop);

static void
LTREE_init_leaf_alignments(LTREE_tree *tree);

static void
LTREE_recursive_init_leaf_alignments(LTREE_treenode *node,
				     ALIGN_seqset *seqset);
static void
LTREE_handle_options(int argc, 
		     char *argv[], 
		     LTREE_cl_args *cl_args,
		     ALIGN_scoring_system *ss);

static void
LTREE_print_options(LTREE_cl_args *cl_args);

static void
LTREE_set_seed(LTREE_cl_args *cl_args);

static void
LTREE_auto_output_filename(LTREE_cl_args *cl_args);

static void
LTREE_assign_gap_penalties(ALIGN_scoring_system *ss, 
			   LTREE_cl_args *cl_args);

static int
LTREE_autodetect_input_format(LTREE_cl_args *cl_args);

static int
LTREE_determine_input_type(LTREE_cl_args *cl_args);

static int 
LTREE_check_population_order(LTREE_population *population);

static void
LTREE_build_taxa_list(LTREE_treenode *node);

static void
LTREE_add_taxa_node(LTREE_taxa_list *head,
		    LTREE_taxa_list *node);

static void
LTREE_print_taxa_list(LTREE_taxa_list *head);

static void
LTREE_free_taxa_list(LTREE_taxa_list *head);

static double
LTREE_align_node(LTREE_treenode *node,
		 ALIGN_scoring_system *ss);

static int
LTREE_compare_fitness(const void *a,
		      const void *b);

static LTREE_treenode *
LTREE_node_crossover(LTREE_treenode *parent1,
		     LTREE_treenode *parent2);

static LTREE_tree *
LTREE_dup_tree(LTREE_tree *tree);

static LTREE_treenode *
LTREE_dup_treenodes(LTREE_treenode *node);

static void
LTREE_tag_tree(LTREE_tree *tree,
	       LTREE_taxa_list *taxa_list);

static void
LTREE_recursive_tag_nodes(LTREE_treenode *node,
			  LTREE_taxa_list *taxa_list);

static void
LTREE_remove_tagged_nodes(LTREE_treenode *node);

static unsigned int
LTREE_in_taxa_list(LTREE_taxa_list *taxa_list,
		   short int seq_id);

static void
LTREE_clear_internal_alignments(LTREE_tree *tree);

static void
LTREE_recursive_clear_internal_alignments(LTREE_treenode *node);




/***********************************
 *                                 *
 * END STATIC FUNCTION PROTOTYPES  *
 *                                 *
 ***********************************/





/* GLOBAL VARIABLES */
unsigned int LTREE_vector_index;
unsigned long int LTREE_global_node_count;
LTREE_taxa_list *LTREE_global_taxa_list;





/*
 * main() - 
 *
 */
int
main(int argc,
     char *argv[]) {

  LTREE_cl_args cl_args;
  LTREE_population *pop;
  LTREE_tree *best_tree;
  MT_state old_prng_state;

  ALIGN_scoring_system ss;
  ALIGN_seqset *seqset;

  short int type;
  double fitness;


  /* here, we handle all of the options passed to evalyn */
  LTREE_handle_options(argc, argv, &cl_args, &ss);

  /* print the settings */
  LTREE_print_options(&cl_args);

  /* parse the substitution matrix file here */
  ss.submat = 
    ALIGN_read_submat(cl_args.matrix_filename);
  if(!ss.submat) {
    fprintf(stderr, "EVALYN: Fatal error parsing substitution matrix file\n");
    exit(-1);
  }

  /* for DEBUGGING */
  // ALIGN_print_submat(ss.submat);
  
  /* allocate our set of sequences */
  seqset = ALIGN_alloc_seqset(ss.submat->symtab);
  
  /* treat the input file properly (FASTA?  MSF?) */
  switch(cl_args.infiletype) {
  case LTREE_MSF:

    /* read the MSF file and get the type specified in the file (DNA, Protein, etc) */
    type = 
      ALIGN_read_MSF(seqset, 
		     cl_args.input_filename);
    if(type == -1) {
      fprintf(stderr, "EVALYN: Failed to read MSF file %s\n", cl_args.input_filename);
      exit(-1);
    }

    break;

  case LTREE_FASTA:
    
    /* read the specified file and ascertain the type */
    type = 
      ALIGN_read_FASTA(seqset, 
		       cl_args.input_filename);
    if(type == -1) {
      fprintf(stderr, "EVALYN: Failed to read FASTA file %s\n", cl_args.input_filename);
      exit(-1);
    }

    break;

  default:
    printf("EVALYN: Unknown input file format\n");
    exit(-1);
  }
  

  /* 
   * Lets make sure that our substitution matrix type (DNA/PROTEIN) match
   * the type inferred from reading the sequences, if user was not explicit
   * about types on the command-line.
   */
  if(!cl_args.dna_flag && !cl_args.protein_flag) {
    if(type != ss.submat->matrix_type) {
      fprintf(stderr, "Error: Discrepancy between the sequence type and the sub. matrix type.\n");
      if(ss.submat->matrix_type == ALIGN_PROTEIN) {
	if(type == ALIGN_DNA) {
	  printf("Matrix is PROTEIN.  Sequences are DNA.\n");
	}
      } else if(ss.submat->matrix_type == ALIGN_DNA) {
	if(type == ALIGN_PROTEIN) {
	  printf("Matrix is DNA.  Sequences are PROTEIN.\n");
	}
      }
      exit(-1);
    }
  }
  
  /* Iterate and evolve (GA) */
  pop = 
    LTREE_evolve(&cl_args,
		 seqset,
		 &ss);

  best_tree = pop->trees[0];

  /* print the best tree to stdout */
  printf("\n");
  LTREE_print_tree(best_tree); 
  printf("\n\n");

  /* write out the tree */
  LTREE_write_tree_DND(best_tree, cl_args.tree_filename);
  
  /* 
     We discarded all alignments after computing them, 
     so now re-compute the final alignment using the same 
     PRNG state that we used when we computed the alignment
     for the best tree.
  */
  MT_get_state(&old_prng_state);    /* get current PRNG state */
  MT_set_state(&(pop->best_prng_state));    /* set PRNG state */
  
  fitness = 
    LTREE_eval_fitness(seqset,
		       &ss,
		       best_tree,
		       SAVE_ALIGNMENT);
  
  MT_set_state(&old_prng_state);    /* restore PRNG state */

  /* print the alignment to stdout here */
  ALIGN_print_alignment(best_tree->root->alignment, seqset);

  /* write the best alignment to disk here */
  LTREE_write_alignment(best_tree->root->alignment,
			seqset,
			&cl_args);
  
  printf("\nSum-of-Pairs Score (SPS): %lf\n\n", fitness);
  
  LTREE_dump_population("popdump.txt", pop, -1, NULL, NULL, NULL);
  
  /* free everything to make valgrind extra happy */
  LTREE_free_population(pop);

  ALIGN_free_seqset(seqset);
  ALIGN_free_submat(ss.submat);
  ALIGN_free_matrix();
  ALIGN_free_alignment_cache();
  
  exit(0);
}





/*
 * LTREE_align_node() -
 *
 * when enough information is at a particular node, align that node
 *
 */
static double
LTREE_align_node(LTREE_treenode *node,
		 ALIGN_scoring_system *ss) {
  
  /* align the left alignment to the right alignment */
  node->alignment = 
    ALIGN_alignment_alignment(ss,
			      node->left->alignment,
			      node->right->alignment);

  return(node->alignment->score);
}





/*
 * LTREE_pop_stats() - 
 * 
 * Collect fitness stats for the population 
 *
 */
void
LTREE_pop_stats(LTREE_population *population,
		double *best_fitness,
		double *mean_fitness) {
  
  unsigned int i;

  *mean_fitness = 0.0;
  *best_fitness = -(DBL_MAX);

  for(i=0;i<population->size;i++) {
    
    *mean_fitness += population->trees[i]->fitness;

    if(population->trees[i]->fitness > *best_fitness) {
      *best_fitness = population->trees[i]->fitness;
    }
  }
  
  /* compute the mean fitness here */
  *mean_fitness /= (double)(population->size);
  
  return;
}




/*
 * LTREE_bsearch_compare()
 *
 * The comparison routine used by bsearch()
 *
 */
static int
LTREE_bsearch_compare(const void *a,
		      const void *b) {

  double fitness;
  LTREE_tree **tree, **next;


  fitness = *((double *)a);
  tree = (LTREE_tree **)b;
  next = (LTREE_tree **)((int)tree+(int)sizeof(LTREE_tree *));

  /* if fitnesses are equal */
  if(fitness == (*tree)->fitness) {
    return(0);
  }

  if( (fitness < (*tree)->fitness) &&
      (fitness > (*next)->fitness ) ) {
    return(0);
  }

  if(fitness > (*tree)->fitness) {
    return(-1);
  } else {
    return(1);
  }
}






/*
 * LTREE_bsearch() - 
 *
 * Perform a binary search to look for the proper 
 * insertion point
 *
 */
static int
LTREE_bsearch(LTREE_population *population,
	      double fitness) {

  LTREE_tree **tree;

  if( (fitness < population->trees[population->size-1]->fitness) ||
      (fitness == population->trees[population->size-1]->fitness) ) {
    return(population->size);
  } else if (fitness > population->trees[0]->fitness) {
    return(0);
  } else if ( (fitness > population->trees[0]->fitness) ) {
    return(1);
  } else {

    tree = 
      bsearch((void *)&fitness,
	      (void *)&(population->trees[0]),
	      population->size+1,
	      sizeof(LTREE_tree **),
	      LTREE_bsearch_compare);

    /* pointer arithmetic to determine the index into the tree population */
    if(fitness < (*tree)->fitness) {
      return(((int)tree-(int)(population->trees))/sizeof(LTREE_tree *)+1);
    } else {
      return(((int)tree-(int)(population->trees))/sizeof(LTREE_tree *));
    }
  }
}





/*
 * LTREE_check_population_order()
 *
 * For debugging.  Check the order of the individuals in the population
 * to make sure that they are sorted post-insertion.  They should 
 * always be sorted in descending order after an insertion.
 *
 * RETURNS:
 * --------
 *
 *   0 - Ordered
 *   1 - Disordered
 *
 */
static int
LTREE_check_population_order(LTREE_population *population) {

  int i;
  double prev;
  
  prev = population->trees[0]->fitness;
  for(i=1;i<population->size;i++) {
    if(prev < population->trees[i]->fitness) {

      printf("In LTREE_check_population_order() - Item %d (%lf) is smaller than its successor %d (%lf)\n", 
	     i-1, population->trees[i-1]->fitness,
	     i, population->trees[i]->fitness);

      return(1);
    }
    
    prev = population->trees[i]->fitness;
  }
 
  return(0);
}






/*
 * LTREE_update_population() -
 *
 * Update the population fitness table here
 *
 * This involves performing a binary search to locate the
 * point to insert the new tree, and doing some memmove()
 * calls to move things around to replace the tree marked
 * for removal via our selection routine
 *
 */
void
LTREE_update_population(LTREE_population *population,
			LTREE_tree *child_tree,
			unsigned short int child,
			double fitness) {
  
  int insertion_point;

  /* perform binary search to look for an insertion point in 
     the population which will maintain ordering
  */
  insertion_point = 
    LTREE_bsearch(population, fitness);

  /* free the child tree which will now be replaced */
  LTREE_free_tree(population->trees[child]);

  /* 
     The simple (but probably rare) case where the child to be 
     replaced is where the new tree needs to be to keep the 
     population sorted 
  */
  if(insertion_point == child) {

    population->trees[insertion_point] = child_tree;

  } else if(insertion_point < child) {

    /*
     * The case where the insertion point is higher up than 
     * the child to be replaced
     */
    memmove(&population->trees[insertion_point+1],
	    &population->trees[insertion_point],
	    (child-insertion_point)*sizeof(LTREE_tree *));

    population->trees[insertion_point] = child_tree;

  } else {

    memmove(&population->trees[child],
	    &population->trees[child+1],
	    (insertion_point-child-1)*sizeof(LTREE_tree *));

    population->trees[insertion_point-1] = child_tree;

  }
  
  return;
}






/*
 *
 * LTREE_evolve() - 
 *
 * The main GA loop
 *
 */
LTREE_population *
LTREE_evolve(LTREE_cl_args *cl_args,
	     ALIGN_seqset *seqset,
	     ALIGN_scoring_system *ss) {

  LTREE_population *pop;

  LTREE_tree *child_tree;
  LTREE_tree *rnj_tree;
  MT_state prng_state, best_prng_state;
  double fitness, best_fitness, mean_fitness, prob, bfit;
  double stationarity_fitness;
  int i, maxiter, stationarity_count;
  unsigned short int parent1, parent2, child;

  LTREE_population *neighborhood;  /* blah */

  /* initialize our entire population */
  pop = LTREE_init_population(seqset, cl_args->population_size);


#ifdef NEIGHBORHOOD_EXPERIMENT
  /****************************************************/
  printf("Finding LXO Neighborhood for the following tree:\n");
  LTREE_print_tree(pop->trees[0]);
  printf("\n\n");

  printf("TREE DESCRIPTION:\n");
  LTREE_print_tree_desc(LTREE_build_tree_desc(pop->trees[0]));

  neighborhood = LTREE_lxo_step(pop->trees[0]);
  exit(0);
  /****************************************************/

  for(i=0;i<300;i++) {
    LTREE_print_tree(pop->trees[i]);
    printf(";\n");
  }
  exit(0);
#endif /* NEIGHBORHOOD EXPERIMENT */

  /* compute the fitness for our entire population */
  printf("Computing initial fitness scores for population...\n");
  bfit = -DBL_MAX;
  for(i=0;i<cl_args->population_size;i++) {

    /* record PRNG state */
    MT_get_state(&prng_state);

    fitness = 
      LTREE_eval_fitness(seqset,
			 ss,
			 pop->trees[i],
			 PURGE_ALIGNMENT);

    if(fitness > bfit) {
      /* if this is the best fitness so far, save of PRNG state */
      bfit = fitness;
      MT_copy_state(&best_prng_state, &prng_state);
    }

    printf("%c%c%c%c%c%c%0.1f%%", 8,8,8,8,8,8,(double)i/(double)(cl_args->population_size)*100.0);

    fflush(stdout);
  }
  printf("%c%c%c%c%c%c%0.1f%%\n\n", 8,8,8,8,8,8,100.0);
  
  /* if use specified using RNJ tree, then inject that into the population */
  if(cl_args->rnj_flag) {

    rnj_tree = LTREE_create_rnj_tree(ss, pop, seqset);
    fitness = LTREE_eval_fitness(seqset,
				 ss,
				 rnj_tree,
				 PURGE_ALIGNMENT);

    // INSERT RNJ TREE INTO POPULATION HERE!!  sheneman
  }
  
  /* perform the initial sort of the population */
  LTREE_sort_population(pop);

  /* determine the initial best */
  stationarity_fitness = pop->trees[0]->fitness;

  LTREE_pop_stats(pop,
		  &best_fitness,
		  &mean_fitness);

  printf("INITIAL BEST:  %lf\n", best_fitness);
  printf("INITIAL WORST: %lf\n", pop->trees[pop->size-1]->fitness);
  printf("INITIAL MEAN:  %lf\n", mean_fitness);
  if(cl_args->rnj_flag) {
    printf("RNJ TREE FITNESS: %lf\n", fitness);
  }

  if(cl_args->converge_flag) {
    maxiter = INT_MAX;
  } else {
    maxiter = cl_args->iterations;
  }

  stationarity_count = 0;  
  for(i=0;i<maxiter;i++) {    

    /* print out results every 10 iterations */
    if(i%10 == 0) {

      /* get stats on the population */
      LTREE_pop_stats(pop,
		      &best_fitness,
		      &mean_fitness);

      printf("ITER: %d, BEST = %lf, MEAN = %lf\n", i, best_fitness, mean_fitness); 
      fflush(stdout);
    }
    
    /* select some trees here for crossover */
    LTREE_selection(pop,
		    &parent1,
		    &parent2,
		    &child);

    /* perform crossover here */
    child_tree = 
      LTREE_crossover(pop->trees[parent1],
		      pop->trees[parent2]);

    /* perform mutation here */
    prob = MT_get_real1();
    if(prob < cl_args->mutation_rate) {

      child_tree =
	LTREE_mutate(child_tree);

    }

    MT_get_state(&prng_state);

    /* evaluate the child's fitness */
    fitness = 
      LTREE_eval_fitness(seqset, 
			 ss,
			 child_tree,
			 PURGE_ALIGNMENT);
    
    if(fitness > bfit) {
      bfit = fitness;
      MT_copy_state(&best_prng_state, &prng_state);
    }
    
    /* put the child in the right place to keep things in order */
    LTREE_update_population(pop, child_tree, child, fitness);

    /* some debugging */
    LTREE_dump_population("popdump.txt", pop, i, pop->trees[parent1], pop->trees[parent2], child_tree);
    
    
    /* track the best fitness and how many iterations we go without it changing */
    if(cl_args->converge_flag) {
      if(fitness > stationarity_fitness) {
	stationarity_fitness = fitness;
	stationarity_count = 0;
      } else {

	stationarity_count++;

	if(stationarity_count >= cl_args->stationarity) {
	  printf("Reached stationarity (%d static steps) in %d iterations.\n",
		 cl_args->stationarity,
		 i+1);
	  break;
	}
      }
    }

  }

  /* Record the PRNG state associated with the best tree in the pop */
  MT_copy_state(&(pop->best_prng_state), 
		&best_prng_state);

  return(pop);
}












/*
 * LTREE_print_fitnesses() -
 *
 * Print all of the fitnesses here
 *
 */
static void
LTREE_print_fitnesses(LTREE_population *pop) {
  
  unsigned short int i;

  printf("\nPOPULATION FITNESS VALUES\n\n");
  for(i=0;i<pop->size;i++) {
    printf("%d> -- %lf\n", i, pop->trees[i]->fitness);
  }
  
  printf("\n");
  
  return;
}




/*
 * LTREE_mutate() -
 *
 * Mutate a tree 
 *
 */
LTREE_tree *
LTREE_mutate(LTREE_tree *child) {

  LTREE_tree *ret_tree;
  
 
  /* perform crossover in which we swap branches on the same tree */
  ret_tree = 
    LTREE_crossover(child,
		    child);


  return(ret_tree);
}



/*
 * TREE_selection() -
 *
 * Select two parents and a child
 *
 */
void
LTREE_selection(LTREE_population *population,
		unsigned short int *parent1,
		unsigned short int *parent2,
		unsigned short int *child) {
  
  /* we choose parents with high scores */

  /* choose two distinct parents for crossover */
  *parent1 = (unsigned short int)(MT_get_bounded_exp()*population->size);
  *parent2 = (unsigned short int)(MT_get_bounded_exp()*population->size);
  while(*parent2 == *parent1) {
    *parent2 = (unsigned short int)(MT_get_bounded_exp()*population->size);
  }

    
  /* choose a child which is neither parent and which also implements elitism */
  *child   = (unsigned short int)((1.0 - MT_get_bounded_exp())*population->size);
  while(*child == *parent1 || *child == *parent2 || *child == 0) {
    *child = (unsigned short int)((1.0 - MT_get_bounded_exp())*population->size);
  }

  return;
}





/*
 * LTREE_eval_fitness()
 * 
 * Take a chromosome and evaluate its fitness 
 *
 */
double
LTREE_eval_fitness(ALIGN_seqset *seqset,
		   ALIGN_scoring_system *ss,
		   LTREE_tree *tree,
		   char delete_alignment) {

  /* free and clear all of the internal alignments in the tree */
  LTREE_clear_internal_alignments(tree);
  
  LTREE_init_leaf_alignments(tree);
  
  /* recursively align the tree */
  tree->fitness = LTREE_align_tree(tree, ss);
  
  /* free and clear all of the internal alignments in the tree */
  if(delete_alignment == PURGE_ALIGNMENT) {
    LTREE_clear_internal_alignments(tree);
  }

  return(tree->fitness);
}




/*
 * LTREE_align_tree() - 
 *
 * Calls LTREE_recursive_align() to perform depth-first
 * recursive traversal of the guide tree to align 
 * left and right children of nodes in order to perform
 * the progressive alignment 
 *
 */
double
LTREE_align_tree(LTREE_tree *tree,
		 ALIGN_scoring_system *ss) {

  double fitness;

  fitness =
    LTREE_recursive_align_tree(tree->root, ss);

  return(fitness);
}




/*
 * ltree_recursive_align_tree() -
 *
 * A recursive function for computing fitness of an evaluation tree 
 *
 */
static double
LTREE_recursive_align_tree(LTREE_treenode *node,
			   ALIGN_scoring_system *ss) {

  double fitness = 0.0;

  if(node->left && node->right) {

    LTREE_recursive_align_tree(node->left,  ss);
    LTREE_recursive_align_tree(node->right, ss);

    fitness = LTREE_align_node(node, ss);
  }

  return(fitness);
}







/*
 * LTREE_init_population()
 * 
 * Initialize the population here
 *
 */
LTREE_population *
LTREE_init_population(ALIGN_seqset *seqset,
		      unsigned long int pop_size) {

  unsigned long int i;
  LTREE_population *pop;


  printf("Initializing Population:  pop->size = %ld, seq->num = %d\n", pop_size, seqset->num);

  pop = (LTREE_population *)calloc(1, sizeof(LTREE_population));
  pop->size = pop_size;
  
  pop->trees = (LTREE_tree **)calloc(pop->size+1, sizeof(LTREE_tree *));

  /* get the initial, random trees here */
  for(i=0;i<pop_size;i++) {
    pop->trees[i] = LTREE_init_tree(seqset);
  }
  
  pop->trees[pop_size] = NULL;

  return(pop);
}





 


/*
 * LTREE_compare_fitness() - 
 *
 * Actually compare the relative fitnesses of two individuals
 *
 */
static int 
LTREE_compare_fitness(const void *a,
		      const void *b) {
  
  LTREE_tree **x, **y;
  LTREE_tree *w, *z;
  
  x = (LTREE_tree **)a;
  y = (LTREE_tree **)b;
  
  w = *x;
  z = *y;
  
  if(w->fitness > z->fitness) {
    return -1;
  } else if (w->fitness < z->fitness) {
    return 1;
  } else {
    return 0;
  }
  
  return(0);
}
 


/*
 * LTREE_sort_population() -
 *
 * Sort the trees in the population based on fitness
 *
 */
void
LTREE_sort_population(LTREE_population *population) {

  qsort(population->trees,
	population->size, 
	sizeof(LTREE_tree *), 
	LTREE_compare_fitness);

  return;
}




/*
 * LTREE_crossover() - 
 *
 * A wrapper function for implementing crossover.  Here, we take first-class
 * tree function, and call another function which operates only on tree nodes.
 *
 */
LTREE_tree *
LTREE_crossover(LTREE_tree *parent1,
		LTREE_tree *parent2) {

  LTREE_tree *child;

  child = (LTREE_tree *)calloc(1, sizeof(LTREE_tree));
  child->seqset = parent1->seqset;
  
  child->root = 
    LTREE_node_crossover(parent1->root,
			 parent2->root);
  
  return(child);
}



/*
 * LTREE_node_crossover() - 
 *
 * Crossover the trees as Lewis has recommended
 *
 */
static LTREE_treenode *
LTREE_node_crossover(LTREE_treenode *parent1,
		     LTREE_treenode *parent2) {

  
  LTREE_treenode *cp, *cp_tree, *tmp, *insertion_point, *new_node, *op;
  LTREE_taxa_list *tlist_cp;


  /* find a crossover point in the first parent which is not the root */
  cp = LTREE_find_cp(parent1);  
  while(cp == parent1) {
    cp = LTREE_find_cp(parent1);
  }

  /* construct a duplicate subtree */
  cp_tree = LTREE_dup_treenodes(cp);

  /* construct a duplicate of parent2 */
  tmp = LTREE_dup_treenodes(parent2);

  /* construct a taxa list for all taxa under cp in child */
  LTREE_global_taxa_list = NULL;
  LTREE_build_taxa_list(cp_tree);
  tlist_cp = LTREE_global_taxa_list;

  /* remove all nodes from tmp which occur under cp in child */
  LTREE_recursive_tag_nodes(tmp, tlist_cp);
  LTREE_free_taxa_list(tlist_cp);
  LTREE_global_taxa_list = NULL;
  LTREE_remove_tagged_nodes(tmp);
  LTREE_recursive_compress_tree(tmp);

  /* pick insertion branch into tmp */
  insertion_point = LTREE_find_cp(tmp);  
  
  op = insertion_point->parent;

  /* allocate a new node which will join the trees */
  new_node              = (LTREE_treenode *)calloc(1, sizeof(LTREE_treenode));
  new_node->left        = insertion_point;
  new_node->right       = cp_tree;
  new_node->parent      = insertion_point->parent;
  new_node->alignment   = NULL;
  new_node->seq_id      = -1;
  new_node->remove_flag = 0;

  insertion_point->parent = new_node;
  
  if(insertion_point == tmp) {

    return(new_node);

  } else {

    if(op->left == insertion_point) {
      op->left = new_node;
    } else if(op->right == insertion_point) {
      op->right = new_node;
    } else {
      exit(0);
    }
  }

  return(tmp);
}





/*
 * LTREE_find_cp()
 *
 * Given a tree, determine a random crossover point
 *
 */
static LTREE_treenode *
LTREE_find_cp(LTREE_treenode *node) {
  
  LTREE_treenode *cp;
  unsigned long int size, depth, target_depth, nodes_at_depth;
  unsigned long int cp_node;
  LTREE_treenode **vector;
  
  depth = LTREE_tree_depth(node, 0);
  size  = LTREE_tree_size(node);

  /* randomly pick a target depth for crossover */
  //  target_depth = rand()%depth;
  target_depth = MT_get_range_uint32(depth);

  nodes_at_depth = LTREE_nodes_at_depth(node, target_depth, 0);
  
  vector = (LTREE_treenode **)calloc(nodes_at_depth, sizeof(LTREE_treenode *));

  // cp_node = rand()%nodes_at_depth;
  cp_node = MT_get_range_uint32(nodes_at_depth);

  LTREE_vector_index  = 0;
  LTREE_get_vector(node, target_depth, 0, vector);

  cp = vector[cp_node];
  free(vector);

  return(cp);
}






/*
 * LTREE_print_tree() - 
 *
 * Print the specified tree
 *
 */
void
LTREE_print_tree(LTREE_tree *tree) {

  LTREE_recursive_print_tree(tree->root,
			     tree->seqset);

  return;
}




/*
 * LTREE_recursive_print_tree() -
 *
 * Recursively (depth-first) traverse the tree nodes and print out details
 *
 */
static void 
LTREE_recursive_print_tree(LTREE_treenode *node,
			   ALIGN_seqset *seqset) {
  

  if(node->seq_id>=0) {

    if(seqset) {
      printf("%s", seqset->seq[node->seq_id].title);
    } else {
      printf("%d", node->seq_id);
    }

  } else {
    

    if(node->left && node->right) {
      printf("( ");
    }
    if(node->left) {
      LTREE_recursive_print_tree(node->left, seqset);
    }

    if(node->left && node->right) {
      printf(", ");
    }
    if(node->right) {
      LTREE_recursive_print_tree(node->right, seqset);
    }
    if(node->left && node->right) {
      printf(") ");
    }
  }

  return;
}






/*
 * LTREE_write_tree() - 
 *
 * Write a tree to disk.
 *
 */
void
LTREE_write_tree(LTREE_tree *tree,
		 FILE *fp) {

  LTREE_recursive_write_tree(tree->root,
			     tree->seqset,
			     fp);
  
  return;
}





/*
 * LTREE_recursive_write_tree() -
 *
 * Recursively write the tree to disk
 *
 */
static void 
LTREE_recursive_write_tree(LTREE_treenode *node,
			   ALIGN_seqset *seqset,
			   FILE *fp) {

  if(node->seq_id>=0) {
    fprintf(fp, "%s", seqset->seq[node->seq_id].title);
  } else {
    if(node->left && node->right) {
      fprintf(fp, "( ");
    }
    if(node->left) {
      LTREE_recursive_write_tree(node->left, seqset, fp);
    }

    if(node->left && node->right) {
      fprintf(fp, ", ");
    }
    if(node->right) {
      LTREE_recursive_write_tree(node->right, seqset, fp);
    }
    if(node->left && node->right) {
      fprintf(fp, ") ");
    }
  }

  return;
}




/*
 * LTREE_init_tree() -
 *
 * Allocate and construct a valid, random ltree and return it
 *
 */
LTREE_tree *
LTREE_init_tree(ALIGN_seqset *seqset) {

  unsigned int i, nodes, pos1, pos2, index;
  LTREE_treenode *node, *root;
  LTREE_treenode **buffer, **tmp_buffer;
  LTREE_tree *tree;

  
  nodes = seqset->num;
  buffer = (LTREE_treenode **)calloc(nodes, sizeof(LTREE_treenode *));

  /* we need nodes for all the terminals */
  for(i=0;i<nodes;i++) {

    buffer[i] = (LTREE_treenode *)calloc(1, sizeof(LTREE_treenode));
    buffer[i]->seq_id = i;
    buffer[i]->remove_flag = 0;
    buffer[i]->left   = NULL;
    buffer[i]->right  = NULL;
    buffer[i]->parent = NULL;
  }
  
  while(nodes>1) {

    /* determine the nodes to join */
    pos1 = MT_get_range_uint32(nodes);
    pos2 = MT_get_range_uint32(nodes);
    while(pos1 == pos2) {
      pos2 = MT_get_range_uint32(nodes);
    }

    node = (LTREE_treenode *)calloc(1, sizeof(LTREE_treenode));

    node->remove_flag   = 0;
    node->seq_id        = -1;
    node->left          = buffer[pos1];
    node->right         = buffer[pos2];
    node->parent        = NULL;
    node->left->parent  = node;
    node->right->parent = node;

    tmp_buffer = (LTREE_treenode **)calloc((nodes-1), sizeof(LTREE_treenode *));
    tmp_buffer[0] = node;

    index = 1;
    for(i=0;i<nodes;i++) {
      if(i!=pos1 && i!=pos2) {
	tmp_buffer[index++]=buffer[i];
      }
    }
    
    free(buffer);
    buffer = tmp_buffer;

    nodes--;
  }

  root = buffer[0];
  free(buffer);
  
  /* now lets allocate a tree */
  tree         = (LTREE_tree *)calloc(1, sizeof(LTREE_tree));
  tree->root   = root;
  tree->seqset = seqset;

  return(tree);
}




/*
 * LTREE_test_alignment() - 
 *
 * A simple function to test the speed, efficiency, 
 * and correctness of the alignment library
 *
 */
void
LTREE_test_alignment(void) {
  


  ALIGN_alignment *s1, *s2, *s3;
  ALIGN_alignment *t, *final;
  ALIGN_scoring_system ss;
  char aln_a[8] = {0, 0, 0, 1, 2, 0, 0, 0};
  char aln_b[2] = {1, 2};
  char aln_c[1] = {1};
  double fitness;

  ss.submat = 
  ALIGN_read_submat(DEFAULT_DNA_MATRIX);
  if(!ss.submat) {
    fprintf(stderr, "EVALYN: Fatal error parsing substitution matrix file\n");
    exit(-1);
  }
  
  ALIGN_print_submat(ss.submat);
  
  ss.gap_open   = -1.0;
  ss.gap_extend = -0.25;
  
  s1 = ALIGN_init_alignment(aln_a, ss.submat->symtab, 8, 0);
  s2 = ALIGN_init_alignment(aln_b, ss.submat->symtab, 2, 1);
  s3 = ALIGN_init_alignment(aln_c, ss.submat->symtab, 1, 2);
  
  ALIGN_print_alignment(s1, NULL);
  ALIGN_print_profile(s1);
  
  ALIGN_print_alignment(s2, NULL);
  ALIGN_print_profile(s2);

  ALIGN_print_alignment(s3, NULL);
  ALIGN_print_profile(s3);


  t = ALIGN_alignment_alignment(&ss,
				s1, 
				s2);

  
  ALIGN_print_alignment(t, NULL);
  ALIGN_print_profile(t);
  

  final = ALIGN_alignment_alignment(&ss,
				    t, 
				    s3);

  
  
  printf("\nFINAL ALIGNMENT\n");

  ALIGN_print_alignment(final, NULL);
  ALIGN_print_profile(final);

  exit(0);
  
  fitness = ALIGN_compute_sp(&ss,
			     t);
  
  //  ALIGN_print_profile(final);

  printf("FITNESS OF FINAL: %lf\n", fitness);
  
  /* free everything */
  ALIGN_free_alignment(s1);
  ALIGN_free_alignment(s2);
  
  ALIGN_free_alignment(final);

  return;
}










/*
 * LTREE_free_population() -
 *
 * Free an entire population here
 *
 */
void
LTREE_free_population(LTREE_population *population) {

  unsigned long int i;
  
  if(!population) {
    return;
  }

  if(population->trees) {
    for(i=0;i<population->size;i++) {
      LTREE_free_tree(population->trees[i]);
    }
 
    free(population->trees);
  }

  free(population);

  return;
}




/*
 * LTREE_compress_tree() - 
 *
 * Compress a tree
 *
 */
void
LTREE_compress_tree(LTREE_tree *tree) {

  LTREE_recursive_compress_tree(tree->root);

  return;
}



/*
 * LTREE_recursive_compress_tree() - 
 *
 * Recursive function to compress the tree here
 *
 */
static void
LTREE_recursive_compress_tree(LTREE_treenode *node) {
  
  LTREE_treenode *tmp;

  /* traverse left */
  if(node->left) {
    LTREE_recursive_compress_tree(node->left);
  }

  /* traverse right */
  if(node->right) {
    LTREE_recursive_compress_tree(node->right);
  }
  
  /* if internal node with one branch */
  if( ( (node->left)  && (!node->right) ) ||
      ( (!node->left) && (node->right)  ) ) {
    
    if(node->left) {
      tmp = node->left;
    } else {
      tmp = node->right;
    }
    
    if(node->alignment) {
      ALIGN_free_alignment(node->alignment);
    }

    node->seq_id      = tmp->seq_id;
    node->remove_flag = tmp->remove_flag;
    node->left        = tmp->left;
    node->right       = tmp->right;
    
    if(tmp->alignment) {
      ALIGN_free_alignment(tmp->alignment);
    }
    free(tmp);
    
    if(node->left) {
      node->left->parent = node;
    }
    
    if(node->right) {
      node->right->parent = node;
    }

  }

  return;
}



/*
 * LTREE_print_population()
 *
 * Print the whole population for debugging purposes 
 *
 */
static void
LTREE_print_population(LTREE_population *population) {
  
  unsigned short int i;

  printf("\n\n");
  printf("Population: %d\n", population->size);

  for(i=0;i<population->size;i++) {
    LTREE_print_tree(population->trees[i]);
    printf("\n");
  }
  
  printf("\n\n");
  
  return;
}





/*
 * LTREE_tree_depth() -
 *
 * Determine the depth of the specified tree here
 *
 */
unsigned int
LTREE_tree_depth(LTREE_treenode *node,
		 unsigned int depth) {

  unsigned int d_left, d_right;

  if(node->left) {
    d_left = LTREE_tree_depth(node->left, depth+1);
  } else {
    d_left = depth+1;
  }

  if(node->right) {
    d_right = LTREE_tree_depth(node->right, depth+1);
  } else {
    d_right = depth+1;
  }
  
  if(d_left > d_right) {
    return(d_left);
  } else {
    return(d_right);
  }
}





/*
 * LTREE_tree_size() - 
 *
 * return the size of the tree here
 *
 */
unsigned int
LTREE_tree_size(LTREE_treenode *node) {

  unsigned int s=0;

  s++;

  if(node->left) {
    s += LTREE_tree_size(node->left);
  }

  if(node->right) {
    s += LTREE_tree_size(node->right);
  }

  return(s);
}




/*
 * LTREE_get_vector() -
 *
 * Get all of the nodes at one level and put them in a vector
 *
 */
static void
LTREE_get_vector(LTREE_treenode *node,
		 unsigned int target_depth,
		 unsigned int current_depth,
		 LTREE_treenode **vector) {

  if(current_depth == target_depth) {
    
    vector[LTREE_vector_index++] = node;
    return;
  }
  
  if(node->left) {
    LTREE_get_vector(node->left,
		     target_depth,
		     current_depth+1,
		     vector);
  }

  if(node->right) {
    LTREE_get_vector(node->right,
		     target_depth,
		     current_depth+1, 
		     vector);
  }
}




/*
 * LTREE_nodes_at_depth() - 
 *
 * Determine the number of nodes at a given level in the tree
 *
 */
static unsigned int
LTREE_nodes_at_depth(LTREE_treenode *node,
		     unsigned int target_depth,
		     unsigned int current_depth) {


  unsigned int num_nodes = 0;
  
  if(node) {
    if(current_depth == target_depth) {
      return(1);
    }
  
    num_nodes += LTREE_nodes_at_depth(node->left,
				      target_depth,
				      current_depth+1);

    num_nodes += LTREE_nodes_at_depth(node->right,
				      target_depth,
				      current_depth+1);
    return(num_nodes);

  } else {
    return(0);
  }
}







/*
 * LTREE_build_taxa_list() - 
 *
 * Build a linked list of taxa (by seq_id).  This will be used 
 * during crossover
 *
 */ 
static void
LTREE_build_taxa_list(LTREE_treenode *node) {


  LTREE_taxa_list *h;

  /* if we are at a terminal, lets allocate a taxa node */
  if(!(node->left) && !(node->right)) {
    h = (LTREE_taxa_list *)calloc(1, sizeof(LTREE_taxa_list));
    h->seq_id = node->seq_id;
    h->next   = NULL;
    
    /* we allocated a node, now stick it on the end of the global linked list */
    LTREE_add_taxa_node(LTREE_global_taxa_list, h);
    
    return;
  }
  
  if(node->left) {
    LTREE_build_taxa_list(node->left);
  }
  
  if(node->right) {
    LTREE_build_taxa_list(node->right);
  }

  return;
}




/*
 * LTREE_add_taxa_node() - 
 *
 * Add a taxa node to a linked list of taxa nodes
 * 
 */
static void
LTREE_add_taxa_node(LTREE_taxa_list *head,
		    LTREE_taxa_list *node) {

  LTREE_taxa_list *n;
  
  if(head == NULL) {
    LTREE_global_taxa_list = node;
    return;
  }

  n = head;
  while(n->next) {
    n = n->next;
  }

  n->next = node;
  
  return;
}




/*
 * LTREE_print_taxa_list() -
 *
 * Well, it prints a taxa list
 * 
 */
static void
LTREE_print_taxa_list(LTREE_taxa_list *head) {
  
  LTREE_taxa_list *n;

  n = head;
  while(n) {
    printf("%d, ", n->seq_id);
    n = n -> next;
  }

  return;
}



/*
 * LTREE_count_taxa_list() - 
 * 
 */
long int
LTREE_count_taxa_list(LTREE_taxa_list *head) {
  
  LTREE_taxa_list *n;
  long int i = 0;

  n = head;
  while(n) {
    i++;
    n = n -> next;
  }

  return(i);
}







/*
 * LTREE_free_taxa_list() -
 *
 * Free a taxa list here
 *
 */
static void
LTREE_free_taxa_list(LTREE_taxa_list *head) {

  LTREE_taxa_list *n, *f;

  n = head;
  while(n) {
    if(n->next) {
      f = n->next;
      free(n);
      n = f;
    } else {
      free(n);
      n = NULL;
    }
  }
  
  return;
}






/*
 * LTREE_dup_tree() - 
 *
 * Allocate and duplicate a tree 
 *
 */
static LTREE_tree *
LTREE_dup_tree(LTREE_tree *tree) {
  
  LTREE_tree *ret_tree;
  
  ret_tree = (LTREE_tree *)calloc(1, sizeof(LTREE_tree));
  
  ret_tree->seqset     = tree->seqset;
  ret_tree->prng_state = tree->prng_state;
  ret_tree->fitness    = tree->fitness;
  
  ret_tree->root = LTREE_dup_treenodes(tree->root);
  
  return(ret_tree);
}





/*
 * LTREE_dup_treenodes() - 
 *
 * Recursively duplicate the node structure in a tree
 *
 */
static LTREE_treenode *
LTREE_dup_treenodes(LTREE_treenode *node) {
  
  LTREE_treenode *ret_node;
  
  if(node) {
    ret_node = (LTREE_treenode *)calloc(1, sizeof(LTREE_treenode));
    ret_node->seq_id       = node->seq_id;
    ret_node->remove_flag  = 0;
  } else {
    return(NULL);
  }
  
  if(node->left) {
    ret_node->left = LTREE_dup_treenodes(node->left);
    ret_node->left->parent = ret_node;
  } else {
    ret_node->left = NULL;
  }

  if(node->right) {
    ret_node->right = LTREE_dup_treenodes(node->right);
    ret_node->right->parent = ret_node;
  } else {
    ret_node->right = NULL;
  }
  
  return(ret_node);
}






/*
 * LTREE_tag_tree() -
 *
 * Mark all terminal nodes which have seq_ids which match 
 * entries in the taxa_list for later removal in stage 2
 *
 */
static void
LTREE_tag_tree(LTREE_tree *tree,
	       LTREE_taxa_list *taxa_list) {

  LTREE_recursive_tag_nodes(tree->root,
			    taxa_list);
  
  return;
}






/*
 * LTREE_recursive_tag_nodes() -
 *
 * Mark all terminal nodes which have seq_ids which match 
 * entries in the taxa_list for later removal in stage 2
 *
 */
static void
LTREE_recursive_tag_nodes(LTREE_treenode *node,
			  LTREE_taxa_list *taxa_list) {
  
  /* descend left */
  if(node->left) {
    LTREE_recursive_tag_nodes(node->left,
			      taxa_list);
  }

  /* descend right */
  if(node->right) {
    LTREE_recursive_tag_nodes(node->right,
			      taxa_list);
  }

  /* if we are at a leaf */
  if(node->seq_id >= 0) {
    
    if(LTREE_in_taxa_list(taxa_list,
			  node->seq_id)) {

      /* mark for removal here */
      node->remove_flag = 1;
    }
  }


  /* mark non-terminal nodes which will need to be removed */
  if(node->left && node->right) {
    
    if(node->left->remove_flag &&
       node->right->remove_flag) {

      node->remove_flag = 1;
    }
  }
  
  return;
}





/*
 * LTREE_remove_tagged_nodes() - 
 *
 * Remove all nodes which are marked for removal 
 *
 */
static void
LTREE_remove_tagged_nodes(LTREE_treenode *node) {
  
  if(node->left) {
    if(node->left->remove_flag) {
      LTREE_free_subtree(node->left);
      node->left = NULL;
    } else {
      LTREE_remove_tagged_nodes(node->left);
    }
  }

  if(node->right) {
    if(node->right->remove_flag) {
      LTREE_free_subtree(node->right);
      node->right = NULL;
    } else {
      LTREE_remove_tagged_nodes(node->right);
    }
  }
  
  return;
}



/*
 * 
 * LTREE_in_taxa_list() - 
 *
 * Return 1 if the sequence is contained within the taxa list or 0 otherwise
 *
 */
static unsigned int
LTREE_in_taxa_list(LTREE_taxa_list *taxa_list,
		   short int seq_id) {

  LTREE_taxa_list *n;
  
  n = taxa_list;
  while(n) {
    if(n->seq_id == seq_id) {
      return(1);
    }
    n = n->next;
  }
  
  return(0);
}





/*
 * LTREE_free_tree() - 
 *
 * Free an entire tree, starting with a recursive
 * freeing of the treenodes
 *
 */
void
LTREE_free_tree(LTREE_tree *tree) {


  LTREE_free_subtree(tree->root);
  free(tree);
  
  return;
}




/*
 * LTREE_free_subtree() - 
 *
 * Free a tree recursively
 * 
 */
static void
LTREE_free_subtree(LTREE_treenode *node) {
  
  if(!node) {
    return;
  }

  if(node->left) {
    LTREE_free_subtree(node->left);
  }
  
  if(node->right) {
    LTREE_free_subtree(node->right);
  }

  /* free the alignment at this node, if there is one */
  if(node->alignment) {
    ALIGN_free_alignment(node->alignment);
  }
  
  free(node);

  return;
}






/*
 * LTREE_clear_internal_alignments() - 
 *
 * Clear all internal alignments here
 *
 */
static void
LTREE_clear_internal_alignments(LTREE_tree *tree) {
  
  LTREE_recursive_clear_internal_alignments(tree->root);
  
  return;
}




/*
 * LTREE_recursive_clear_internal_alignents() -
 *
 * Free all of the alignments which are not on leaf nodes
 *
 */
static void
LTREE_recursive_clear_internal_alignments(LTREE_treenode *node) {
  
  if(node->left) {
    LTREE_recursive_clear_internal_alignments(node->left);
  }
    
  if(node->right) {
    LTREE_recursive_clear_internal_alignments(node->right);
  }

  if(node->alignment) {
    ALIGN_free_alignment(node->alignment);
    node->alignment = NULL;
  }
  
  return;
}




/*
 * LTREE_init_leaf_alignments() - 
 *
 */
static void
LTREE_init_leaf_alignments(LTREE_tree *tree) {

  LTREE_recursive_init_leaf_alignments(tree->root,
				       tree->seqset);
  
  return;
}



/*
 * 
 * LTREE_recursive_init_leaf_alignments() - 
 *
 * Create the initial sequence alignments for terminal nodes in a given tree
 *
 */
static void
LTREE_recursive_init_leaf_alignments(LTREE_treenode *node,
				     ALIGN_seqset *seqset) {

  /* traverse left */
  if(node->left) {
    LTREE_recursive_init_leaf_alignments(node->left, seqset);
  }
  
  /* traverse right */
  if(node->right) {
    LTREE_recursive_init_leaf_alignments(node->right, seqset);
  }

  /* create initial alignment */
  if(node->seq_id >= 0) {
    node->alignment = 
      ALIGN_init_alignment(seqset->seq[node->seq_id].data, 
			   seqset->symtab,
			   seqset->seq[node->seq_id].length,
			   node->seq_id);
  }

  return;
}




/*
 *
 * LTREE_handle_options() - 
 *
 * Here, we handle all of the command-line options sent
 * to the program
 *
 */
static void
LTREE_handle_options(int argc,
		     char *argv[],
		     LTREE_cl_args *cl_args,
		     ALIGN_scoring_system *ss) {

  int option_index, c;
  int alg_count = 0;

  /* define our valid command-line options */
  struct option LTREE_long_options[] = {

    /* These options set a flag */
    {"verbose",    no_argument, &(cl_args->verbose_flag),        1},
    {"quiet",      no_argument, &(cl_args->quiet_flag),          1},
    {"debug",      no_argument, &(cl_args->debug_flag),          1},
    {"dna",        no_argument, &(cl_args->dna_flag),            1},
    {"test",       no_argument, &(cl_args->test_flag),           1},
    {"protein",    no_argument, &(cl_args->protein_flag),        1},
    {"parsimony",  no_argument, &(cl_args->parsimony_flag),      1},
    {"sps",        no_argument, &(cl_args->sps_flag),            1},
    {"rnj",        no_argument, &(cl_args->rnj_flag),            1},
    {"ga",         no_argument, &(cl_args->search_algorithm),    LTREE_ALGORITHM_GA},
    {"randomwalk", no_argument, &(cl_args->search_algorithm),    LTREE_ALGORITHM_RANDOM_WALK},
    {"hillclimb",  no_argument, &(cl_args->search_algorithm),    LTREE_ALGORITHM_HILL_CLIMB},
    {"anneal",     no_argument, &(cl_args->search_algorithm),    LTREE_ALGORITHM_SIMULATED_ANNEALING},

    /* these options don't set a flag */
    {"infile",      required_argument, NULL, 'f'},
    {"outfile",     required_argument, NULL, 'o'},
    {"treefile",    required_argument, NULL, 'T'},
    {"infiletype",  required_argument, NULL, 'I'},
    {"outfiletype", required_argument, NULL, 'O'},
    {"matrix",      required_argument, NULL, 'm'},
    {"iterations",  required_argument, NULL, 'i'},
    {"population",  required_argument, NULL, 'p'},
    {"mutation",    required_argument, NULL, 'M'},
    {"gapopen",     required_argument, NULL, 'g'},
    {"gapextend",   required_argument, NULL, 'G'},
    {"seed",        required_argument, NULL, 's'},
    {"converge",    required_argument, NULL, 'c'},
    {0, 0, 0, 0}

  };
  
  /* initialize CL args here */
  cl_args->verbose_flag             = 0;
  cl_args->quiet_flag               = 0;
  cl_args->debug_flag               = 0;
  cl_args->protein_flag             = 0;
  cl_args->dna_flag                 = 0;
  cl_args->seq_type                 = 0;
  cl_args->test_flag                = 0;
  cl_args->parsimony_flag           = 0;
  cl_args->sps_flag                 = 0;
  cl_args->rnj_flag                 = 0;
  cl_args->search_algorithm         = LTREE_DEFAULT_ALGORITHM;
  cl_args->converge_flag            = 0;
  cl_args->stationarity             = DEFAULT_CONVERGE_STATIONARITY;
  cl_args->seed                     = 0;
  cl_args->iterations               = DEFAULT_ITERATIONS;
  cl_args->mutation_rate            = DEFAULT_MUTATION_RATE;
  cl_args->population_size          = DEFAULT_POPULATION_SIZE;
  cl_args->gap_open                 = 0;
  cl_args->gap_extend               = 0;
  cl_args->gap_open_set             = 0;
  cl_args->gap_extend_set           = 0;
  cl_args->input_filename[0]        = '\0';
  cl_args->output_filename[0]       = '\0';
  cl_args->tree_filename[0]         = '\0';
  cl_args->matrix_filename[0]       = '\0';
  cl_args->infiletype               = DEFAULT_INFILETYPE;
  cl_args->outfiletype              = DEFAULT_OUTFILETYPE;

  /* read and process arguments here */
  while(1) {
    c = getopt_long(argc,
		    argv,
		    "f:o:T:I:O:m:i:M:p:g:G:s:c:dqvtPSrelRha",
		    LTREE_long_options,
		    &option_index);
    
    if(c == -1) {
      break;
    }
 
    switch(c) {
    case 0:
      if(LTREE_long_options[option_index].flag) {
	break;
      }
      
      printf("option %s", LTREE_long_options[option_index].name);
      if(optarg) {
	printf(" with arg %s", optarg);
      }
      printf("\n");
      break;
      
    case 'i':

      cl_args->iterations = atoi(optarg);
      break;
      
    case 'M':
      
      cl_args->mutation_rate = atof(optarg);
      break;

    case 'f':

      strncpy(cl_args->input_filename, optarg, MAXPATHLEN);
      break;

    case 'o':

      strncpy(cl_args->output_filename, optarg, MAXPATHLEN);
      break;

    case 'T':

      strncpy(cl_args->tree_filename, optarg, MAXPATHLEN);
      break;
      
    case 'p':
      
      cl_args->population_size = atoi(optarg);
      break;

    case 's':
      
      cl_args->seed = atoi(optarg);
      break;

    case 'm':

      strncpy(cl_args->matrix_filename, optarg, MAXPATHLEN);
      break;

    case 'g':

      cl_args->gap_open = atof(optarg);
      cl_args->gap_open_set = 1;

      break;

    case 'G':

      cl_args->gap_extend = atof(optarg);
      cl_args->gap_extend_set = 1;
      break;
      
    case 'P':

      cl_args->parsimony_flag = 1;
      
      break;
      
    case 'S':

      cl_args->sps_flag = 1;

      break;

    case 'r':

      cl_args->rnj_flag = 1;
      
      break;

    case 'e':
      
      cl_args->search_algorithm = LTREE_ALGORITHM_GA;
      alg_count++;
      
      break;
      
    case 'R':
      
      cl_args->search_algorithm = LTREE_ALGORITHM_RANDOM_WALK;
      alg_count++;
      
      break;
      
    case 'h':
      
      cl_args->search_algorithm = LTREE_ALGORITHM_HILL_CLIMB;
      alg_count++;
      
      break;

    case 'a':
      
      cl_args->search_algorithm = LTREE_ALGORITHM_SIMULATED_ANNEALING;
      alg_count++;
      
      break;
      
    case 'c':

      cl_args->converge_flag = 1;
      cl_args->stationarity = atoi(optarg);
      if(cl_args->stationarity <= 0) {
	cl_args->stationarity = DEFAULT_CONVERGE_STATIONARITY;
      }
      
      break;

    case 'q':
      cl_args->verbose_flag = 0;
      cl_args->quiet_flag   = 1;
      
      break;

    case 'v':
      cl_args->verbose_flag = 1;
      cl_args->quiet_flag   = 0;
      
      break;

    case 'I':
      if(!strcmp(optarg, "MSF") ||
	 !strcmp(optarg, "msf") ||
	 !strcmp(optarg, "Msf")) {

	cl_args->infiletype = LTREE_MSF;
      } else {
	cl_args->infiletype = LTREE_FASTA;
      }
      
      break;

    case 'O':
      if(!strcmp(optarg, "MSF") ||
	 !strcmp(optarg, "msf") ||
	 !strcmp(optarg, "Msf")) {

	cl_args->outfiletype = LTREE_MSF;

      } else if(!strcmp(optarg, "ALN") ||
		!strcmp(optarg, "aln") ||
		!strcmp(optarg, "Aln")) {

	cl_args->outfiletype = LTREE_ALN;

      } else {
	fprintf(stderr, "EVALYN: Unknown specified output file type.\n");
	exit(-1);
      }
       
      break;

    case 'd':
      cl_args->debug_flag = 1;
      
      break;
      
    case 't':
      cl_args->test_flag = 1;

      break;
      
    case '?':
      exit(-1);
      
    default:
      /*      abort(); */
      break;
    }
  }
  
  /* see if there remains any unprocessed portion of the command-line */
  if(optind < argc) {
    fprintf(stderr, "EVALYN: Unknown command-line argument:\n --> %s\n", argv[optind]);
    exit(-1);
  }

  /* if a test has not been specified */
  if(!cl_args->test_flag) {

    /* make sure all mandatory command-line args are provided */
    if(cl_args->input_filename[0] == '\0') {
      fprintf(stderr, "EVALYN: No input file specified\n");
      exit(-1);
    }

    /* if there is no explicit output filename, make one */
    LTREE_auto_output_filename(cl_args);

    /* make sure that user does not specify both DNA and Protein input */
    if(cl_args->protein_flag && cl_args->dna_flag) {
      fprintf(stderr, "EVALYN: Protein and DNA are mutually exclusive\n");
      exit(-1);
    } 
    
    /* make sure user doesn't specify both verbose and quiet */
    if(cl_args->verbose_flag && cl_args->quiet_flag) {
      fprintf(stderr, "EVALYN: Verbose mode and quiet mode are mutually exclusive\n");
      exit(-1);
    }

    /* the selection routine requires a population */
    if(cl_args->population_size < LTREE_MIN_POPSIZE) {
      fprintf(stderr, "EVALYN: Minimum population size is %d\n", LTREE_MIN_POPSIZE);
      exit(-1);
    }
    
    /* make sure that ambiguous, mutually exclusive search algorithms have not been specified */
    if(alg_count == 0) {
      cl_args->search_algorithm = LTREE_DEFAULT_ALGORITHM;
    } else if(alg_count > 1) {
      fprintf(stderr, "EVALYN: You must specify only one search algorithm (--ga, --random_walk, --hill_climb, --anneal)\n");
      exit(-1);
    }

  } else {  /* TEST */
    cl_args->dna_flag     = 1;
    cl_args->protein_flag = 0;
    cl_args->verbose_flag = 1;
    cl_args->debug_flag   = 1;

    LTREE_test_alignment();
    exit(0);
  }
  

  /* handle the command-line setting of the scoring method */
  if(cl_args->parsimony_flag && cl_args->sps_flag) {
    fprintf(stderr, 
	    "EVALYN: Parsimony and SPS are mutually exclusive scoring mechanisms.\n");
    exit(-1);
  } else {

    if( (!cl_args->parsimony_flag && !cl_args->sps_flag) ||
	(!cl_args->parsimony_flag && cl_args->sps_flag) ) {
      cl_args->scoring_method = LTREE_SPS;
      cl_args->sps_flag       = 1;
      cl_args->parsimony_flag = 0;
    } else {
      cl_args->sps_flag       = 0;
      cl_args->parsimony_flag = 1;
      cl_args->scoring_method = LTREE_PARSIMONY;
    }
  }
  
  if(cl_args->protein_flag) {
    cl_args->seq_type = ALIGN_PROTEIN;
  }
  
  if(cl_args->dna_flag) {
    cl_args->seq_type = ALIGN_DNA;
  }

  /* set the random seed here */
  LTREE_set_seed(cl_args);

  LTREE_autodetect_input_format(cl_args);  /* FASTA or MSF */

  LTREE_determine_input_type(cl_args);     /* PROTEIN OR DNA */
  
  LTREE_assign_gap_penalties(ss, cl_args);  
  
  if(cl_args->matrix_filename[0] == '\0') {
    if(cl_args->seq_type == ALIGN_PROTEIN) {
      strcpy(cl_args->matrix_filename, DEFAULT_PRO_MATRIX);
    } else {
      strcpy(cl_args->matrix_filename, DEFAULT_DNA_MATRIX);
    }
  }
  
  return;
}




/*
 * LTREE_print_options() -
 *
 * Pretty-print the command-line arguments
 *
 */
static void
LTREE_print_options(LTREE_cl_args *cl_args) {
  
  printf("-->> EVALYN COMMAND-LINE ARGUMENTS <<--\n");

  printf("INPUT FILE: %s\n",  cl_args->input_filename);
  printf("OUTPUT FILE: %s\n", cl_args->output_filename);
  printf("TREE FILE: %s\n",   cl_args->tree_filename);

  if(cl_args->infiletype == LTREE_MSF) {
    printf("INPUT FILETYPE: MSF\n");
  } else if(cl_args->infiletype == LTREE_FASTA) {
    printf("INPUT FILETYPE: FASTA\n");
  } else {
    printf("INPUT FILETYPE: <UNKNOWN>\n");
  }


  if(cl_args->outfiletype == LTREE_MSF) {
    printf("OUTPUT FILETYPE: MSF\n");
  } else if(cl_args->outfiletype == LTREE_ALN) {
    printf("OUTPUT FILETYPE: ALN\n");
  } else {
    printf("OUTPUT FILETYPE: <UNKNOWN>\n");
  }

  
  if(cl_args->matrix_filename) {
    printf("SUBSTITUTION MATRIX: %s\n", cl_args->matrix_filename);
  }
  
  printf("VERBOSE: [%c]\n",    (cl_args->verbose_flag?'Y':'N'));
  printf("QUIET: [%c]\n",      (cl_args->quiet_flag?'Y':'N'));
  printf("DEBUG: [%c]\n",      (cl_args->debug_flag?'Y':'N'));
  printf("SEQ TYPE: [%s]\n",   (cl_args->seq_type==ALIGN_PROTEIN?"PROTEIN":"DNA"));
  printf("TEST: [%c]\n",       (cl_args->test_flag?'Y':'N'));

  printf("SCORING: [%s]\n",    (cl_args->scoring_method==LTREE_PARSIMONY?"PARSIMONY":"SPS"));
  
  if(cl_args->rnj_flag) {
    printf("SEEDING INITIAL POPULATION WITH A RELAXED NEIGHBOR-JOINING TREE\n");
  }
  
  if(cl_args->converge_flag) {
    printf("CONVERGENCE: [STATIONARITY = %d]\n", cl_args->stationarity);
  } else {
    printf("CONVERGENCE: [FIXED ITERATIONS = %d]\n", cl_args->iterations);
  }
  
  switch(cl_args->search_algorithm) {
  case LTREE_ALGORITHM_GA:
    printf("SEARCH ALGORITHM: GENETIC ALGORITHM WITH LEWIS CROSSOVER\n");
    break;
  case LTREE_ALGORITHM_RANDOM_WALK:
    printf("SEARCH ALGORITHM: RANDOM WALK\n");
    break;
  case LTREE_ALGORITHM_HILL_CLIMB:
    printf("SEARCH ALGORITHM: SPR HILL CLIMB\n");
    break;
  case LTREE_ALGORITHM_SIMULATED_ANNEALING:
    printf("SEARCH ALGORITHM: SIMULATED ANNEALING\n");
    break;
  default:
    fprintf(stderr, "*** NO SEARCH ALGORITHM ***\n");
 
  }

  printf("PRNG SEED: [%ld]\n", cl_args->seed);

  printf("ITERATIONS: [%d]\n",      cl_args->iterations);
  printf("MUTATION RATE: [%lf]\n",   cl_args->mutation_rate);
  printf("POPULATION SIZE: [%d]\n", cl_args->population_size);
  
  printf("GAP OPEN COST: [%lf]\n",      cl_args->gap_open);
  printf("GAP EXTENSION COST: [%lf]\n", cl_args->gap_extend);

  
  return;
}







/*
 * LTREE_set_seed()
 *
 * Set the PRNG seed here
 *
 */
static void
LTREE_set_seed(LTREE_cl_args *cl_args) {

  int pid;

  /* set random seed here */
  if(!cl_args->seed) {
    pid = getpid();
    cl_args->seed = (time(0)^(pid<<8));
  }

  /* set the Mersenne Twister seed */
  MT_init_genrand(cl_args->seed);
  
  return;
}




/*
 * LTREE_auto_output_filename() - 
 *
 * If there is no explicit output filenames specified, then construct 
 * defaults based on the infile and other command-line parameters 
 *
 * We do this for trees and alignment output
 *
 */
static void
LTREE_auto_output_filename(LTREE_cl_args *cl_args) {

  int i, len;
  char new_extension[5];


  /* OUTPUT FILENAME */
  if(cl_args->output_filename[0] == '\0') {
  
    switch(cl_args->outfiletype) {
    case LTREE_MSF:
    
      strcpy(new_extension, ".msf");

      break;

    case LTREE_ALN:
    
      strcpy(new_extension, ".aln");
    
      break;

    default:
      fprintf(stderr, "Output file format unsupported\n");
      exit(-1);
    }
  
    strcpy(cl_args->output_filename, cl_args->input_filename);
    
    /* shorten the string at the last dot */
    len = strlen(cl_args->output_filename);
    for(i=len-1;i>=0;i--) {
      if(cl_args->output_filename[i] == '.') {
	cl_args->output_filename[i] = '\0';
	break;
      }
    }

    /* append the new extension to the base of the filename */
    strcat(cl_args->output_filename, new_extension);
    
    if(!strcmp(cl_args->input_filename, cl_args->output_filename)) {
      sprintf(cl_args->output_filename, "out_%s", cl_args->input_filename);
    }
  }


  /* TREE FILENAME */
  if(cl_args->tree_filename[0] == '\0') {
      strcpy(new_extension, ".dnd");
  
      strcpy(cl_args->tree_filename, cl_args->input_filename);

      /* shorten the string at the last dot */
      len = strlen(cl_args->tree_filename);
      for(i=len-1;i>=0;i--) {
	if(cl_args->tree_filename[i] == '.') {
	  cl_args->tree_filename[i] = '\0';
	  break;
	}
      }

      /* append the new extension to the base of the filename */
      strcat(cl_args->tree_filename, new_extension);
  }

  return;
}






/*
 * LTREE_determine_input_type() -
 *
 * A function which determines the input type
 * (protein or DNA)
 *
 * RETURNS:
 * --------
 *
 *  ALIGN_DNA
 *  ALIGN_PROTEIN
 *  -1 if things are ambiguous
 *
 */
static int
LTREE_determine_input_type(LTREE_cl_args *cl_args) {

  
  FILE *fp;
  unsigned int i, j, k;
  char line[MAX_SEQ_LENGTH];
  char uniq[MAX_SYMS] = "", uniq_t[MAX_SYMS], tmp[2], *c=NULL;
  int found, flag;


  switch(cl_args->infiletype) {
  case LTREE_MSF:

    /* open specified fasta file here */
    fp = fopen(cl_args->input_filename, "r");
    if(!fp) {
      fprintf(stderr, "EVALYN:  Failed to open input MSF file: %s\n", cl_args->input_filename);
      exit(-1);
    }

    /* read file and grab sequence names until you get to the demarcation (double slash) */
    flag = 0;
    while(!feof(fp)) {

      fgets(line, MAX_SEQ_LENGTH, fp);
      
      if(line[0] == '/' && line[1] == '/') {
	fgets(line, MAX_SEQ_LENGTH, fp);
        flag = 1;
      }

      if(flag) {
      
	for(i=0;i<strlen(line);i++) {
	  if(line[i] == ' ' ||
	     line[i] == '\t') {
	    c = &(line[i]);
	    break;
	  }
	}

	if(c) {
	  
	  for(i=0;i<strlen(c);i++) {
	
	    found = 0;
	    for(k=0;k<strlen(uniq) && !found;k++) {
	      if(uniq[k] == toupper(c[i])) {
		found = 1;
		break;
	      }
	    }

	    if(!found) {
	      sprintf(tmp, "%c", toupper(c[i]));
	      strcat(uniq, tmp);
	    }
	  }
	}
      }
    }

    /* clean up uniq */
    k = 0;
    for(i=0;i<strlen(uniq);i++) {
      if(uniq[i] >= 'A' && uniq[i] <= 'Z') {
	uniq_t[k++] = uniq[i];
      }
    }
    uniq_t[k] = '\0';
    
    strcpy(uniq, uniq_t);

    printf("strlen(uniq) --> %d\n", strlen(uniq));
    printf("UNIQ: %s\n", uniq);
 
    if(strlen(uniq) > ALIGN_NUM_DNA_SYMS) {
      fclose(fp); 
      if(!cl_args->seq_type) {
	cl_args->seq_type = ALIGN_PROTEIN;
      }
      return(ALIGN_PROTEIN);

    } else {
      fclose(fp); 
      if(!cl_args->seq_type) {
	cl_args->seq_type = ALIGN_DNA;
      }
      return(ALIGN_DNA);

    }

    break;

  case LTREE_FASTA:

    /* open specified fasta file here */
    fp = fopen(cl_args->input_filename, "r");
    if(!fp) {
      fprintf(stderr, "EVALYN:  Failed to open input FASTA file: %s\n", cl_args->input_filename);
      exit(-1);
    }

    /* read every line of the file and parse */
    i = 0;
    while (!feof(fp)) {                      
    
      fgets(line, sizeof(line), fp);
      if(feof(fp)) {
	break;
      }

      line[strlen(line) - 1] = '\0';  /* remove the CR */

      /* parse the sequence titles */
      if(line[0] != '>') {

	if(strlen(line)>2) {

	  for(j=0;j<strlen(line);j++) {

	    /* try to add every character to the string of uniq symbols */
	    /* if it is not there, add it                               */
	    /* we do this to count the number of unique symbols         */
	    found = 0;
	    for(k=0;k<strlen(uniq) && !found;k++) {
	      if(uniq[k] == toupper(line[j])) {
		found = 1;
		break;
	      }
	    }

	    sprintf(tmp, "%c", toupper(line[j]));
	    if(!found) {
	      strcat(uniq, tmp);
	    }
	  }
	  i++;
	}
      }
    }

    /* determine type of sequences based on types of symbols */
    if(strlen(uniq) > ALIGN_NUM_DNA_SYMS) {
      fclose(fp); 
      if(!cl_args->seq_type) {
	cl_args->seq_type = ALIGN_PROTEIN;
      }
      return(ALIGN_PROTEIN);

    } else {
      fclose(fp); 
      if(!cl_args->seq_type) {
	cl_args->seq_type = ALIGN_DNA;
      }
      return(ALIGN_DNA);

    }

    break;

  default:
    fprintf(stderr, "EVALYN: Unknown infiletype in LTREE_determine_input_type()\n");
    exit(-1);
  }
}






/*
 * 
 * LTREE_autodetect_input_format()
 *
 * Detect the input format (MSF or FASTA)
 *
 */
static int
LTREE_autodetect_input_format(LTREE_cl_args *cl_args) {
  
  FILE *fp;
  int i;
  char extension[16] = "";
  char line[MAX_SEQ_LENGTH];

  
  /* lets look at file extensions for a clue */
  for(i=strlen(cl_args->input_filename)-1;i>0;i--) {
    if(cl_args->input_filename[i] == '.') {
      strcpy(extension, &(cl_args->input_filename[i+1]));
    }
  }

  for(i=0;i<strlen(extension);i++) {
    extension[i] = toupper(extension[i]);
  }
  
  if(!strcmp(extension, "MSF")) {
    cl_args->infiletype = LTREE_MSF;
    return(LTREE_MSF);
  }

  /* extension told us nothing, lets look deeper */
  fp = fopen(cl_args->input_filename, "r");
  if(!fp) {
    fprintf(stderr, "EVALYN: Could not open input file %s\n", cl_args->input_filename);
    exit(-1);
  }
  
  /* MSF files have this wacky separator line with // in the middle of the file
   * so clue off of that
   */
  while(!feof(fp)) {
    fgets(line, MAX_SEQ_LENGTH, fp);

    if(line[0] == '/' && line[1] == '/') {
      /* this indicates an MSF file */
      fclose(fp);
      cl_args->infiletype = LTREE_MSF;
      return(LTREE_MSF);
    }
  }
  
  fclose(fp);
  cl_args->infiletype = LTREE_FASTA;

  /* Its not MSF, so its probably FASTA */
  return(LTREE_FASTA);
}






/*
 * LTREE_assign_gap_penalties() -
 *
 * Assign the proper gap penalties here
 *
 */
static void
LTREE_assign_gap_penalties(ALIGN_scoring_system *ss,
			   LTREE_cl_args *cl_args) {
  
  /* handle the gap open penalties */
  if(cl_args->gap_open_set) {
    ss->gap_open = cl_args->gap_open;
  } else {
    if(cl_args->seq_type == ALIGN_PROTEIN) {
      ss->gap_open = cl_args->gap_open = DEFAULT_PRO_GAPOPEN;
    } else if(cl_args->seq_type == ALIGN_DNA) {
      ss->gap_open = cl_args->gap_open = DEFAULT_DNA_GAPOPEN;
    }
  }

  /* handle the gap extend penalties */
  if(cl_args->gap_extend_set) {
    ss->gap_extend = cl_args->gap_extend;
  } else {
    if(cl_args->seq_type == ALIGN_PROTEIN) {
      ss->gap_extend = cl_args->gap_extend = DEFAULT_PRO_GAPEXTEND;
    } else if(cl_args->seq_type == ALIGN_DNA) {
      ss->gap_extend = cl_args->gap_extend = DEFAULT_DNA_GAPEXTEND;
    }
  }

  return;
}




/*
 * LTREE_write_alignment()
 *
 * write the alignment to disk
 *
 */
void
LTREE_write_alignment(ALIGN_alignment *alignment,
		      ALIGN_seqset *seqset,
		      LTREE_cl_args *cl_args) {

  switch(cl_args->outfiletype) {
  case LTREE_ALN:

    /* write CLUSTAL-like alignment output */
    ALIGN_write_alignment_ALN(alignment, 
			      seqset, 
			      cl_args->output_filename);

    break;

  case LTREE_MSF:
    
    /* write GCG wisconsin-package output here */
    ALIGN_write_alignment_MSF(alignment, 
			      seqset, 
			      cl_args->output_filename,
			      cl_args->seq_type);
    break;
    
  default:
    fprintf(stderr, "EVALYN: Unknown output file type\n");
    exit(0);
  }

  return;
}






/*
 * LTREE_dump_population()
 *
 * Dumps an entire population of trees.  Used for debugging and analysis.
 *
 */
void
LTREE_dump_population(char *filename,
		      LTREE_population *population,
		      int iteration,
		      LTREE_tree *parent1,
		      LTREE_tree *parent2,
		      LTREE_tree *child) {

  FILE *fp;
  int i;
  double meanfit, accumulated_fitness;


  /* open trace file here */
  fp = fopen(filename, "a");
  if(!fp) {
    fprintf(stderr, "EVALYN: Could not open population trace file %s\n", filename);
    return;
  }
  
  /* print some header info here */
  if(iteration == 0) {
    fprintf(fp, "\n\n**********************************************\n");
    fprintf(fp, "**********************************************\n");
    fprintf(fp, "INITIAL POPULATION\n");
    fprintf(fp, "**********************************************\n");
    fprintf(fp, "**********************************************\n\n");


  } else if(iteration < 0) {

    fprintf(fp, "\n\n**********************************************\n");
    fprintf(fp, "**********************************************\n");
    fprintf(fp, "FINAL POPULATION\n");
    fprintf(fp, "**********************************************\n");
    fprintf(fp, "**********************************************\n\n");

  } else {

    fprintf(fp, "\n\n**********************************************\n");
    fprintf(fp, "**********************************************\n");
    fprintf(fp, "ITERATION: %d\n", iteration);
    fprintf(fp, "**********************************************\n");
    fprintf(fp, "**********************************************\n\n");

    fprintf(fp, "  POPULATION SIZE: %d\n",        population->size);
    fprintf(fp, "  FITNESS OF BEST TREE:  %lf\n", population->trees[0]->fitness);
    fprintf(fp, "  FITNESS OF WORST TREE: %lf\n", population->trees[population->size-1]->fitness);
    
    accumulated_fitness = 0.0;
    for(i=0;i<population->size;i++) {
      accumulated_fitness += population->trees[i]->fitness;
    }

    meanfit = accumulated_fitness/(double)(population->size);
    
    fprintf(fp, "  MEAN TREE FITNESS: %f\n", meanfit);

    if(parent1 && parent2 && child) {
      fprintf(fp, "PARENT 1 (fitness = %lf):\n", parent1->fitness);
      LTREE_write_tree(parent1, fp);
      fprintf(fp, " ;\n");
      fprintf(fp, "PARENT 2 (fitness = %lf):\n", parent2->fitness);
      LTREE_write_tree(parent2, fp);
      fprintf(fp, " ;\n");
      fprintf(fp, "CHILD (fitness = %lf):\n", child->fitness);
      LTREE_write_tree(child, fp);
      fprintf(fp, " ;\n");
    }
  }

  
  /* loop through the population and print tree and associated fitness */
  for(i=0;i<population->size;i++) {
    fprintf(fp, "*************\n");
    fprintf(fp, "TREE #%d\n", i);
    fprintf(fp, "TREE FITNESS: %lf\n", population->trees[i]->fitness);
    LTREE_write_tree(population->trees[i], fp);
    fprintf(fp, ";\n");
  }
  
  fclose(fp);

  return;
}






/*
 * LTREE_lxo_step() -  
 * 
 * Perform a lewis crossover step by duplicating a tree and performing
 * a crossover at every possible crossover point and every possible 
 * insertion point.  Constructs a population of children.
 *
 * This should not require any direct recursion except to build
 * exhaustive lists of crossover points and insertin points (to be
 * handled in a different function).
 *
 */
LTREE_population *
LTREE_lxo_step(LTREE_tree *source_tree) {

  LTREE_tree *parent1=NULL, *parent2=NULL;
  LTREE_tree *child = NULL;
  LTREE_population *neighborhood;
  LTREE_taxa_list *tlist_node;

  //  LTREE_treenode **inodes;  // all internal nodes
  LTREE_treenode **cps;     // crossover points 
  
  ALIGN_seqset *seqset;

  LTREE_taxa_list *tlist_cp;

  long int i, j;
  long int ntaxa, n_inodes, hoodsize;
  long int tlist_size;
  LTREE_tree_desc *parent2_desc;

  
  // determine the number of taxa in the source tree
  ntaxa    = source_tree->seqset->num;
  n_inodes   = 2*ntaxa-2;
  hoodsize = n_inodes*2;  // for now
  
  // assign seqset alias
  seqset = source_tree->seqset;
  
  // allocate the neighborhood 
  neighborhood        = (LTREE_population *)calloc(1, sizeof(LTREE_population));
  neighborhood->size  = hoodsize;
  neighborhood->trees = (LTREE_tree **)calloc(hoodsize, sizeof(LTREE_tree *));

  // create a copy of the source tree to be parent1
  parent1 = LTREE_dup_tree(source_tree);

  // recursively traverse parent1 and return all of the internal crossover points
  cps = LTREE_enumerate_cps(parent1, n_inodes);

  // iterate here over all crossover points in parent1
  for(i=0;i<n_inodes;i++) {
    
    // make backup of source tree to be parent2
    parent2 = LTREE_dup_tree(source_tree);  

    // make a list of all taxa under this crossover point
    LTREE_global_taxa_list = NULL;
    LTREE_build_taxa_list(cps[i]);
    tlist_cp   = LTREE_global_taxa_list;
    tlist_size = LTREE_count_taxa_list(tlist_cp);

#ifdef FOO
    // for debugging, print all of the leaf titles here
    tlist_node = tlist_cp;
    while(tlist_node != NULL) {
      printf("%s\n", seqset->seq[tlist_node->seq_id].title);
      tlist_node = tlist_node->next;
    }
#endif /* FOO */

    // remove all nodes from parent2 which occur under cp in parent1
    LTREE_recursive_tag_nodes(parent2->root, tlist_cp);
    LTREE_free_taxa_list(tlist_cp);
    LTREE_global_taxa_list = NULL;
    LTREE_remove_tagged_nodes(parent2->root);
    LTREE_compress_tree(parent2);
    
    // build a pointer-independent representation of the compressed parent2
    parent2_desc = LTREE_build_tree_desc(parent2);
    //    LTREE_print_tree_desc(parent2_desc);
    
    // for every crossover point, crossover
    for(j=0;j<parent2_desc->num;j++) {
      child = LTREE_nondestructive_join(cps[i], parent2, parent2_desc->desc[j]);
    }
    
    // free the tree description for parent2
    LTREE_free_tree_desc(parent2_desc);

  }

  return(neighborhood);
}





/*
 * LTREE_enumerate_cps() - 
 *
 * Take the root of a tree and traverse the tree, building an
 * exhaustive list of valid crossover points in the tree.
 *
 */
LTREE_treenode **
LTREE_enumerate_cps(LTREE_tree *tree,
		    long int ncps) {

  LTREE_treenode **cps;

  cps = (LTREE_treenode **)calloc(ncps, sizeof(LTREE_treenode *));
  LTREE_vector_index = -1;
  LTREE_recursively_enumerate_cps(tree->root, cps);
  
  return(cps);
}





/*
 * LTREE_enumerate_ips() - 
 *
 * Take the root of a tree and traverse the tree, building an
 * exhaustive list of valid insertion points in the tree.
 *
 */
LTREE_treenode **
LTREE_enumerate_ips(LTREE_tree *tree,
		    long int nips) {

  LTREE_treenode **ips;

  ips = (LTREE_treenode **)calloc(nips, sizeof(LTREE_treenode *));
  LTREE_vector_index = 0;
  LTREE_recursively_enumerate_ips(tree->root, ips);
  
  return(ips);
}




/*
 * LTREE_recursively_enumerate_cps() - 
 *
 * Recursive function for traversing the tree and building a list of 
 * all internal nodes, which are the valid crossover points (except the root)
 *
 */
void
LTREE_recursively_enumerate_cps(LTREE_treenode *node,
				LTREE_treenode **cps) {

  /* if this is the root node and we don't add it to the cps list */
  if(LTREE_vector_index == -1) {

    if(node->left) {
      LTREE_vector_index++;
      LTREE_recursively_enumerate_cps(node->left, cps);
    }
    
    if(node->right) {
      LTREE_vector_index++;
      LTREE_recursively_enumerate_cps(node->right, cps);
    }

  } else {

    cps[LTREE_vector_index] = node;

    if(node->left) {
      LTREE_vector_index++;
      LTREE_recursively_enumerate_cps(node->left, cps);
    }
    
    if(node->right) {
      LTREE_vector_index++;
      LTREE_recursively_enumerate_cps(node->right, cps);
    }
  }
  
  return;
}







/*
 * LTREE_recursively_enumerate_ips() - 
 *
 * Recursive function for traversing the tree and building a list of 
 * all internal nodes, which are the valid insertion points
 *
 */
void
LTREE_recursively_enumerate_ips(LTREE_treenode *node,
				LTREE_treenode **ips) {

  ips[LTREE_vector_index] = node;

  if(node->left) {
    LTREE_vector_index++;
    LTREE_recursively_enumerate_cps(node->left, ips);
  }
    
  if(node->right) {
    LTREE_vector_index++;
    LTREE_recursively_enumerate_cps(node->right, ips);
  }

  return;
}



/*
 * LTREE_print_nodelist() - 
 *
 * Print the crossover points or insertion point lists (for debugging)
 *
 */
void
LTREE_print_nodelist(LTREE_treenode **list,
		     int size) {
  
  int i;
  
  for(i=0;i<size;i++) {
    printf("%d --> %p\n", i, list[i]);
  }
  
  return;
}





/*
 * LTREE_nondestructive_join() - 
 *
 * Take a subtree from parent1 (cp),
 * the compressed parent2 tree, and an insertion point in parent2
 * and construct a new tree without manipulating either parent1 or parent2
 */
LTREE_tree *
LTREE_nondestructive_join(LTREE_treenode *cp,
			  LTREE_tree *parent2,
			  char *desc) {
  
  LTREE_treenode *parent2_dup;
  LTREE_treenode *cp_dup;
  LTREE_treenode *ip;
  LTREE_treenode *new_node, *op;
  LTREE_tree *child;
  
  // duplicate the parent1 subtree and collapsed parent2 subtree
  cp_dup      = LTREE_dup_treenodes(cp);
  parent2_dup = LTREE_dup_treenodes(parent2->root);

  // find the pointer to the insertion point in the context of the dup'd parent2
  ip = LTREE_find_node_from_desc(parent2_dup, desc);

  // allocate child
  child         = (LTREE_tree *)calloc(1, sizeof(LTREE_tree));
  child->seqset = parent2->seqset;
  child->root   = parent2_dup;
  
  // now we have copies of the parent1 subtree, the collapsed parent2 subtree, 
  // and a pointer to where to insert the parent1 subtree into parent2.
  
  op = ip->parent;

  // allocate a new node which will join the trees
  new_node              = (LTREE_treenode *)calloc(1, sizeof(LTREE_treenode));
  new_node->left        = ip;
  new_node->right       = cp_dup;
  new_node->parent      = op;
  new_node->alignment   = NULL;
  new_node->seq_id      = -1;
  new_node->remove_flag = 0;

  ip->parent     = new_node;
  cp_dup->parent = new_node;
  
  if(ip == parent2_dup) {
    child->root = new_node;
  } else {
    if(op->left == ip) {
      op->left = new_node;
    } else if(op->right == ip) {
      op->right = new_node;
    } else {
      fprintf(stderr, "FUCKED\n");
      exit(-1);
    }
  }

  LTREE_print_tree(child);
  printf(" ;\n");

  return(child);
}





/*
 * LTREE_build_tree_desc() -
 *
 * Take a pointer-based binary tree and summarize it into 
 * a pointer-independent representation
 *
 */
LTREE_tree_desc * 
LTREE_build_tree_desc(LTREE_tree *tree) {

  LTREE_tree_desc *desc;
  char *buf;
  long int num_nodes, ntaxa;
  int i;
  
  desc       = (LTREE_tree_desc *)calloc(1, sizeof(LTREE_tree_desc));
  num_nodes  = LTREE_count_nodes(tree);
  ntaxa      = tree->seqset->num;
  desc->desc = (char **)calloc(num_nodes, sizeof(char *));
  desc->num  = 0;
  for(i=0;i<num_nodes;i++) {
    desc->desc[i] = (char *)calloc(ntaxa+1, sizeof(char));
  }
  buf = (char *)calloc(num_nodes+1, sizeof(char));
  
  LTREE_recursively_build_tree_desc(tree->root, desc, "");
  desc->num++;

  return(desc);
}



/*
 * LTREE_recursively_build_tree_desc() - 
 *
 * Recursively traverse a binary tree and assemble a 
 * pointer-independent tree description for every
 * node in the tree.
 *
 */
void
LTREE_recursively_build_tree_desc(LTREE_treenode *node,
				  LTREE_tree_desc *desc,
				  char *level_str) {
  if(node->left) {
    desc->num++;
    strcat(desc->desc[desc->num], level_str);
    strcat(desc->desc[desc->num], "L");
    LTREE_recursively_build_tree_desc(node->left, desc, desc->desc[desc->num]);
  }

  if(node->right) {
    desc->num++;
    strcat(desc->desc[desc->num], level_str);
    strcat(desc->desc[desc->num], "R");
    LTREE_recursively_build_tree_desc(node->right, desc, desc->desc[desc->num]);
  }

  return;
}



/*
 * LTREE_free_tree_desc() - 
 * 
 * Take a tree description structure and free it
 *
 */
void
LTREE_free_tree_desc(LTREE_tree_desc *desc) {

  unsigned long int i;
  
  if(desc) {
    if(desc->desc) {
      for(i=0;i<desc->num;i++) {
	if(desc->desc[i]) {
	  free(desc->desc[i]);
	}
      }
    
      free(desc->desc);
    }
    
    free(desc);
  }
  
  return;
}


/*
 * LTREE_count_nodes() - 
 *
 * Take a tree and count the number of nodes in the tree
 * including root, leaves, and all internal nodes
 *
 */
unsigned long int
LTREE_count_nodes(LTREE_tree *tree) {

  LTREE_global_node_count = 0;
  LTREE_recursively_count_nodes(tree->root);
  
  return(LTREE_global_node_count);
}




/*
 * LTREE_recursively_count_nodes() - 
 *
 * Recursively traverse a binary tree and return the
 * number of nodes in the tree
 *
 */
void
LTREE_recursively_count_nodes(LTREE_treenode *node) {

  if(node) {
    LTREE_global_node_count++;
  }
  
  if(node->left) {
    LTREE_recursively_count_nodes(node->left);
  }
  
  if(node->right) {
    LTREE_recursively_count_nodes(node->right);
  }
  
  return;
}



/*
 * LTREE_print_tree_desc() - 
 *
 * Pretty-print a tree description
 *
 */
void
LTREE_print_tree_desc(LTREE_tree_desc *desc) {
  
  unsigned long int i;
  
  for(i=0;i<desc->num;i++) {
    printf("%ld] %s\n", i, desc->desc[i]);
  }
  
  return;
}





/*
 * LTREE_find_node_from_desc() - 
 *
 * Find the pointer to the node described by the description
 *
 */
LTREE_treenode *
LTREE_find_node_from_desc(LTREE_treenode *treenode,
			  char *desc) {

  LTREE_treenode *ptr;
  char c;
  long int i;
  
  ptr = treenode;
  i = 0;
  c = desc[0];
  while(c) {
    
    switch(c) {
    case 'L':
      ptr = ptr->left;
      break;

    case 'R':
      ptr = ptr->right;
      break;
      
    default:
      fprintf(stderr, "ERROR FINDING NODE FROM DESC\n");
      exit(-1);
    }
    
    i++;
    c = desc[i];
  }

  return(ptr);
}
			
			

/* 
 * LTREE_canonize_tree() - 
 *
 * Take a tree and re-order the branches in order to canonicalize it
 *
 */
void
LTREE_canonize_tree(LTREE_tree *tree) {

  LTREE_recursively_canonize_tree(tree->root);
  LTREE_recursively_clear_internal_tree_seqid(tree->root);
  
  return;
}



/* 
 * LTREE_recursively_canonize_tree() - 
 *
 * The recursive node-based version of LTREE_canonize_tree()
 *
 */
void
LTREE_recursively_canonize_tree(LTREE_treenode *node) {
  
  return;
}



/*
 * LTREE_recursively_clear_internal_tree_seqid() - 
 *
 * Recursively traverse a tree and reset all internal nodes
 * and clear seqids = -1.  Leave terminal nodes alone.
 *
 */
void
LTREE_recursively_clear_internal_tree_seqid(LTREE_treenode *node) {
  
  if(node) {
    if ((node->left) && (node->right)) {
      node->seq_id = -1;
    }
    
    if(node->left) {
      LTREE_recursively_clear_internal_tree_seqid(node->left);
    }
    
    if(node->right) {
      LTREE_recursively_clear_internal_tree_seqid(node->right);
    }
  }
  
  return;
}
