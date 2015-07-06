/*
 * ltree.h
 *
 * $Id: ltree.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2005 Luke Sheneman.  All rights reserved.
 *
 */

#ifndef _INC_LTREE_H_
#define _INC_LTREE_H_ 1


#include <sys/param.h>

#include <prng.h> /* mersenne twister stuff */



#define SAVE_ALIGNMENT  100
#define PURGE_ALIGNMENT 101

/* some defaults */
#define DEFAULT_DNA_MATRIX "./def_dna_matrix.txt"
#define DEFAULT_PRO_MATRIX "./def_pro_matrix.txt"
#define DEFAULT_ITERATIONS      1000
#define DEFAULT_MUTATION_RATE   0.001
#define DEFAULT_POPULATION_SIZE 100

#define DEFAULT_CONVERGE_STATIONARITY 1000 

#define DEFAULT_DNA_GAPOPEN      (-15.0)
#define DEFAULT_DNA_GAPEXTEND    (-6.66)

#define DEFAULT_PRO_GAPOPEN      (-10.0)
#define DEFAULT_PRO_GAPEXTEND    (-0.1)
   
#define LTREE_MSF     100
#define LTREE_FASTA   101
#define LTREE_ALN     102

#define LTREE_SPS       100
#define LTREE_PARSIMONY 101

#define LTREE_MIN_POPSIZE 4

/* the possible search algorithms */
#define LTREE_ALGORITHM_GA                  100
#define LTREE_ALGORITHM_RANDOM_WALK         101
#define LTREE_ALGORITHM_HILL_CLIMB          102
#define LTREE_ALGORITHM_SIMULATED_ANNEALING 103

#define LTREE_DEFAULT_ALGORITHM LTREE_ALGORITHM_GA


#define DEFAULT_INFILETYPE  LTREE_FASTA
#define DEFAULT_OUTFILETYPE LTREE_ALN




/* some datatypes */
typedef struct LTREE_STRUCT_cl_args {
  
  int verbose_flag;
  int quiet_flag;
  int debug_flag;

  int dna_flag;
  int protein_flag;

  int test_flag;
  int parsimony_flag;
  int sps_flag;

  /* define our search algorithm */
  int search_algorithm;

  int rnj_flag;
  
  int converge_flag;
  int stationarity;

  int scoring_method;

  int infiletype;
  int outfiletype;

  int gap_extend_set;
  int gap_open_set;

  long int seed;
  
  int seq_type;

  char matrix_filename[MAXPATHLEN];
  char input_filename[MAXPATHLEN];
  char output_filename[MAXPATHLEN];
  char tree_filename[MAXPATHLEN];
  
  int population_size;
  int iterations;

  double mutation_rate;
  double gap_open;
  double gap_extend;
  
  
} LTREE_cl_args;



typedef struct LTREE_STRUCT_treenode {

  struct LTREE_STRUCT_treenode *left;
  struct LTREE_STRUCT_treenode *right;
  struct LTREE_STRUCT_treenode *parent;
  
  ALIGN_alignment *alignment;

  short int seq_id;
  char remove_flag;

  /* parsimony related elements */
  unsigned int parsimony; 
  char pbuf[MAX_SYMS];  /* buffer used for parsimony calculations */
  char pbuf_p;          /* number of actual symbols in pbuf */

} LTREE_treenode;



typedef struct LTREE_STRUCT_tree {
  
  LTREE_treenode *root;

  ALIGN_seqset *seqset;

  MT_state *prng_state;
  
  double fitness;

} LTREE_tree;



typedef struct LTREE_STRUCT_population {

  unsigned short int size;  

  LTREE_tree **trees;      /* pop of trees */
  
  MT_state best_prng_state;

} LTREE_population;



typedef struct LTREE_STRUCT_taxa_list {
  
  struct LTREE_STRUCT_taxa_list *next;
  short int seq_id;

} LTREE_taxa_list;




typedef struct LTREE_STRUCT_tree_desc {

  char **desc;
  long int num;

} LTREE_tree_desc;





/* some prototypes for inline functions */

/* some prototypes for exposed functions */
LTREE_population *
LTREE_evolve(LTREE_cl_args *cl_args,
	     ALIGN_seqset *seqset,
	     ALIGN_scoring_system *ss);

LTREE_population *
LTREE_init_population(ALIGN_seqset *seqset, 
		      unsigned long int pop_size);

LTREE_tree *
LTREE_init_tree(ALIGN_seqset *seqset);

void
LTREE_sort_population(LTREE_population *population);


LTREE_tree *
LTREE_crossover(LTREE_tree *parent1,
		LTREE_tree *parent2);

LTREE_tree *
LTREE_mutate(LTREE_tree *child);

void
LTREE_selection(LTREE_population *population,
		unsigned short int *parent1,
		unsigned short int *parent2,
		unsigned short int *child);

void
LTREE_print_tree(LTREE_tree *tree);

void
LTREE_free_tree(LTREE_tree *tree);

double
LTREE_eval_fitness(ALIGN_seqset *seqset, 
		   ALIGN_scoring_system *ss,
		   LTREE_tree *tree,
		   char delete_alignment);
double
LTREE_align_tree(LTREE_tree *tree,
		 ALIGN_scoring_system *ss);

void
LTREE_pop_stats(LTREE_population *population,
		double *best_fitness,
		double *mean_fitness);

void
LTREE_update_population(LTREE_population *population,
			LTREE_tree *child_tree,
			unsigned short int child,
			double fitness);


void
LTREE_test_alignment(void);

void 
LTREE_free_population(LTREE_population *population);

void
LTREE_compress_tree(LTREE_tree *tree);

void
LTREE_write_tree(LTREE_tree *tree,
		 FILE *fp);


unsigned int 
LTREE_tree_depth(LTREE_treenode *node,
		 unsigned int depth);

unsigned int
LTREE_tree_size(LTREE_treenode *node);

void
LTREE_inject_rnj_tree(ALIGN_scoring_system *ss,
		      LTREE_population *pop,
		      ALIGN_seqset *seqset);
void
LTREE_write_alignment(ALIGN_alignment *alignment,
		      ALIGN_seqset *seqset,
		      LTREE_cl_args *cl_args);

void
LTREE_dump_population(char *filename,
		      LTREE_population *population,
		      int index,
		      LTREE_tree *parent1,
		      LTREE_tree *parent2,
		      LTREE_tree *child);

LTREE_population *
LTREE_lxo_step(LTREE_tree *source_tree);

LTREE_treenode **
LTREE_enumerate_cps(LTREE_tree *tree,
		    long int ncps);

LTREE_treenode **
LTREE_enumerate_ips(LTREE_tree *tree,
		    long int nips);

void
LTREE_recursively_enumerate_cps(LTREE_treenode *node,
				LTREE_treenode **cps);

void
LTREE_recursively_enumerate_ips(LTREE_treenode *node,
				LTREE_treenode **ips);

void
LTREE_print_nodelist(LTREE_treenode **list,
		     int size);

long int
LTREE_count_taxa_list(LTREE_taxa_list *head);

LTREE_tree *
LTREE_nondestructive_join(LTREE_treenode *cp,
			  LTREE_tree *parent2,
			  char *desc);

LTREE_tree_desc *
LTREE_build_tree_desc(LTREE_tree *tree);

void
LTREE_recursively_build_tree_desc(LTREE_treenode *node,
				  LTREE_tree_desc *desc,
				  char *level_str);

void
LTREE_print_tree_desc(LTREE_tree_desc *desc);

unsigned long int
LTREE_count_nodes(LTREE_tree *tree);

void
LTREE_recursively_count_nodes(LTREE_treenode *node);

LTREE_treenode *
LTREE_find_node_from_desc(LTREE_treenode *treenode,
			  char *desc);

void
LTREE_free_tree_desc(LTREE_tree_desc *desc);

void
LTREE_canonize_tree(LTREE_tree *tree);

void
LTREE_recursively_canonize_tree(LTREE_treenode *node);

void
LTREE_recursively_clear_internal_tree_seqid(LTREE_treenode *node);

#endif   /* _INC_LTREE_H_ */









