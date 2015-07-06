/*
 * nj_main.c
 *
 * $Id: nj_main.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *****************************************************************************
 *
 * Copyright (c) 2005,  Luke Sheneman
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions 
 * are met:
 * 
 *  + Redistributions of source code must retain the above copyright 
 *    notice, this list of conditions and the following disclaimer. 
 *  + Redistributions in binary form must reproduce the above copyright 
 *    notice, this list of conditions and the following disclaimer in 
 *    the documentation and/or other materials provided with the 
 *    distribution. 
 *  + The names of its contributors may not be used to endorse or promote 
 *    products derived  from this software without specific prior 
 *    written permission. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.  
 *
 *****************************************************************************
 *
 * An implementation of the Relaxed Neighbor-Joining algorithm 
 *  of Evans, J., Sheneman, L., and Foster, J., 2005
 *
 *
 * AUTHOR:
 * 
 *   Luke Sheneman
 *   sheneman@cs.uidaho.edu
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include <prng.h>

#include "nj_main.h"

#include "nj_dist.h"
#include "nj_dmat.h"
#include "nj_fasta.h"
#include "nj_cmdargs.h"
#include "nj_common.h"
#include "nj_bootstrap.h"
#include "nj_clearcut.h"




/*
 * main() - 
 *
 * The entry point to the program.
 *
 */
int
main(int argc,
     char *argv[]) {

  NJ_DMAT *dmat = NULL;          /* The working distance matrix */
  NJ_DMAT *dmat_backup=NULL;     /* A backup distance matrix    */
  NJ_TREE *tree;                 /* The phylogenetic tree       */
  NJ_ARGS *nj_args;              /* Structure for holding command-line arguments */
  NJ_alignment *alignment=NULL;  /* Structure containing multiple sequence alignment */
  long int i;
  unsigned long int maxiter = 1; /* The number of iterations */


  /* some variables for tracking time */
  struct timeval tv;
  unsigned long long startUs, endUs;
  

  /* check and parse supplied command-line arguments */
  nj_args = NJ_handle_args(argc, argv);
  if(!nj_args) {
    fprintf(stderr, "Clearcut: Error processing command-line arguments.\n");
    exit(-1);
  }

  /* for verbose reporting, print the random number seed to stdout */
  if(nj_args->verbose_flag) {
    printf("PRNG SEED: %d\n", nj_args->seed);
  }

  /* Initialize Mersenne Twister PRNG */
  MT_init_genrand(nj_args->seed);


  switch(nj_args->input_mode) {

    /* If the input type is a distance matrix */
  case NJ_INPUT_MODE_DISTANCE:

    /* parse the distance matrix */
    dmat = NJ_parse_distance_matrix(nj_args);
    if(!dmat) {
      exit(-1);
    }

    break;

    /* If the input type is a multiple sequence alignment */
  case NJ_INPUT_MODE_ALIGNED_SEQUENCES:

    /* build a distance matrix from a multiple sequence alignment */
    if(!nj_args->bootstrap) {
      dmat = NJ_build_distance_matrix(nj_args);
      if(!dmat) {
	fprintf(stderr, "Clearcut: Failed to build distance matrix from alignment.\n");
	exit(-1);
      }
    } else {
      alignment = NJ_read_fasta(nj_args);
      if(!alignment) {
	fprintf(stderr, "Clearcut: Failed to parse alignment.\n");
	exit(-1);
      }
    }
    
    break;

  default:

    fprintf(stderr, "Clearcut: Could not determine how to process input\n");
    exit(-1);
  }

  /*
   * Output the computed distance matrix,
   *  if the user specified to do so.
   */
  if(nj_args->matrixout && !nj_args->bootstrap) {
    NJ_output_matrix(nj_args, dmat);
  }

  /* 
   * If we are going to generate multiple trees from
   * the same distance matrix, we need to make a backup 
   * of the original distance matrix.
   */
  if(nj_args->ntrees > 1) {

    dmat_backup = NJ_dup_dmat(dmat);
    maxiter = nj_args->ntrees;

  } else if(nj_args->bootstrap) {
    dmat = NJ_bootstrap_alignment(nj_args, alignment);
    maxiter = nj_args->bootstrap_replicates;
  }
  
  /* process <maxiter> trees */
  for(i=0;i<maxiter;i++) {
    
    /* 
     * If the user has specified matrix shuffling, we need 
     * to randomize the distance matrix
     */
    if(nj_args->shuffle) {
      NJ_shuffle_distance_matrix(dmat);
    }

    /* RECORD THE PRECISE TIME OF THE START OF THE NEIGHBOR-JOINING */
    gettimeofday(&tv, NULL);
    startUs = ((unsigned long long) tv.tv_sec * 1000000ULL)
      + ((unsigned long long) tv.tv_usec);

    
    /* 
     * Invoke either the Relaxed Neighbor-Joining algorithm (default)
     * or the "traditional" Neighbor-Joining algorithm 
     */
    if(nj_args->neighbor) {
      tree = NJ_neighbor_joining(nj_args, dmat);
    } else {
      tree = NJ_relaxed_nj(nj_args, dmat);
    }
  
    if(!tree) {
      fprintf(stderr, "Clearcut: Failed to construct tree.\n");
      exit(0);
    }

    /* RECORD THE PRECISE TIME OF THE END OF THE NEIGHBOR-JOINING */
    gettimeofday(&tv, NULL);
    endUs = ((unsigned long long) tv.tv_sec * 1000000ULL)
      + ((unsigned long long) tv.tv_usec);

    /* print the time taken to perform the neighbor join */
    if(nj_args->verbose_flag) {
      if(nj_args->neighbor) {
	fprintf(stderr, "NJ tree built in %llu.%06llu secs\n",
		(endUs - startUs) / 1000000ULL,
		(endUs - startUs) % 1000000ULL);
      } else { 
	fprintf(stderr, "RNJ tree built in %llu.%06llu secs\n",
		(endUs - startUs) / 1000000ULL,
		(endUs - startUs) % 1000000ULL);
      }
    }

    /* Output the neighbor joining tree here */
    NJ_output_tree(nj_args, tree, dmat, i);
    
    NJ_free_tree(tree);  /* Free the tree */
    NJ_free_dmat(dmat);  /* Free the working distance matrix */

    /* 
     * If we need to do another iteration, lets re-initialize 
     * our working distance matrix.
     */
    if(nj_args->bootstrap) {
      if(i<(nj_args->bootstrap_replicates-1)) {
	dmat = NJ_bootstrap_alignment(nj_args, alignment);
      }
    } else {
      if(nj_args->ntrees > 1 && i<(nj_args->ntrees-1) ) {
	dmat = NJ_dup_dmat(dmat_backup);
      }
    }
  }
  
  /* Free the backup distance matrix */
  if(nj_args->ntrees > 1) {
    NJ_free_dmat(dmat_backup);
  }

  /* If verbosity, describe where the tree output is */
  if(nj_args->verbose_flag) {
    if(nj_args->neighbor) {
      printf("NJ tree(s) in %s\n", nj_args->outfilename);
    } else {
      printf("Relaxed NJ tree(s) in %s\n", nj_args->outfilename);
    }
  }

  exit(0);
}


