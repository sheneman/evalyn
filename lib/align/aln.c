/*
 * aln.c
 *
 * $Id: aln.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2003 Luke Sheneman.  All rights reserved.
 *
 *
 * Routines for writing alignments in CLUSTAL W formatted ALN format
 *
 * 
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */


#include <stdio.h>
#include <string.h>

#include "align.h"
#include "aln.h"



/*
 * ALIGN_write_alignment_ALN() -
 *
 * Write a ClustalW-style ALN file
 *
 */
void
ALIGN_write_alignment_ALN(ALIGN_alignment *alignment,
			  ALIGN_seqset *seqset,
			  char *filename) {
  
  FILE *fp;
  unsigned int i, j, c;
  unsigned int segments;

  
  fp = fopen(filename, "w");
  if(!fp) {
    fprintf(stderr, "Error:  Could not open output file %s\n", filename);
    return;
  }
  
  fprintf(fp, "Alignment From Evalyn.  Evalyn Copyright (c) 2003, Luke Sheneman. \n\n\n");
  
  segments = (alignment->n)/60 + 1;

  for(c=0;c<segments;c++) {

    for(j=0;j<alignment->k;j++) {   /* sequences */
      
      if(seqset) {
	fprintf(fp, "%s       ", seqset->seq[alignment->seq_ids[j]].title);
      } else {
	fprintf(fp, "%d       ", j);
      }

      for(i=c*60;i<c*60+60;i++) {     /* length    */      
	if(i<alignment->n) {
	  if(alignment->text[ALIGN_index(alignment, i, j)] == GAP_CHAR) {
	    fprintf(fp, "%c", GAP_CHAR);
	  } else {
	    fprintf(fp, "%c", alignment->symtab->syms[(int)(alignment->text[ALIGN_index(alignment, i, j)])]);
	  }
	}
      }
      
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  return;
}


