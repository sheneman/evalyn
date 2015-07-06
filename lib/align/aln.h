/*
 * aln.h
 *
 * $Id: aln.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2003 Luke Sheneman.  All rights reserved.
 *
 *
 * Routines for writing alignments in CLUSTAL W formatted ALN format
 * 
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */

#ifndef _INC_ALN_H_
#define _INC_ALN_H_ 1

#include "align.h" 


/* some function prototypes */
void
ALIGN_write_alignment_ALN(ALIGN_alignment *alignment,
			  ALIGN_seqset *seqset,
			  char *filename);

#endif /* _INC_ALN_H_ */





