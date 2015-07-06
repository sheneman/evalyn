/*
 * rnjhooks.h
 *
 * $Id: rnjhooks.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2005 Luke Sheneman.  All rights reserved.
 *
 *
 */

#ifndef _INC_RNJHOOKS_H_
#define _INC_RNJHOOKS_H_ 1

#include <align.h>
#include <nj_fasta.h>

#include "ltree.h"


/* some function prototypes */
LTREE_tree *
LTREE_create_rnj_tree(ALIGN_scoring_system *ss,
		      LTREE_population *pop,
		      ALIGN_seqset *seqset);

LTREE_tree *
LTREE_rnjtree_to_ltree(NJ_TREE *rnj_tree,
		       NJ_DMAT *dmat,
		       ALIGN_seqset *seqset);
LTREE_treenode *
LTREE_recurse_rnjtree_to_ltree(NJ_TREE *rnj_tree,
			       NJ_DMAT *dmat,
			       ALIGN_seqset *seqset);

NJ_alignment *
LTREE_ev_alignment_to_rnj(ALIGN_seqset *seqset,
			  ALIGN_alignment *ev_alignment);


#endif /* _INC_RNJHOOKS_H_ */



