/*
 * sp.h
 *
 * $Id: sp.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2003 Luke Sheneman.  All rights reserved.
 * 
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */

#ifndef _INC_SP_H_
#define _INC_SP_H_ 1


#include "align.h"

/* function prototypes */
float
ALIGN_compute_sp(ALIGN_scoring_system *ss,
		 ALIGN_alignment *alignment);


#endif /* _INC_SP_H_ */











