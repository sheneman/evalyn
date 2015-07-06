/*
 * msf.h
 *
 * $Id: msf.h,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2003 Luke Sheneman.  All rights reserved.
 * 
 *
 * A routine for outputting MSF-formatted alignments
 *
 * The specification for the MSF format can be found here:
 *
 *   
 *  http://www.embl-heidelberg.de/predictprotein/Dexa/optin_msfDes.html
 *
 * 
 * Luke Sheneman
 * sheneman@cs.uidaho.edu
 *
 */

#ifndef _INC_MSF_H_
#define _INC_MSF_H_ 1

#define ALIGN_MSF_INITIAL_SEQNUM 2


#include "align.h"



/* some function prototypes */

int
ALIGN_MSF_checksum(char *sequence);

int
ALIGN_read_MSF(ALIGN_seqset *seqset,
	       char *filename);

void
ALIGN_write_alignment_MSF(ALIGN_alignment *alignment,
			  ALIGN_seqset *seqset,
			  char *outfile,
			  int type);

#endif /* _INC_MSF_H_ */











