/*
 * msf.c
 *
 * $Id: msf.c,v 1.1.1.1 2006/05/09 21:00:46 cvsuser Exp $
 *
 *
 * Copyright (c) 2005 Luke Sheneman.  All rights reserved.
 *
 *
 * A routine for outputting MSF-formatted alignments.  The MSF format
 * is a common alignment format, used by BAliBASE and the GCG (Wisconsin)
 * package.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "align.h"
#include "msf.h"


/*************************************
 *                                   *
 * BEGIN STATIC FUNCTION PROTOTYPES  *
 *                                   *
 *************************************/


static int
ALIGN_get_seq_index(ALIGN_seqset *seqset,
		    char *title);

static void
ALIGN_clean_seq(char *sequence);

/***********************************
 *                                 *
 * END STATIC FUNCTION PROTOTYPES  *
 *                                 *
 ***********************************/





/*
 * ALIGN_read_alignment_MSF() - 
 *
 * Read an MSF alignment into memory and convert it to a seqset
 * For now, don't worry about checksums or other annotation, just
 * grab sequence names and strings
 *
 */
int
ALIGN_read_MSF(ALIGN_seqset *seqset,
	       char *filename) {
  
  FILE *fp;
  char line[MAX_SEQ_LENGTH], buf[MAX_SEQ_LENGTH];
  char **seq_names;
  char *string, *c=NULL, ck;
  int i, j, index;
  int seqnum;

  
  /* open MSF alignment file here */
  fp = fopen(filename, "r");
  if(!fp) {
    fprintf(stderr, "ALIGN: Could not read MSF file %s\n", filename);
    exit(-1);
  }
  
  seqnum = ALIGN_MSF_INITIAL_SEQNUM;
  seq_names = (char **)calloc(seqnum, sizeof(char *));
  for(i=0;i<seqnum;i++) {
    seq_names[i] = (char *)calloc(MAX_SEQ_TITLE_LENGTH, sizeof(char));
  }

  /* read file and grab sequence names until you get to the demarcation (double slash) */
  seqset->num = i = 0;
  while(!feof(fp)) {

    fgets(line, MAX_SEQ_LENGTH, fp);

    string = (char *)strtok(line, " \t");

    /* read lines here */
    if(!strcmp(string, "Name:") ||
       !strcmp(string, "name:") ||
       !strcmp(string, "NAME:") ) {

      string = (char *)strtok(NULL, " \t");
      if(string) {
	if(i >= seqnum) {
	  seqnum *= 2;
	  seq_names = (char **)realloc(seq_names, seqnum*sizeof(char *));
	  for(j=i;j<seqnum;j++) {
	    seq_names[j] = (char *)calloc(j, MAX_SEQ_TITLE_LENGTH);
	  }
	}
	strcpy(seq_names[i++], string);
      }
    }

    if(line[0] == '/' && line[1] == '/') {
      break;
    }
  }
  seqset->num = i;
  
  /* if there are no sequences, this is not an MSF formatted file */
  if(!seqset->num) {
    fprintf(stderr, "ALIGN: File %s does not appear to be a valid MSF file\n", filename);
    exit(-1);
  }

  /* allocate the sequences in the seqset */
  seqset->seq = (ALIGN_seq *)calloc(seqset->num, sizeof(ALIGN_seq));
  for(i=0;i<seqset->num;i++) {

    /* allocate the data */
    seqset->seq[i].data = (char *)calloc(MAX_SEQ_LENGTH, sizeof(char));

    /* assign names */
    strcpy(seqset->seq[i].title, seq_names[i]);
  }
  
  /* parse the rest of the file here and put the sequences in the right place */
  i = 0;
  while(!(feof(fp))) {

    fgets(line, sizeof(line), fp);
    
    strcpy(buf, line);
    string = (char *)strtok(line, " \t");
    if(string && string[0] != 10 && string[0] != 13) {
      index = ALIGN_get_seq_index(seqset, string);

      /* find first whitespace here */
      for(i=0;i<strlen(buf);i++) {
	if(buf[i] == ' ' ||
	   buf[i] == '\t') {
	  c = &buf[i];
	  break;
	}
      }

      /* remove whitespace, gaps, periods, blah */
      ALIGN_clean_seq(c);
      
      /* concatenate the cleaned up sequence to the seqset */
      strcat(seqset->seq[index].data, c);
    }
  }
  
  /* lets assign lengths to the sequences and convert to indices */
  for(i=0;i<seqset->num;i++) {
    seqset->seq[i].length = strlen(seqset->seq[i].data);
    
    /* convert to the canonical internal representation */
    for(j=0;j<seqset->seq[i].length;j++) {
      ck = seqset->seq[i].data[j];
      seqset->seq[i].data[j] = ALIGN_seqset_lookup_sym(seqset, ck);
      if(seqset->seq[i].data[j] == -1) {
	fprintf(stderr, "ALIGN: Seq. file contains symbol \"%c\" that is not in substitution matrix.\n", ck);
	exit(-1);
      }
    }
  }
  
  /* free our sequence names here */
  if(seq_names) {
    for(i=0;i<seqnum;i++) {
      if(seq_names[i]) {
	free(seq_names[i]);
      }
    }
    
    free(seq_names);
  }

  /* close the msf file pointer here */
  fclose(fp);
  
  //  return(ALIGN_seqset_type(seqset));
  return(ALIGN_PROTEIN);  /* for now */
}






/*
 *
 * ALIGN_clean_seq() - 
 *
 * Remove whitespace, gaps, and other codes and uppercase sequences
 *
 */
static void
ALIGN_clean_seq(char *sequence) {

  int i, j;
  char *buf, c;
  
  buf = (char *)calloc(strlen(sequence)+1, sizeof(char));
  
  j=0;
  for(i=0;i<strlen(sequence);i++) {

    c = toupper(sequence[i]);

    if(c >= 'A' && c <= 'Z') {
      buf[j++] = c;
    }
  }
  buf[j]='\0';

  strcpy(sequence, buf);

  free(buf);

  return;
}







/* 
 *
 * ALIGN_get_seq_index() - 
 *  
 * Look up the sequence index given a title
 * it should be okay to do a simple linear search here (for now)
 * 
 * -1 == BAD
 *
 */
static int
ALIGN_get_seq_index(ALIGN_seqset *seqset,
		    char *title) {
  
  int i;

  /* simple linear search for now */
  for(i=0;i<seqset->num;i++) {
    if(strcmp(seqset->seq[i].title, title) == 0 ) {
      return(i);
    }
  }
  
  return(-1);
}




/*
 *
 * ALIGN_write_alignment_MSF() -
 *
 * Output the given alignment in MSF format
 *
 * the argument "type" should be either ALIGN_DNA or ALIGN_PROTEIN
 *
 */
void
ALIGN_write_alignment_MSF(ALIGN_alignment *alignment,
			  ALIGN_seqset *seqset,
			  char *outfile,
			  int type) {
  
  FILE *fp;
  time_t t;
  int i, j, segments, checksum;
  char c, datestr[64];
  char *buf;


  printf("In ALIGN_write_alignment_MSF()\n");

  /* open outfile for writing */
  fp = fopen(outfile, "w");
  if(!fp) {
    fprintf(stderr, "ALIGN: Could not write MSF alignment %s\n", outfile); 
    return;
  }
  

  /* get current localtime and convert to a usable string */
  time(&t);
  sprintf(datestr, "%s", (char *)asctime(localtime(&t)));
  datestr[strlen(datestr)-1] ='\0';

  /* compute the project checksum here */
  checksum = 0;
  for(i=0;i<seqset->num;i++) {
    checksum += ALIGN_MSF_checksum(seqset->seq[i].data);
  }
  checksum = checksum % 10000L;
  
  /* print MSF header here */
  fprintf(fp, "MSF of: ALIGNMENT from: 1 to %d\n", alignment->n);
  fprintf(fp, "%s  MSF: %d  Type: %c  %s  Check: %d ..\n\n", 
	  outfile,
	  alignment->n,
	  ((type==ALIGN_PROTEIN)?'P':'D'),
	  datestr,
	  checksum);

  
  /* allocate a tmp buffer */
  buf = (char *)calloc(alignment->n+1, sizeof(char));
  
  /* print the summary of the sequences  */
  printf("ALIGNMENT->N = %d, ALIGNMENT->K = %d\n", alignment->n, alignment->k);
  for(i=0;i<seqset->num;i++) {
    
    /* construct alignment here */
    for(j=0;j<alignment->n;j++) {
      buf[j] = alignment->text[ALIGN_index(alignment, j, i)];
    }
    buf[j]='\0';

    fprintf(fp, "Name: %s  Len:  %d  Check:  %d  Weight:  1.00\n",
	    seqset->seq[alignment->seq_ids[i]].title,
	    alignment->n,
	    ALIGN_MSF_checksum(buf));
  }
  
  /* print separator which indicates the beginning of the alignment */
  fprintf(fp, "\n//\n\n");
  
  /* print alignment here */
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
  
  free(buf);

  return;
}




/*
 * ALIGN_MSF_checksum() - 
 *
 * Compute the sequence checksum for the given sequence 
 *
 * This function derived from GCG documentation about the MSF format
 *
 */
int
ALIGN_MSF_checksum(char *sequence) {

    int len, position, chk = 0;

    len = strlen(sequence);
    for (position = 0; position < len; position ++)
        chk = (chk + (position % 57 + 1) * ( toupper(sequence[position]) ) ) %10000;

    return chk;
}


