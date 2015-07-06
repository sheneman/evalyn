#include <stdio.h>
#include <stdlib.h>
#include "prng.h"

int 
main(int argc, 
     char **argv) {
  
  unsigned int seed;
  int number;
  unsigned int unumber;
  float flt_number;
  int i;

  MT_state saved_state;


  printf("Testing some stuff.\n");
  
  seed = 100;
  MT_init_genrand(seed);
  
  /* generate a number */
  number = MT_get_sint32();
  printf("Number: %d\n", number);
  
    flt_number = MT_get_real1();

  /* save state */  
  printf("SAVINGS STATE\n");
  MT_get_state(&saved_state);

  /* generate a number */
  number = MT_get_sint32();
  printf("Number: %d\n", number);


  /* generate a number */
  number = MT_get_sint32();
  printf("Number: %d\n", number);


  /* generate a number */
  number = MT_get_sint32();
  printf("Number: %d\n", number);

  /* generate a float here */
  flt_number = MT_get_real1();
  printf("FLT number: %f\n", flt_number);

  /* generate a float here */
  flt_number = MT_get_real1();
  printf("FLT number: %f\n", flt_number);

  /* generate a float here */
  flt_number = MT_get_real1();
  printf("FLT number: %f\n", flt_number);

  /* generate a float here */
  flt_number = MT_get_real1();
  printf("FLT number: %f\n", flt_number);

  /* generate a float here */
  flt_number = MT_get_real1();
  printf("FLT number: %f\n", flt_number);

  /* restore state here */
  printf("RESTORING STATE\n");
  MT_set_state(&saved_state);

  /* generate a number */
  number = MT_get_sint32();
  printf("Number: %d\n", number);


  /* generate a number */
  number = MT_get_sint32();
  printf("Number: %d\n", number);


  /* generate a number */
  number = MT_get_sint32();
  printf("Number: %d\n", number);

  /* generate a float here */
  flt_number = MT_get_real1();
  printf("FLT number: %f\n", flt_number);


  /* generate a float here */
  flt_number = MT_get_real1();
  printf("FLT number: %f\n", flt_number);

  /* generate a float here */
  flt_number = MT_get_real1();
  printf("FLT number: %f\n", flt_number);

  /* generate a float here */
  flt_number = MT_get_real1();
  printf("FLT number: %f\n", flt_number);

  /* generate a float here */
  flt_number = MT_get_real1();
  printf("FLT number: %f\n", flt_number);

  for(i=0;i<100;i++) {
    unumber = MT_get_range_uint32(5);
    printf("%d\n", unumber);
  }

  exit(0);
}











