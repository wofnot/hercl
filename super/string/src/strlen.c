/** strlen.c */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "mystring.h"


/********************************************************//**
   Input: a string
  Output: the length of that string
*/
int main( int argc, char **argv )
{
  int s[256];
  int tests = NUM_TRAIN_TEST;
  int j,k,n;
  
  if (argc > 1) {
    tests = atoi(argv[1]);
  }

  srandom((argc < 3 ? time(NULL) : atoi(argv[2])));

  printf("1\n"); // inputs_per_item

  for( n=1; n <= 2*tests; n++ ) {

    k = -2 + random_geometric( 4,32 )
           + random_geometric( 4,32 );

    for( j=0; j < k; j++ ) {
      s[j] = random_character();
    }
    if( k > 0 ) {
      print_string( stdout,k,s );
    }
    else {
      printf(",");
    }

    printf("\n");
    printf(",%d\n", k );

    if( n == tests ) {
      printf("\\\n");
    }
  }
  return 0;
}
