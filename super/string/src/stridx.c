/** stridx.c */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "mystring.h"

/********************************************************//**
   Input: a (non-empty) string followed by a list of indices
  Output: the list of characters at the specified indices
*/
int main(int argc, char *argv[])
{
  int s[256];
  int t[256];
  int tests = NUM_TRAIN_TEST;
  int i,j,k,l,n;
  
  if( argc > 1 ){
    tests = atoi(argv[1]);
  }

  srandom((argc < 3 ? time(NULL) : atoi(argv[2])));

  printf("2\n"); // inputs_per_item

  for( n=1; n <= 2*tests; n++ ) {

    j = -1 + random_geometric( 4,32 );
    k =  j + random_geometric( 4,32 );

    for( i=0; i < k; i++ ) {
      s[i] = random_character();
    }
    print_string( stdout,k,s );
    printf("\n");

    if( j == 0 ) {
      printf(",\n");
    }
    else {
      for( i=0; i < j; i++ ) {
        l = random()% k;
        printf(",%d",l);
        t[i] = s[l];
      }
      printf("\n");
    }

    print_string( stdout,j,t );
    printf("\n");

    if( n == tests ) {
      printf("\\\n");
    }
  }

  return 0;
}
