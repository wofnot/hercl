/** idxstr.c */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "mystring.h"

/********************************************************//**
   Input: index followed by (non-empty) string
  Output: character at the specified index
*/
int main(int argc, char *argv[])
{
  int s[256];
  int tests = NUM_TRAIN_TEST;
  int i,j,k,n;
  
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
    printf(",%d\n", j ); // index of character

    print_string( stdout,k,s );
    printf("\n");

    print_character( stdout,s[j] ); // the character at that index
    printf("\n");

    if( n == tests ) {
      printf("\\\n");
    }
  }

  return 0;
}
