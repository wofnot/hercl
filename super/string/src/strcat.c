/** catstr.c */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "mystring.h"

void generate_strings( int j, int k );


/********************************************************//**
   Input: two strings s[] and t[]
  Output: the concatenated string s[] followed by t[]
*/
int main( int argc, char **argv )
{
  int tests = NUM_TRAIN_TEST;
  int j,k,n;

  if (argc > 1) {
    tests = atoi(argv[1]);
  }

  srandom((argc < 3 ? time(NULL) : atoi(argv[2])));

  printf("2\n"); // inputs_per_item

  for( n=1; n <= 2*tests; n++ ) {

    j = -2 + random_geometric( 4,32 )
           + random_geometric( 4,32 );
    k = -2 + random_geometric( 4,32 )
           + random_geometric( 4,32 );

    generate_strings( j,k );

    if( n == tests ) {
      printf("\\\n");
    }
  }
  return 0;
}

/********************************************************//**
   Generate two random strings of printable characetrs,
   of length j and k. Print each string separately,
   and then their concatenation.
*/
void generate_strings( int j, int k )
{
  int s[256];
  int t[256];
  int i;

  for( i=0; i < j; i++ ) {
    s[i] = random_character();
  }
  for( i=0; i < k; i++ ) {
    t[i] = random_character();
  }
  if( j > 0 ) {
    print_string( stdout,j,s );
  }
  else {
    printf(",");
  }
  printf("\n");
  if( k > 0 ) {
    print_string( stdout,k,t );
  }
  else {
    printf(",");
  }
  printf("\n");
  if( j > 0 ) {
    print_string( stdout,j,s );
  }
  if( k > 0 ) {
    print_string( stdout,k,t );
  }
  if( j == 0 && k == 0 ) {
    printf(",");
  }
  printf("\n");
}
