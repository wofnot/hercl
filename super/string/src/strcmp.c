/** strcmp.c */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "mystring.h"

void generate_strings( int j, int k, int l );

/********************************************************//**
   Input: two strings
  Output: difference between first char where strings differ
*/
int main( int argc, char **argv )
{
  int j,k,l,n;
  
  int tests = NUM_TRAIN_TEST;

  if (argc > 1) {
    tests = atoi(argv[1]);
  }

  srandom((argc < 3 ? time(NULL) : atoi(argv[2])));

  printf("2\n");

  generate_strings( 0,0,0 );

  for( n=2; n <= 2*tests; n++ ) {

    do {
      j = -1 + random_geometric(4,32); // substring common to s[],t[]
      k = -1 + random_geometric(4,32); // substring unique to s[]
      l = -1 + random_geometric(4,32); // substring unique to t[]
      //j = random()% 5;
      //k = random()% 5;
      //l = random()% 5;
    } while( j == 0 && k == 0 && l == 0 );

    generate_strings( j,k,l );

    if( n == tests ) {
      printf("\\\n");
    }
  }

  return 0;
}

/********************************************************//**
   Generate two random strings s[] and t[] which share a
   common prefix of length j, followed by a distinct portion
   of length k, l (respectively). Print the two strings,
   then the difference in the first character at which they differ
*/
void generate_strings( int j, int k, int l )
{
  int s[256];
  int t[256];
  int cmp;
  int i;

  for( i = 0; i < j; i++ ) {
    s[i] = random_character();
    t[i] = s[i];
  }
  for( i = j; i < j+k; i++ ) {
    s[i] = random_character();
  }
  s[j+k] = '\0';
  for( i = j; i < j+l; i++ ) {
    t[i] = random_character();
  }
  t[j+l] = '\0';

  if( j+k > 0 ) {
    print_string( stdout,j+k,s );
  }
  else {
    printf(",");
  }
  printf("\n");
  if( j+l > 0 ) {
    print_string( stdout,j+l,t );
  }
  else {
    printf(",");
  }
  printf("\n");

  if( k == 0 && l == 0 ) {
    cmp = 0;
  }
  else if( k == 0 ) {
    cmp = -t[j];
  }
  else if( l == 0 ) {
    cmp = s[j];
  }
  else {
    i = j;
    while( i < j+k && s[j] == t[j] ) {
      j++; // check if additional chars happen to match
    }
    cmp = s[j] - t[j];
  }

  printf(",%d\n", cmp );
}

