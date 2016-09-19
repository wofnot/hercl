/** chrstr.c */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "mystring.h"

/********************************************************//**
   Input: a character followed by a string
  Output: index of first occurrence of char in string
         (or empty message, if char doesn't occur)
*/
int main( int argc, char *argv[] )
{
  int s[256];
  int ch;
  int tests = NUM_TRAIN_TEST;
  int i,j,k,n;
  
  if( argc > 1 ){
    tests = atoi(argv[1]);
  }

  srandom((argc < 3 ? time(NULL) : atoi(argv[2])));

  printf("2\n"); // inputs_per_item

  for( n=1; n <= 2*tests; n++ ) {

    j = -2 + random_geometric( 4,32 );
    k =  j + random_geometric( 4,32 );

    for( i=0; i < k; i++ ) {
      s[i] = random_character();
    }
    if( j == -1 ) {
      ch = random_character(); // character not from string
    }
    else { // otherwise, choose a character in the string (uniformly)
      ch = s[j];
    }
    print_character( stdout,ch );
    printf("\n");

    if( k > 0 ) {
      print_string( stdout,k,s );
    }
    else {
      printf(",");
    }
    printf("\n");

    j=0; // make sure we find the first occurrence
    while( j < k && s[j] != ch ) {
      j++;
    }
    if( j == k ) { // ch does not occur
      printf(",\n"); // empty output
    }
    else {
      printf(",%d\n",j);
    }
    if( n == tests ) {
      printf("\\\n");
    }
  }

  return 0;
}
