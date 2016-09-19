/** \file gen_fold.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../base/step.h"
#include "../base/scan_print.h"
#include "../base/cross_mutate.h"
#include "../base/eval.h"
#include "../base/super.h"

void           usage(char *argv_0);
int        next_fold(int num_fold);
int * generate_folds(Super *sup,int num_fold);


/********************************************************//***/
int main( int argc, char *argv[] )
{
  FILE  *fp;
  Super *sup;
  int   *fold;
  long seed=0;
  int  num_fold;
  int r;

  if( argc < 3 ) {
    usage( argv[0] );
  }
  if( seed == 0 ) {
    seed = time( NULL );
  }
  //printf("seed=%ld\n",seed);
  srandom(seed);

  num_fold = atoi( argv[1] );

  fp = fopen_check( argv[2],"r" );
  sup = scan_super( fp );
  fclose( fp );

  fold = generate_folds( sup,num_fold );

  for( r=1; r <= sup->total_items; r++ ) {
    printf("fold[%d] = %d\n",r,fold[r]-1);
  }

  return 0;
}

/********************************************************//**
   Print a usage message and exit.
*/
void usage( char *argv_0 )
{
  fprintf( stderr,"Usage: %s <num> <super.in>\n",argv_0);
  exit(1);
}

/********************************************************//**
   
*/
int * generate_folds( Super *sup, int num_fold )
{
  Channel *target;
  int     *fold;
  int     *item;
  int same_target;
  int num_assigned=0,num_same;
  int i,j,k,r,s;

  target = sup->target;

  fold=(int *)malloc((sup->total_items+1)*sizeof(int));
  check_null(fold);
  memset(fold,0,(sup->total_items+1)*sizeof(int));
  item=(int *)malloc((sup->total_items+1)*sizeof(int));
  check_null(item);

  while( num_assigned < sup->total_items ) {
    do {
      s = 1 + random()% sup->total_items;
    } while( fold[s] != 0 );
    j=0;
    for( r=1; r <= sup->total_items; r++ ) {
      if( target->length[r] == target->length[s] ) {
        same_target = TRUE;
        for( i=0; i < target->length[r]; i++ ) {
          if(   target->val[target->index[r]+i]
             != target->val[target->index[s]+i] ) {
            same_target = FALSE;
          }
        }
        if( same_target ) {
          item[j++] = r;
        }
      }
    }
    num_same = j;
    k=0;
    while( k < num_same ) {
      j = random()% num_same;
      if( fold[item[j]] == 0 ) {
        fold[item[j]] = next_fold(num_fold);
        //printf("fold[%2d] = %2d\n",item[j],fold[item[j]]);
        k++;
      }
    }
    num_assigned += num_same;
  }

  return( fold );
}

/********************************************************//**
   
*/
int next_fold( int num_fold )
{
  static int this_fold=0;

  this_fold++;
  if( this_fold > num_fold ) {
    this_fold = 1;
  }
  return( this_fold );
}
