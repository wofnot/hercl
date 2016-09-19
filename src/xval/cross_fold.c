/** \file cross_fold.c
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

#define  PRINT_TEST  1

void        usage(char *argv_0);
void    scan_text(FILE *fp,int num_lines,char *text[]);
int  * scan_folds(FILE *fp, int num_items);
void extract_name(char name[], char*argv_2);
void  print_folds(int num_items,char *text[],char  name[],
                  int num_fold,int fold[]);

/********************************************************/
int main( int argc, char *argv[] )
{
  FILE  *fp;
  Super *sup;
  int   *fold;
  char **text;
  char name[256];
  int  num_fold;

  if( argc < 4 ) {
    usage( argv[0] );
  }

  num_fold = atoi( argv[1] );

  fp = fopen_check( argv[2],"r" );
  sup = scan_super( fp );
  fclose( fp );

  fp = fopen_check( argv[3],"r" );
  fold = scan_folds( fp,sup->total_items );
  fclose( fp );

  text = (char **)malloc(2*sup->total_items*sizeof(char *));
  fp = fopen_check( argv[2],"r" );
  scan_text( fp, 2*sup->total_items, text );
  fclose( fp );

  extract_name( name,argv[2] );
  print_folds( sup->total_items,text,name,num_fold,fold );

  return 0;
}

/********************************************************//**
   Print a usage message and exit.
*/
void usage( char *argv_0 )
{
  fprintf( stderr,"Usage: %s <num> <super.in> <task.fold>\n",argv_0);
  exit(1);
}

/********************************************************//**
   scan the specified number of lines from fp into text[]
*/
void scan_text( FILE *fp, int num_lines, char *text[] )
{
  char line[1024];
  int n;
  for( n=0; n < num_lines; n++ ) {
    fgets( line,1024,fp );
    text[n]=(char *)malloc((strlen(line)+1)*sizeof(char));
    strncpy( text[n],line,strlen(line)+1);
  }
}

/********************************************************//**
   scan the fold information for num_items, from fp
*/
int * scan_folds( FILE *fp, int num_items )
{
  int *fold;
  int a,r,s;

  fold=(int *)malloc((num_items+1)*sizeof(int));
  check_null(fold);

  for( r=1; r <= num_items; r++ ) {
    fscanf( fp,"fold[%d] = %d\n", &s, &a );
    if( r != s ) {
      printf("HEY!\n");
      exit( 1 );
    }
    else {
      fold[r] = a+1;
    }
  }
  return( fold );
}

/********************************************************//**
   extract the name of the dataset from argv_2,
   by stripping the path from beginning and .in from the end,
   and store the result to name[]
*/
void extract_name( char name[], char*argv_2 )
{
  int i=0,j,k=0;

  while( argv_2[k] != '\0' ) {
    k++;
  }
  j = k;
  while( j > 0 && argv_2[j-1] != '/' ) {
    j--;
  }
  while ( j < k ) {
    name[i++] = argv_2[j++];
  }
  name[i] = '\0';
}

/********************************************************//**
   print the training and test items for each fold into
   an appropriately named file.
*/
void print_folds(
                 int num_items,
                 char *text[],
                 char  name[],
                 int num_fold,
                 int fold[]
                )
{
  FILE *fp;
  int i=0,k=0,a,n;

  while( name[k] != '\0' ) {
    k++;
  }
  i = k;
  while( i > 0 && name[i] != '.' ) {
    i--;
  }
  if( i == 0 ) {
    i = k;
  }
  while( k >= i ) {
    name[k+1] = name[k];
    k--;
  }

  for( a=1; a <= num_fold; a++ ) {
    name[i] = '0'+a-1;
    fp = fopen( name,"w" );
    for( n=1; n <= num_items; n++ ) {
      if( fold[n] != a ) {
        fprintf( fp,"%s",text[2*n-2]);
        fprintf( fp,"%s",text[2*n-1]);
      }
    }
    fprintf( fp,"\\\n");
#ifdef PRINT_TEST
    for( n=1; n <= num_items; n++ ) {
      if( fold[n] == a ) {
        fprintf( fp,"%s",text[2*n-2]);
        fprintf( fp,"%s",text[2*n-1]);
      }
    }
#endif
    fclose( fp );
  }
}
