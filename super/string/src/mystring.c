/** mystring.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#define  MAX_INT_PLUS_ONE  2147483648.0

#define  TRUE    1
#define  FALSE   0

int random_gaussian_phase=0;

/********************************************************//**
   Return a random number chosen uniformly between 0 and 1
*/
double random_uniform()
{
  int n = random();
  return( n/MAX_INT_PLUS_ONE );
}

/********************************************************//**
   Return random variable between 1 and max (inclusive)
   from a geometric distribution with the specified mean
*/
int random_geometric( double mean, int max )
{
  double x;
  int n;
  if( mean <= 0.0 ) {
    return 0;
  }
  do {
    x = 1.0 - random_uniform();
    n = 1+(int)floor(log(x)/log(1.0-1.0/mean));
  } while( n > max );

  return( n );
}

int int_cauchy( int min, int max )
{
  double x,z;
  int n;
  do {
    x = random_uniform();
    if( min == 0 ) {
      x = 0.4 + 0.6 * x;
    }
    else {
      x = 0.5 + 0.5 * x;
    }
    z = tan( x * 3.14159265359/2.0);
    n = (int)floor(z);
  } while( n < min || n > max );

  return( n );
}

/********************************************************//**
   Return random variable from a Gaussian distribution
*/
double random_gaussian()
{
  static double V1, V2, S;
  double X;
  if(random_gaussian_phase == 0) {
    do {
      double U1 = (double) random() / MAX_INT_PLUS_ONE;
      double U2 = (double) random() / MAX_INT_PLUS_ONE;
      V1 = 2.0 * U1 - 1.0;
      V2 = 2.0 * U2 - 1.0;
      S = V1 * V1 + V2 * V2;
    } while( S >= 1.0 || S == 0.0 );

    X = V1 * sqrt(-2.0 * log(S) / S);
  }
  else {
    X = V2 * sqrt(-2.0 * log(S) / S);
  }
  random_gaussian_phase = 1 - random_gaussian_phase;

  return X;
}

/********************************************************//**
   Return a random (printable) character
*/
int random_character()
{
  int ch;
  do {
    ch = random()% 256;
  } while( !isprint( ch ));
  return( ch );
}

/********************************************************//**
   Print string in a form readable by scan_inputs()
*/
void print_string( FILE *fp, int length, int s[] )
{
  int preceded_by_char = TRUE;
  int ch;
  int k;

  if( length == -1 ) {
    fprintf(fp,"\\");
  }
  else if( length == 0 ) {
    fprintf(fp,",");
  }
  for( k=0; k < length; k++ ) {
    ch = s[k];
    if( ch >  0 && ch < 256 && isprint( ch )) {
      if( !preceded_by_char ) {
        putc(' ',fp);
      }
      preceded_by_char = TRUE;
      if( ch == ',' ) {
        fprintf(fp,"\\,");
      }
      else if( ch == '\\' ) {
        fprintf(fp,"\\\\");
      }
      else {
        putc( ch, fp );
      }
    }
    else {
      fprintf(fp,",%d",ch);
      preceded_by_char = FALSE;
    }
  }
}

/********************************************************//**
   Print character in a form readable by scan_inputs()
*/
void print_character( FILE *fp, int ch )
{
  if( ch > 0 && ch < 256 && isprint( ch )) {
    if( ch == '\\' || ch == ',' ) {
      putc( '\\', fp );
    }
    putc( ch, fp );
  }
  else {
    fprintf( fp,",%d",ch );
  }
}

double random_cauchy( int min )
{
  double x = random_uniform();
  if( min == 1 ) {
    x = 0.5 *( x + 1.0 );
  }
  return( tan( x * 3.14159265359/2.0));
}

double signed_cauchy()
{
  double x = random_uniform();
  return( tan( 3.14159265359 *( x - 0.5 )));
}

/********************************************************//**
   Return random variable from a pseudo-Cauchy distribution
*/
double pseudo_cauchy( double median )
{
  double x = 1.0 - random_uniform();
  return( median/x - median );
}

/********************************************************//**
   Divide data into training and test sets
*/
void divide_data(
                 int train_items,
                 int test_items,
                 int item_in[],
                 int item_out[]
                )
{
  int j,k;
  j = train_items;
  k = test_items;
  while( j+k > 0 ) {
    if( random()%(j+k) < j ) {
      item_out[j] = item_in[j+k];
      j--;
    }
    else {
      item_out[train_items + k] = item_in[j+k];
      k--;
    }
  }
}
