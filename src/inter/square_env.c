/** \file square_env.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../base/step.h"
#include "../base/eval.h"
#include "square_env.h"

#define TOLERANCE 0.001

/********************************************************//**
   Allocate a new EnvState
*/
void * new_env_state()
{
  EnvState *es = (EnvState *)malloc(sizeof(EnvState));
  check_null( es );
  memset(es,0,sizeof(EnvState));
  return((void *)es);
}

/********************************************************//**
   Free state occupied by env_state
*/
void free_env_state( void *void_es )
{
  EnvState *es = (EnvState *)void_es;
  free( es );
}

/********************************************************//**
   Initialize the environment
*/
Boolean env_reset(
                  void    *void_es,
                  Channel *in,
                  Channel *out,
                  int r
                 )
{
  EnvState *es = (EnvState *)void_es;
  es->score.penalty = FALSE;
  
  double min_input;
  double max_input;
  
  ;
  if(! (  fetch_input( in ) 
        && scan_value( in,&min_input )
        && scan_value( in,&max_input ))) {

    // if no cfg
    min_input = -32.0;
    max_input = 32.0;
  }
  
  // input to agent (A number that should be squared)
  double squareThis = (random_uniform()*(max_input-min_input)+min_input);
  // expected output
  es->expected = squareThis * squareThis;
  
  //printf("input %lf, expected %lf\n", squareThis, es->expected);

  write_value( out, squareThis );
  output_message( out );

  return( TRUE );
}

/********************************************************//**
   Continue running environment.
   Return TRUE if the trial is still running, FALSE otherwise
*/
Boolean env_continue(
                     void    *void_es,
                     Channel *in,
                     Channel *out
                    )
{
  EnvState *es = (EnvState *)void_es;
  double actual = -1.0;
  
  es->score.successful = FALSE;
  
  // if agent fails to output anything, penalise
  if( !fetch_input(in)) {
    es->score.cost = 10000.0;
    es->score.penalty = 1.0;
    es->score.reject = FALSE;    
  
  } else {
  // if agent output has no items, penalise (slightly less)
    if( !scan_value(in,&actual)) {
      es->score.cost = 5000.0;
      es->score.penalty = 1.0;
      es->score.reject = FALSE;
  
    // Otherwise, calculate the score
    } else {
      double raw_cost = d_abs( es->expected - actual ); 
      if(raw_cost < TOLERANCE){
        es->score.cost = 0.0;
        es->score.successful = TRUE;
      } else {
        es->score.cost = raw_cost;
      }
      es->score.penalty = 0.0;
      es->score.reject = FALSE;
    }
  }
  // The trial is always finished after one round.
  //printf("expected: %f \n", es->expected);
  //printf("actual: %f \n", actual);
  return( FALSE );
}

/********************************************************//**
   Return the score assigned to the agent by the environment
*/
Score *env_get_score( void *void_es )
{
  EnvState *es = (EnvState *)void_es;
  // printf("-----score returned-----\n");
  //printf("cost: %f \n", es->score.cost);
  //printf("penalty: %f \n", es->score.penalty);
  // printf("reject: %d \n", es->score.reject);
  return(&(es->score));
}

/********************************************************//**
   return the absolute value of a double
*/
double d_abs( double n )
{
  if(n < 0.0){
    n *= -1.0;
  }
  return n;
}



