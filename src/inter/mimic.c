/** \file mimic.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../base/step.h"
#include "../base/eval.h"
#include "../base/scan_print.h"
#include "mimic.h"

#define TOLERANCE 0.001

#define MAX_INPUT 25
#define P         2.43290200817664E+18 // penalty co-efficient (25!)

// HERCL program for the agent to try and mimic
#define PROGRAM "./lib/hrc/factorial.hrc"

//#define DEBUG 1

/********************************************************//**
   Allocate a new EnvState
*/
void * new_env_state()
{
  EnvState *es = (EnvState *)malloc(sizeof(EnvState));
  check_null( es );
  
  memset(es,0,sizeof(EnvState));
  
  // Create an agent containing the code to mimic
  es->agent = new_agent();
  FILE *fp = fopen_check( PROGRAM ,(char *)"r" );

  es->agent->cd = scan_code( fp );
  
  fclose(fp);
  
  // Allocate channels for this agent
  es->chan = new_channel();
  
  return((void *)es);
}

/********************************************************//**
   Free state occupied by env_state
*/
void free_env_state( void *void_es )
{
  EnvState *es = (EnvState *)void_es;  
  
  // Free agent
  free_code( es->agent->cd );
  compress_agent( es->agent );
  free( es->agent );
  
  // Free channels
  free_channel(es->chan);
  
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
  es->score.reject = FALSE;
  
  double input = random_uniform() * MAX_INPUT + 1; //min of 1
  #ifdef DEBUG
  printf("Input: %f\n", input);
  #endif
  
  // run the agent to mimic  
  reset_agent( es->agent );
  
  // Give agent-to-mimic input
  clear_message(es->chan, 1);
  write_value(es->chan, input);
  output_message( es->chan );
  cue_message(es->chan, es->chan->om);
  
  int op = NONE;
  while( op != OUT && es->agent->running ){ 
    //Note that this will loop indefinitely if agent provided does not output/terminate
    op = step( es->agent, es->chan, es->chan );
  }
  if(op == OUT){
    cue_message(es->chan, es->chan->om);
  }
  
  //get response from agent-to-mimic and set target to this value
  if(!fetch_input(es->chan)){
    printf("Error: agent-to-mimic failed to output message\n");
  }
  if(!scan_value(es->chan,&(es->target))){
    printf("Error: agent-to-mimic failed to output item\n");
  }
  
  #ifdef DEBUG
  printf("Target: %f\n", es->target);
  #endif
  
  // Give input to evolving agent
  write_value(out, input);
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
  
  es->score.successful = FALSE;
  
  double fguess;
  
  // if agent fails to output anything, penalise and stop testing
  if( !fetch_input(in) ) {
    es->score.cost = 20.0*P;
    es->score.penalty = 3.0;
    es->score.reject = FALSE;
    return( FALSE ); 
  }
  // if agent output has no items, penalise (slightly less) and stop
  if( !scan_value(in,&fguess)) {
    es->score.cost = 10.0*P;
    es->score.penalty = 2.0;
    es->score.reject = FALSE;
    return( FALSE );
  }  
  // if agent has two or more items in output, penalise (too many items)
  if( !input_exhausted( in ) ) {
    es->score.cost = 5.0*P;
    es->score.penalty = 2.0;
    es->score.reject = FALSE;
    return( FALSE );
  }
  
  es->score.cost = fabs(fguess - es->target)/ es->target; //cost by percentage error
  es->score.penalty = 0.0;
  es->score.reject = FALSE;
  
  if(es->score.cost < TOLERANCE){
    es->score.successful = TRUE;
  }
  
  return(FALSE); 
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

