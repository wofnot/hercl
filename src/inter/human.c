/** \file human.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../base/step.h"
#include "../base/eval.h"
#include "../base/scan_print.h"

#include "human.h"

int user_message(Channel *out);

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
  es->score.cost = 0;
  es->score.reject = FALSE;
  es->score.successful = TRUE;
  
  user_message(out);

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
  //EnvState *es = (EnvState *)void_es;
  double actual = -1.0;
  
  // if agent fails to output anything, penalise
  if( !fetch_input(in)) {
    printf("(no input)\n"); 
  } else {
    printf("input: ");
    while( scan_value(in,&actual)) {
      printf("%0.2lf ", actual);
    }
    printf("\n");
  }
  if (!out) {
    printf("Agent finished\n");
    return FALSE;
  }
  if (user_message(out)) {
    return TRUE; //keep going
  } else {
    return( FALSE ); // stop
  }
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
   Take single line from stdin and output it as a message to agent
   Return 0 and output nothing if EOF is received
*/
int user_message(Channel *out) {
  printf(" > ");
  char line[256];
  gets(line);
  if (feof(stdin)) {
    return 0;
  }
  scan_next_input(out, line);
  out->om--; // scan_next_input pushes om forward 1 place too many
  
  output_message( out );
  print_channel(out, stdout);
  return 1;
}


