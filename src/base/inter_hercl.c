/** \file env_hercl.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "step.h"
#include "eval.h"
#include "inter_hercl.h"

#define  BANG_BANG  1

Code *env_code = NULL;

int max_env_step;

long noise_seed;

/********************************************************//**
   Set the (global) variable env_code to the provided code 
*/
void set_env_code( Code *cd )
{
  env_code = cd;
}

/********************************************************//**
   Assume env_code has already been set by set_env_code()
   Create a new candidate, with the pre-defined env_code
   and return it.
*/
void * new_env_state()
{
  Candidate *task = new_candidate();
  task->agt->cd = env_code;
  return((void *)task);
}

/********************************************************//**
   Free state occupied by env_state
*/
void free_env_state( void *void_es )
{
  Candidate *task = new_candidate();
  free_candidate( task );
}

/********************************************************//**
   Initialize the task environment, and let it run until<br>
   (a) it produces its first output (return TRUE), or<br>
   (b) it exceeds the maximum number of steps, or halts
   without producing output (update task->score.nsteps,
   set task->score.reject to TRUE and return FALSE).
*/
Boolean env_reset(
                  void *void_env_state,
                  Channel *in,
                  Channel *out,
                  int r
                 )
{
  Candidate *task = (Candidate *)void_env_state;
  Agent *env = task->agt;
  int op = NONE;

  task->score.num_trials++;
  max_env_step = 10 * ENV_STEPS_PER_MOVE;
  //task_cfg = in; // save for later use

  reset_agent( env );

  while(   op != OUT && env->running
        && env->step <= max_env_step ) {
    op = step( env,in,out );
  }
  if( !env->running || env->step > max_env_step ) {
    env->running = FALSE;
    task->score.reject = TRUE;
    task->score.num_steps += env->step;
    //printf("Environment failed to produce output\n");
    return FALSE; // environment has failed
  }
  else {
    return TRUE;
  }
}

/********************************************************//**
   Continue running the task environment until it<br>
   (a) produces output (return TRUE)<br>
   (b) halts (update task->score.nsteps ad return FALSE), or <br>
   (b) exceeds the maximum number of steps
      (update task->score.nsteps,
       set task->score.reject=TRUE and return FALSE).
*/
Boolean env_continue(
                     void    *void_env_state,
                     Channel *in,
                     Channel *out
                    )
{
  Candidate *task = (Candidate *)void_env_state;
  Agent      *env = task->agt;
  int op;

  max_env_step += ENV_STEPS_PER_MOVE; 
  op = NONE;
  while( op != OUT && env->running
                   && env->step <= max_env_step ) {
    op = step( env,in,out );
  }
  if( env->step > max_env_step ) {
    task->score.reject = TRUE;
    env->running = FALSE;
  }
  if( !env->running ) {
    task->score.num_steps += env->step;
    if( env->stack_items > 0 ) {
      task->score.cost = top_of_stack( env );
    }
    task->score.successful = (    env->bstack[env->bp]
                              && task->score.cost <= 0.0 );
    return( FALSE ); // interaction has completed
  }
  else {
    return( TRUE );  // interaction still running
  }
}

/********************************************************//**
  Return the score for this task
*/
Score *env_get_score( void *void_es )
{
  Candidate *task = (Candidate *)void_es;
  return(&(task->score));
}

/********************************************************//**
  Return TRUE if environment is still running, FALSE otherwise
int env_running( void *void_es )
{
  Candidate *task = (Candidate *)void_es;
  return( task->agt->running );
}
*/
