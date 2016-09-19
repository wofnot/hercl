/** \file pole_cart.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../base/step.h"
#include "../base/eval.h"
#include "pole_cart.h"

//#define  HALT -1
//#define  X_COST  1

#define  BANG_BANG  1

//#define  ADD_NOISE  1

EnvState es1,es2,es3;

long noise_seed;

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
  double x0=0,x1=0,dx0=0,dx1=0,H0=0,H1=0,dH0=0,dH1=0;
  //`[x0|x1|x'0|x'1|H0|H1|H'0|H'1|xmin|xmax|Hmin|Hmax|Tmax|tau|mc|mp|l|g]
#ifdef ADD_NOISE
  noise_seed = random();
#endif
  es->running = TRUE;
  es->score.penalty = FALSE;

  fetch_input( in );
  scan_value(in,&x0);
  scan_value(in,&x1);
  scan_value(in,&dx0);
  scan_value(in,&dx1);
  scan_value(in,&H0);
  scan_value(in,&H1);
  scan_value(in,&dH0);
  scan_value(in,&dH1);
  scan_value(in,&es->xmin);
  scan_value(in,&es->xmax);
  scan_value(in,&es->Hmin);
  scan_value(in,&es->Hmax);
  scan_value(in,&es->Tmax);
  scan_value(in,&es->tau);
  scan_value(in,&es->mc);
  scan_value(in,&es->mp);
  scan_value(in,&es->l);
  scan_value(in,&es->g);

  // copy fixed parameters to es1,es2,es3
  es1 = *es;
  es2 = *es;
  es3 = *es;

  es->x  =  x0 + random_uniform()*(  x1 -  x0 );
  es->dx = dx0 + random_uniform()*( dx1 - dx0 );
  es->H  =  H0 + random_uniform()*(  H1 -  H0 );
  es->dH = dH0 + random_uniform()*( dH1 - dH0 );
  es->T = 0.0;

  write_value( out,es->x );
  write_value( out,es->dx );
  write_value( out,es->H );
  write_value( out,es->dH );
  output_message( out );

  return( TRUE );
}

/********************************************************//**

*/
void compute_ddxH( void *void_es, double F )
{
  EnvState *es = (EnvState *)void_es;
  double dH2sinH;
  dH2sinH = es->dH * es->dH * sin(es->H);
  es->ddH = (es->g * sin(es->H) + cos(es->H) *
                ( -F - es->mp * es->l*dH2sinH)/(es->mc + es->mp))
       /(es->l*(4.0/3.0 -(es->mp*cos(es->H)*cos(es->H))/(es->mc+es->mp)));
  es->ddx = (F + es->mp * es->l *(dH2sinH - es->ddH*cos(es->H)))/(es->mc + es->mp);
}

/********************************************************//**

*/
void env_set_cost( EnvState *es, double cost )
{
  es->score.cost = cost;
#ifdef X_COST
  if( es->x > 1.0 ) {
    es->score.cost += ( es->x - 1.0 )*( es->x - 1.0 );
  }
  else if( es->x < -1.0 ) {
    es->score.cost += ( es->x + 1.0 )*( es->x + 1.0 );
  }
#endif
  es->running = FALSE;
}

/********************************************************//**

*/
void  env_compute( EnvState *es, double F )
{
#ifdef ADD_NOISE
  long next_seed;
#endif
#ifdef BANG_BANG
  if( F >= 0.0 ) {
    F =  10.0;
  }
  else {
    F = -10.0;
  }
#else
  if( F > 10.0 ) {
    F =  10.0;
  }
  else if( F < -10.0 ) {
    F = -10.0;
  }
#endif
#ifdef ADD_NOISE
  next_seed = random();
  my_srandom( noise_seed );
  F += random_gaussian();
  noise_seed = random();
  my_srandom( next_seed );
#endif
  compute_ddxH( es,F );

  es1.x  = es->x  + 0.5*es->tau * es->dx;
  es1.dx = es->dx + 0.5*es->tau * es->ddx;
  es1.H  = es->H  + 0.5*es->tau * es->dH;
  es1.dH = es->dH + 0.5*es->tau * es->ddH;

  compute_ddxH(&es1,F );

  es2.x  = es->x  + 0.5*es->tau * es1.dx;
  es2.dx = es->dx + 0.5*es->tau * es1.ddx;
  es2.H  = es->H  + 0.5*es->tau * es1.dH;
  es2.dH = es->dH + 0.5*es->tau * es1.ddH;

  compute_ddxH(&es2,F );

  es3.x  = es->x  + es->tau * es2.dx;
  es3.dx = es->dx + es->tau * es2.ddx;
  es3.H  = es->H  + es->tau * es2.dH;
  es3.dH = es->dH + es->tau * es2.ddH;

  compute_ddxH(&es3,F );

  es->x  += es->tau*(es->dx + 2*es1.dx + 2*es2.dx + es3.dx )/6.0;
  es->dx += es->tau*(es->ddx+ 2*es1.ddx+ 2*es2.ddx+ es3.ddx)/6.0;
  es->H  += es->tau*(es->dH + 2*es1.dH + 2*es2.dH + es3.dH )/6.0;
  es->dH += es->tau*(es->ddH+ 2*es1.ddH+ 2*es2.ddH+ es3.ddH)/6.0;
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
  double F;
  if( !fetch_input(in)) {
    env_set_cost( es, es->Tmax - es->T );
    es->score.successful = FALSE; // agent not successful
    return( FALSE );
  }
  es->T += 0.5 * es->tau;
  if( !scan_value(in,&F)) {
    env_set_cost( es, es->Tmax - es->T );
    es->score.successful = FALSE; // agent not successful
    return( FALSE );
  }
  es->T += 0.5 * es->tau;

  env_compute( es,F );

  if( es->T >= es->Tmax ) {
    env_set_cost( es,0.0 );
    es->score.successful = TRUE; // agent is  successful
    return( FALSE );
  }
  if(   es->x < es->xmin || es->x > es->xmax
     || es->H < es->Hmin || es->H > es->Hmax ) {
    env_set_cost( es, es->Tmax - es->T );
    es->score.successful = FALSE; // agent not successful
    return( FALSE );
  }

  write_value( out,es->x );
  write_value( out,es->dx );
  write_value( out,es->H );
  write_value( out,es->dH );
  output_message( out );

  return( TRUE );
}

/********************************************************//**
   Return the score assigned to the agent by the environment
*/
Score *env_get_score( void *void_es )
{
  EnvState *es = (EnvState *)void_es;
  return(&(es->score));
}

