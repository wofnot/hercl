/** \file double_pole.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../base/step.h"
#include "../base/eval.h"
#include "double_pole.h"

//#define  HALT -1
//#define  ADD_NOISE  1

#define  BANG_BANG  1

#define  BOTH_POLES_VARY  1

//#define  X_COST  1

EnvState es1,es2,es3;

long noise_seed;

extern long eval_seed;

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
                  void *void_es,
                  Channel *in,
                  Channel *out,
                  int r
                 )
{
  EnvState *es = (EnvState *)void_es;
  double x_lo=0,x_hi=0,dx_lo=0,dx_hi=0,
         H_lo=0,H_hi=0,dH_lo=0,dH_hi=0;
  //`[x0|x1|x'0|x'1|H0|H1|H'0|H'1|xmin|xmax|Hmin|Hmax|Tmax|tau|mc|mp|l|g]
#ifdef ADD_NOISE
  noise_seed = random();
#endif
  es->running = TRUE;
  es->score.penalty = FALSE;

  fetch_input( in );
  scan_value(in,&x_lo);
  scan_value(in,&x_hi);
  scan_value(in,&dx_lo);
  scan_value(in,&dx_hi);
  scan_value(in,&H_lo);
  scan_value(in,&H_hi);
  scan_value(in,&dH_lo);
  scan_value(in,&dH_hi);

  scan_value(in,&es->xmin);
  scan_value(in,&es->xmax);
  scan_value(in,&es->Hmin);
  scan_value(in,&es->Hmax);
  scan_value(in,&es->Tmax);
  scan_value(in,&es->tau);
  scan_value(in,&es->mc);
  scan_value(in,&es->m1);
  scan_value(in,&es->m2);
  scan_value(in,&es->l1);
  scan_value(in,&es->l2);
  scan_value(in,&es->g);

  // copy fixed parameters to es1,es2,es3
  es1 = *es;
  es2 = *es;
  es3 = *es;

  if( eval_seed == 0 ) {
    es->x   = 0.0;
    es->dx  = 0.0;
    es->H1  = 3.14159/180.0;
    es->dH1 = 0.0;
    es->H2  = 0.0;
    es->dH2 = 0.0;
  }
  else {
    es->x   =  x_lo + random_uniform()*(  x_hi -  x_lo );
    es->dx  = dx_lo + random_uniform()*( dx_hi - dx_lo );
    es->H1  =  H_lo + random_uniform()*(  H_hi -  H_lo );
    es->dH1 = dH_lo + random_uniform()*( dH_hi - dH_lo );
#ifdef BOTH_POLES_VARY
    es->H2  =  H_lo + random_uniform()*(  H_hi -  H_lo );
    es->dH2 = dH_lo + random_uniform()*( dH_hi - dH_lo );
#else
    es->H2  = 0.0;
    es->dH2 = 0.0;
#endif
  }
  es->T   = 0.0;

  write_value( out,es->x );
  write_value( out,es->dx );
  write_value( out,es->H1 );
  write_value( out,es->dH1 );
  write_value( out,es->H2 );
  write_value( out,es->dH2 );
  output_message( out );

  return( TRUE );
}

/********************************************************//**

*/
void compute_ddxH( EnvState *es, double F )
{
  double dH2sinH;
  //double ddx,ddH1,ddH2;
  double F1,F2,m1,m2;
  /*
  dH2sinH = es->dH * es->dH * sin(es->H);
  ddH = (es->g * sin(es->H) + cos(es->H) *
                ( -F - es->mp * es->l*dH2sinH)/(es->mc + es->mp))
       /(es->l*(4.0/3.0 -(es->mp*cos(es->H)*cos(es->H))/(es->mc+es->mp)));
  ddx = (F + es->mp * es->l *(dH2sinH - ddH*cos(es->H)))/(es->mc + es->mp);
  */

  dH2sinH = es->dH1 * es->dH1 * sin(es->H1);
  F1 = es->m1*es->l1*dH2sinH - 0.75*es->m1*cos(es->H1)*es->g*sin(es->H1);
  dH2sinH = es->dH2 * es->dH2 * sin(es->H2);
  F2 = es->m2*es->l2*dH2sinH - 0.75*es->m2*cos(es->H2)*es->g*sin(es->H2);

  m1 = es->m1 *( 1.0 - 0.75 * cos(es->H1) * cos(es->H1));
  m2 = es->m2 *( 1.0 - 0.75 * cos(es->H2) * cos(es->H2));

  es->ddx = ( F + F1 + F2 )/( es->mc + m1 + m2 );

  es->ddH1 = -0.75 *( es->ddx * cos(es->H1) - es->g * sin(es->H1))/ es->l1;
  es->ddH2 = -0.75 *( es->ddx * cos(es->H2) - es->g * sin(es->H2))/ es->l2;

  //printf("%1.8f %1.8f %1.8f | %1.8f %1.8f %1.8f\n",
  //       ddx,es->ddx,es->ddx-ddx,ddH,es->ddH,es->ddH-ddH);
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
   Continue running environment.
   Return TRUE if the trial is still running, FALSE otherwise.
*/
Boolean env_continue(
                     void    *void_es,
                     Channel *in,
                     Channel *out
                    )
{
  EnvState *es = (EnvState *)void_es;
#ifdef ADD_NOISE
  long next_seed;
#endif
  double F;
  if( !fetch_input(in)) {
    env_set_cost( es, es->Tmax - es->T );
    es->score.successful = FALSE;
    return( FALSE );
  }
  es->T += 0.5 * es->tau;
  if( !scan_value(in,&F)) {
    env_set_cost( es, es->Tmax - es->T );
    es->score.successful = FALSE;
    return( FALSE );
  }
  es->T += 0.5 * es->tau;
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

  es1.x   = es->x   + 0.5*es->tau * es->dx;
  es1.dx  = es->dx  + 0.5*es->tau * es->ddx;
  es1.H1  = es->H1  + 0.5*es->tau * es->dH1;
  es1.H2  = es->H2  + 0.5*es->tau * es->dH2;
  es1.dH1 = es->dH1 + 0.5*es->tau * es->ddH1;
  es1.dH2 = es->dH2 + 0.5*es->tau * es->ddH2;

  compute_ddxH(&es1,F );

  es2.x   = es->x   + 0.5*es->tau * es1.dx;
  es2.dx  = es->dx  + 0.5*es->tau * es1.ddx;
  es2.H1  = es->H1  + 0.5*es->tau * es1.dH1;
  es2.H2  = es->H2  + 0.5*es->tau * es1.dH2;
  es2.dH1 = es->dH1 + 0.5*es->tau * es1.ddH1;
  es2.dH2 = es->dH2 + 0.5*es->tau * es1.ddH2;

  compute_ddxH(&es2,F );

  es3.x   = es->x   + es->tau * es2.dx;
  es3.dx  = es->dx  + es->tau * es2.ddx;
  es3.H1  = es->H1  + es->tau * es2.dH1;
  es3.H2  = es->H2  + es->tau * es2.dH2;
  es3.dH1 = es->dH1 + es->tau * es2.ddH1;
  es3.dH2 = es->dH2 + es->tau * es2.ddH2;

  compute_ddxH(&es3,F );

  es->x   += es->tau*(es->dx  + 2*es1.dx  + 2*es2.dx  + es3.dx  )/6.0;
  es->dx  += es->tau*(es->ddx + 2*es1.ddx + 2*es2.ddx + es3.ddx )/6.0;
  es->H1  += es->tau*(es->dH1 + 2*es1.dH1 + 2*es2.dH1 + es3.dH1 )/6.0;
  es->H2  += es->tau*(es->dH2 + 2*es1.dH2 + 2*es2.dH2 + es3.dH2 )/6.0;
  es->dH1 += es->tau*(es->ddH1+ 2*es1.ddH1+ 2*es2.ddH1+ es3.ddH1)/6.0;
  es->dH2 += es->tau*(es->ddH2+ 2*es1.ddH2+ 2*es2.ddH2+ es3.ddH2)/6.0;

  /*
  es->x  += es->tau * es->dx;
  es->dx += es->tau * es->ddx;
  es->H  += es->tau * es->dH;
  es->dH += es->tau * es->ddH;
  */

  //es->T += 0.5 * es->tau;
  if( es->T >= es->Tmax ) {
    env_set_cost( es,0.0 );
    es->score.successful = TRUE;
    return( FALSE );
  }
  if(   es->x  < es->xmin || es->x  > es->xmax
     || es->H1 < es->Hmin || es->H1 > es->Hmax
     || es->H2 < es->Hmin || es->H2 > es->Hmax ) {
    env_set_cost( es, es->Tmax - es->T );
    es->score.successful = FALSE;
    return( FALSE );
  }

  write_value( out,es->x );
  write_value( out,es->dx );
  write_value( out,es->H1 );
  write_value( out,es->dH1 );
  write_value( out,es->H2 );
  write_value( out,es->dH2 );
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

