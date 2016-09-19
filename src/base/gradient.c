/** \file gradient.c
*/

#include <stdio.h>
#include <math.h>

#include "step.h"
#include "gradient.h"

#define  EPS    0.01

/********************************************************//**
*/
unsigned int next_grad(
                       Agent *agt,
                       int op,
                       unsigned int arg0_grad,
                       unsigned int arg1_grad,
                       double value
                      )
{
  Gradient *result;
  agt->grad_index++;
  result = &(agt->grad[agt->grad_index]);
  result->op        = op;
  result->arg0_grad = arg0_grad;
  result->arg1_grad = arg1_grad;
  result->value     = value;
  result->delta     = 0.0;

  return( agt->grad_index );
}

/********************************************************//**

*/
void grad_step(
               int  op,
               Agent   *agt,
               Snode   *sf,
               Channel *out,
               Index    d,
               int  mem_index
              )
{
  int sp1;

  switch( op ) {

  case ADD: case MLT:
    sp1 = stack_adjust( agt, 1 );
    //if( isfinite(agt->stack[sp1])) {
    agt->stack_grad[agt->sp]=next_grad(agt,op,agt->stack_grad[agt->sp],
                              agt->stack_grad[sp1],agt->stack[agt->sp]);
    agt->stack_grad[sp1] = 0;
    break;

    // case AND: no action

  case ASN: case EXP: case LOG: case RCP: case RND: case SQT: case TNH:
    agt->stack_grad[agt->sp]=next_grad(agt,op,agt->stack_grad[agt->sp],
                                                 0,agt->stack[agt->sp]);
    break;

    // case BLN: rely on return_from_jump( agt );
    // case BRB: no action
    // case BRF: no action unless return_from_jump()

  case CPY:
    sp1 = stack_adjust( agt, -1 );
    agt->stack_grad[agt->sp] = agt->stack_grad[sp1];
    break;

  case DEC: case INC:
    if( sf->r < 0 ) {            // local
      sf->local_grad=next_grad(agt,op,sf->local_grad,0,sf->local_reg);
    }
    else {                       // global
      agt->reg_grad[sf->r] = next_grad(agt,op,agt->reg_grad[sf->r],0,
                                       agt->reg[sf->r]);
    }
    break;

    // case EQL: no action
    // case GRT: no action
    // case INC: same as DEC
    // case EXP: same as ASN

  case GET:
    if( sf->r < 0 ) {            // local
      agt->stack_grad[agt->sp] = sf->local_grad;
    }
    else {                       // global
      agt->stack_grad[agt->sp] = agt->reg_grad[sf->r];
    }
    break;

    // case INP: no action ??
    // case RND: same as ASN
    // case JSR: rely on push_stack_frame()

  case LOD:
    if( mem_index >= 0 ) {
      agt->stack_grad[agt->sp] = agt->mem_grad[mem_index];
    }
    else {
      agt->stack_grad[agt->sp] = 0;
    }
    break;

  case STO:
    sp1 = stack_adjust( agt, +1 );
    agt->mem_grad[mem_index] = agt->stack_grad[sp1];
    agt->stack_grad[sp1] = 0;
    break;

    // case LOG: same as ASN
    // case MLT: same as ADD ???

  case MOD: case PLR: case TRG:
    sp1 = stack_adjust( agt, -1 );
    agt->stack_grad[sp1]    =next_grad(agt,op,agt->stack_grad[agt->sp],
                              agt->stack_grad[sp1],agt->stack[agt->sp]);
    agt->stack_grad[agt->sp]=next_grad(agt,op,agt->stack_grad[agt->sp],
                              agt->stack_grad[sp1],agt->stack[agt->sp]);
    break;

    // case NOT: no action
    // case OR: no action
    // case OUT: rely on output_message ???
    // case PAY: no action
    // case PLR: same as MOD

  case POP:
    sp1 = stack_adjust( agt, +1 );
    agt->stack_grad[sp1] = 0;
    break;

  case NEG:
    if( sf->k > d+1 ) { // #- combination
      if( is_dot_or_digit(agt->cd->codon[d])) { //  value specified
        agt->stack_grad[agt->sp]=next_grad(agt,op,agt->cd->ival[d+1],
                                 0,agt->cd->fval[agt->cd->ival[d+1]]);
      }
      else {
        agt->stack_grad[agt->sp] = 0;
      }
    }
    else {  // stand-alone negation
      agt->stack_grad[agt->sp]=next_grad(agt,PSH,agt->stack_grad[agt->sp],
                                                    0,agt->stack[agt->sp]);
    }
    break;

  case PSH:
    if( sf->k > d+1 ) {
      agt->stack_grad[agt->sp]=next_grad(agt,op,agt->cd->ival[d+1],
                                0,agt->cd->fval[agt->cd->ival[d+1]]);
    }
    else {
      agt->stack_grad[agt->sp] = 0;
    }
    break;

  case PUT:
    sp1 = stack_adjust( agt, +1 );
    if( sf->r < 0 ) {            // local
      sf->local_grad = agt->stack_grad[sp1];
    }
    else {                       // global
      agt->reg_grad[sf->r] = agt->stack_grad[sp1];
    }
    agt->stack_grad[sp1] = 0;
    break;

    // case RCP: same as ASN

  case RAN:
    agt->stack_grad[agt->sp]=next_grad(agt,op,0,0,agt->stack[agt->sp]);
    break;

  case ROT:
    {
      int sp2;
      unsigned int grad_sp;
      sp1 = stack_adjust( agt, -1 );
      sp2 = stack_adjust( agt, -2 );
      grad_sp = agt->stack_grad[agt->sp];
      agt->stack_grad[agt->sp] = agt->stack_grad[sp1];
      agt->stack_grad[sp1]     = agt->stack_grad[sp2];
      agt->stack_grad[sp2]     = grad_sp;
    }
    break;

  case SCN:
    if( agt->bstack[agt->bp] == TRUE ) {
      agt->stack_grad[agt->sp]=next_grad(agt,op,0,0,agt->stack[agt->sp]);
    }
    break;

    // case SQT: same as ASN

  case SWP:
    {
      unsigned int grad_sp;
      sp1 = stack_adjust( agt, -1 );
      grad_sp = agt->stack_grad[agt->sp];
      agt->stack_grad[agt->sp] = agt->stack_grad[sp1];
      agt->stack_grad[sp1]     = grad_sp;
    }
    break;

    // case TNH: same as ASN
    // case TRG: same as MOD, PLR

  case WRT: // TRICKY BIT ???????
    sp1 = stack_adjust( agt, 1 );
    if( out != NULL &&( out->max_om < 0 || out->max_om > out->om )) {
      agt->write_grad[agt->wgi++] = agt->stack_grad[sp1];
    }
    agt->stack_grad[sp1] = 0;
    break;
  }
}

/********************************************************//**

*/
void grad_backprop( Agent *agt )
{
  Gradient *result;
  Gradient *arg0;
  Gradient *arg1;
  unsigned int gi;

  for( gi = agt->grad_index; gi > 0; gi-- ) {

    result = &agt->grad[gi];
    arg0   = &agt->grad[result->arg0_grad];
    arg1   = &agt->grad[result->arg1_grad];

    if( isfinite(result->value) && isfinite(result->delta)) {

      switch( result->op ) {

      case ADD:
        arg0->delta += result->delta;
        arg1->delta += result->delta;
        break;

        // case AND: no action

      case ASN:
        if( arg0->value < 1.0 && arg0->value > -1.0 ) {
          arg0->delta += result->delta*(EPS+1.0)/(EPS+cos(result->value));
        }
        break;

        // case BLN: no action
        // case BRB: no action
        // case BRF: no action
        // case CPY: no action

      case DEC: case INC:
        arg0->delta += result->delta;
        break;

        // case EQL: no action
        // case GRT: no action
        // case INC: same as DEC

      case EXP:
        arg0->delta += result->delta * result->value;
        break;

        // case GET: no action
        // case INP: no action
        // case RND: no action ????
        // case JSR: no action
        // case LOD: no action
        // case STO: no action

      case LOG:
        if( arg0->value > 0.0 ) {
          arg0->delta += result->delta/(EPS+arg0->value);
        }
        break;

      case MLT:
        if( isfinite(arg0->value) && isfinite(arg1->value)) {
          arg0->delta += result->delta * arg1->value;
          arg1->delta += result->delta * arg0->value;
        }
        break;

      case MOD: // ..(d)(r) -> ..(d/r)(d % r)
        /*             arg0     arg1     result
                     ----------------------------
          grid[gi-1] dividend  divisor  quotient
          grid[gi]   dividend  divisor  remainder
        */
        if( arg1->value != 0.0 ) {
          double dividend,divisor,quotient,remainder;
          dividend  = arg0->value;
          divisor   = arg1->value;
          quotient  =(result-1)->value;
          remainder = result->value;
          gi--;
          if( isfinite( dividend )) {
            if( isfinite( divisor )) {
              double phase  = remainder * 0.5*M_PI/divisor;
              double factor = cos(phase)*cos(phase);
              arg0->delta += result->delta *(1.0 - factor);
              if( quotient != 0.0 ) {
                arg1->delta -= result->delta *(1.0 - factor);
              }
              result--; // quotient
              if( fabs(divisor) > 1.0 ) {
                arg0->delta += result->delta*factor/ divisor;
                arg1->delta -= result->delta*factor/(divisor*divisor);
              }
              else {
                arg0->delta += result->delta*factor*mysign(divisor);
                arg1->delta -= result->delta*factor;
              }
            }
            else {
              arg0->delta += result->delta;
            }
          }
        }
        break;

      case NEG:
        arg0->delta -= result->delta;
        break;

        // case NOT: no action
        // case OR:  no action

        // case OUT: ????????????????
        // output_message( out );
        // break;

        // case PAY: no action

      case PSH:
        agt->delta[result->arg0_grad] += result->delta;
        break;

      case PLR: // ..(y)(x) -> ..(theta)(r)
        /*            arg0 arg1  result
                     ----------------
          grid[gi-1]    y    x   theta
          grid[gi]      y    x     r
        */
        if(isfinite(result->value) && result->value != 0.0){
          double x,y,r,theta,factor=1.0;
          x     = arg1->value;
          y     = arg0->value;
          r     = result->value;
          theta =(result-1)->value;
          arg0->delta += result->delta * y/r;
          arg1->delta += result->delta * x/r;
          gi--;
          result--;
          if( fabs(theta) > 0.5*M_PI ) {
            // attenuate delta near the boundary
            result->delta *= sin(theta)*sin(theta);
          }
          if( fabs(r) > 1.0 ) {
            // delta also attenuated near the origin
            factor = 1.0/(x*x+y*y);
          }
          arg0->delta += result->delta * x * factor;
          arg1->delta -= result->delta * y * factor;
        }
        break;

        // case POP: no action
        // case PSH: no action
        // case PUT: no action

      case RCP:
        // arg0->delta -= result->delta*result->value*result->value;
        arg0->delta -= result->delta/(EPS + arg0->value*arg0->value);
        break;

        // case RAN: no action
        // case ROT: no action
        // case SCN: no action

      case SQT:
        if( arg0->value > 0.0 ) {
          //arg0->delta += result->delta * 0.5/result->value;
          arg0->delta += result->delta/(EPS+2.0*result->value);
        }
        break;

        // case SWP: no action

      case TNH:
        arg0->delta += result->delta *
                (1.0 - result->value*result->value);
        break;

      case TRG: // ..(theta)(r) -> (y)(x)
        /*              arg0 arg1 result
                       -----------------
           grad[gi-1]  theta   r     y
           grad[gi]    theta   r     x
        */
        {
          double x,y,r,theta; //,factor=1.0;
          x = result->value;
          y =(result-1)->value;
          r = arg1->value;
          theta = arg0->value;
          //if( r > 1.0 ) {
          //  factor = (1.0 + EPS)/(1.0 + r*EPS);
          //}
          gi--;
          if( isfinite(r) && isfinite(theta)) {
            //arg0->delta -= result->delta * y * factor;
            arg0->delta -= result->delta * y;
            arg1->delta += result->delta * cos(theta);
            result--; // y
            //arg0->delta += result->delta * x * factor;
            arg0->delta += result->delta * x;
            arg1->delta += result->delta * sin(theta);
          }
        }
        break;

      case WRT:
        /*   NEEDS WORK ??????
    if( out != NULL && out->om < out->max_om
                    && out->wi < out->max_i ) {
      out->val[out->wi++] = agt->stack[agt->sp];
    }
    agt->stack[agt->sp] = 0.0;
    agt->sp = stack_adjust( agt, -1 );
        */
        break;
      }
    }
  }
}

/*
double regularize( double x )
{
  if( fabs(x) <= 1.0 ) {
    return( x );
  }
  else {
    return(( EPS + 1.0)/( EPS + 1.0/x ));
  }
}
*/
