/** \file gradient.h
*/

void grad_step( int op,Agent *agt,Snode *sf,Channel *out,
                Index d,int mem_index );
void grad_backprop( Agent *agt );

