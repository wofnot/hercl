/** \file step.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include "step.h"
#include "gradient.h"

int     get_mem_index(Floating raw_index,int mem_size);
void push_stack_frame(Agent *agt,int c);
void    push_to_stack(Agent *agt,double value);
void return_from_jump(Agent *agt);
Floating     *set_reg(Agent *agt,Snode *sf,int d);


const int power2[17] = {
  1,2,4,8,16,32,64,128,256,512,1024,
  2048,4096,8192,16384,32768,65536};

const char phen[NUM_CODONS] = {
  ADD,AND,ASN,BLN,BRB,BRF,CPY,DEC,EQL,EXP,
  GET,GRT,INC,INP,JSR,LOD,LOG,MLT,MOD,NEG,
  NOT,OR, OUT,PAY,PLR,POP,PSH,PUT,RAN,RCP,
  RND,ROT,SCN,SQT,STO,SWP,TNH,TRG,WRT
};

int is_codon[256];
int op_uses_digits[256];
int is_allowed[256];

// these state variables are reset by init_mutation_state(),
// init_evaluation_state(), restore_mutation_state()
int mutation_state[32];
int evaluation_state[32];
int random_bit_int=0;
int random_gaussian_phase=0;
int mutation_bit_int=0;
int mutation_gaussian_phase=0;

Global_Param globals = {  // global variables
  FALSE,               // verbose
  FALSE,               // multi_process
  FALSE,               // terminate_early
  FALSE,               // no_share
   10,                 // base
    1.0,               // tune_magnitude
    1.0,               // interpolate_rate
    0.0,               // learning_rate
 99900000000           // max_eval
};

/********************************************************//**
   Initialize is_codon[],is_allowed[]
*/
void init_codons()
{
  int i,op;//,type;

  for( i=0; i < 256; i++ ) {
    is_codon[i]       = FALSE;
    is_allowed[i]     = FALSE;
    op_uses_digits[i] = FALSE;
  }
  for( i=0; i < NUM_CODONS; i++ ) {
    is_codon[(int)phen[i]] = TRUE;
    is_allowed[(int)phen[i]] = TRUE;
  }
  op_uses_digits[(int)GET] = TRUE; op_uses_digits[(int)PUT] = TRUE;
  op_uses_digits[(int)LOD] = TRUE; op_uses_digits[(int)STO] = TRUE;
  op_uses_digits[(int)INC] = TRUE; op_uses_digits[(int)DEC] = TRUE;
  op_uses_digits[(int)EQL] = TRUE; op_uses_digits[(int)GRT] = TRUE;
  op_uses_digits[(int)BRB] = TRUE; op_uses_digits[(int)BRF] = TRUE;
  op_uses_digits[(int)PSH] = TRUE; op_uses_digits[(int)JSR] = TRUE;
  op_uses_digits[(int)BLN] = TRUE;

  is_codon[(int)'.'] = TRUE;
  for( i=0; i <= 9; i++ ) {
    op = codon_from_digit(i);
    is_codon[op] = TRUE;
    is_allowed[op] = TRUE;
  }
}

/********************************************************//**
   If pointer is NULL, print error message and exit.
*/
void check_null( void * ptr )
{
  if( ptr == NULL ) {
    fail((char *)"Error: memory allocation failed.");
  }
}

/********************************************************//**
   Print Error message and exit.
*/
void fail( char *msg )
{
  fprintf( stderr,"%s\n",msg );
  exit(1);
}

/********************************************************//**

*/
int mysign( double x )
{
  return( x > 0.0 ? 1 : ( x < 0.0 ? -1 : 0 ));
}

/********************************************************//**
   Ensure identical stream of random variables
*/
void init_mutation_state( long seed )
{
  initstate(seed,(char *)mutation_state,128);
  mutation_bit_int = random_bit_int;
  random_bit_int = 0;
  mutation_gaussian_phase = random_gaussian_phase;
  random_gaussian_phase = 0;
}

/********************************************************//**
*/
void restore_mutation_state()
{
  setstate((char *)mutation_state);
  random_bit_int = mutation_bit_int;
  random_gaussian_phase = mutation_gaussian_phase;
}

/********************************************************//**
*/
void init_evaluation_state( long seed )
{
  initstate(seed,(char *)evaluation_state,128);
  random_bit_int = 0;
  random_gaussian_phase = 0;
}

/********************************************************//**

*/
long my_random()
{
  long n;
  n = random();
  //printf("rand %ld\n",n);
  return( n );
}

/********************************************************//**
   Return 0 or 1, randomly.
*/
int random_bit()
{
  int bit;
  //printf("random_bit()\n");
  if( random_bit_int < 2 ) {
    random_bit_int = random();
  }
  bit = random_bit_int % 2;
  random_bit_int /= 2;
  return( bit );
}

/********************************************************//**
   Return random variable x from a uniform distribution
   with 0.0 <= x < 1.0
*/
double random_uniform()
{
  int n = my_random();
  return( n/MAX_INT_PLUS_ONE );
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
      double U1 = (double) my_random() / MAX_INT_PLUS_ONE;
      double U2 = (double) my_random() / MAX_INT_PLUS_ONE;
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
   Return 1 if the codon is dot or a digit, 0 otherwise.
*/
int is_dot_or_digit( int ch )
{
  return( ch == DOT ||( ch >= '0' && ch <= '9' )
                    ||( ch >= 'A' && ch <= 'Z' ));
}

/********************************************************//**
   Return 1 if the codon is a digit, 0 otherwise.
*/
int codon_is_digit( int ch )
{
  return((ch >= '0' && ch <= '9')||(ch >= 'A' && ch <= 'Z'));
}

/********************************************************//**
   Convert digit to codon by adding '0'.
*/
int codon_from_digit( int digit )
{
  if( digit < 10 ) {
    return( '0' + digit );
  }
  else {
    return( 'A' + digit - 10 );
  }
}

/********************************************************//**
   Convert codon to digit by subtracting '0'.
*/
int digit_from_codon( int codon )
{
  if( codon >= '0' && codon <= '9' ) {
    return( codon - '0' );
  }
  else {
    return( 10 + codon - 'A' );
  }
}

/********************************************************//**
   Assume tp->num_cells has already been set.
   Reset template, ready for copying.
*/
void reset_template( Template *tp )
{
  // tp->value_changed  = FALSE;
  tp->overflow   = FALSE;
  tp->num_focus_bars = 0;
  tp->c=0; tp->b=0; tp->k=0; tp->v=0; tp->p=0;
  tp->cell[0]  =  0;
  tp->bar[0]   =  0;
  tp->codon[0] = '|';
}

/********************************************************//**
   Scan code from fp and return pointer to (new) Code.
*/
Code *code_from_template( Template *tp )
{
  Code *cd;
  unsigned int num,n;
  cd = (Code *)malloc(sizeof(Code));
  check_null( cd );
  cd->num_reg    = tp->num_reg;
  cd->num_cells  = tp->num_cells;
  cd->stack_size = tp->stack_size;
  cd->mem_size   = tp->mem_size;
  cd->num_value  = tp->v;
  cd->num_param  = tp->p;
  cd->level      = tp->level;
  cd->Fc = tp->Fc;
  cd->Nc = tp->Nc;
  cd->num_focus_bars = tp->num_focus_bars;
  for( n=0; n < cd->num_focus_bars; n++ ) {
    cd->Fb[n] = tp->Fb[n];
  }
  num = cd->num_cells;
  cd->cell = (Index *)malloc((num+1)*sizeof(Index));
  check_null( cd->cell );
  memcpy(cd->cell,tp->cell,(num+1)*sizeof(Index));
  num = cd->cell[cd->num_cells]; // num bars
  cd->bar =(Index *)malloc((num+1)*sizeof(Index));
  check_null( cd->bar );
  memcpy(cd->bar,tp->bar,(num+1)*sizeof(Index));
  cd->last_codon = cd->bar[cd->cell[cd->num_cells]];
  cd->codon=(char *)malloc((cd->last_codon+2)*sizeof(char));
  check_null( cd->codon );
  memcpy(cd->codon,tp->codon,(cd->last_codon+1)*sizeof(char));
  cd->codon[cd->last_codon+1] = '\0';
  cd->ival=(Index *)malloc((cd->last_codon+1)*sizeof(Index));
  check_null( cd->ival );
  memcpy(cd->ival,tp->ival,(cd->last_codon+1)*sizeof(Index));
  cd->ival[cd->last_codon] = 0;

  if( cd->num_value > 0 ) {
    cd->fval=(Floating *)malloc(cd->num_value*sizeof(Floating));
    check_null( cd->fval );
    memcpy(cd->fval,tp->fval,cd->num_value*sizeof(Floating));
  }
  else {
    cd->fval = NULL;
  }
  if( cd->num_param > 0 ) {
    cd->param=(Index *)malloc(cd->num_param*sizeof(Index));
    check_null( cd->param );
    memcpy(cd->param,tp->param,(cd->num_param)*sizeof(Index));
  }
  else {
    cd->param = NULL;
  }
#ifdef DEBUG
  //print_code( cd,stdout );
#endif
  return( cd );
}

/********************************************************//**
   Return a new code, with all cells empty.
*/
Code * empty_code(
                  int num_reg,
                  int num_cells,
                  int stack_size,
                  int mem_size
                 )
{
  Template temp;
  temp.num_reg    = num_reg;
  temp.num_cells  = num_cells;
  temp.stack_size = stack_size;
  temp.mem_size   = mem_size;
  temp.overflow   = FALSE;
  temp.level.m    = BLOCK;
  //while( power2[temp.level.m-BLOCK] < num_cells ) {
  while( power2[temp.level.m-BLOCK] <= num_cells ) {
    temp.level.m++;
  }
  temp.level.g = ALIGNED;
  reset_template( &temp );
  temp.Fc = 0;
  temp.Nc = temp.num_cells;
  temp.num_focus_bars = 0;
  while( temp.c < num_cells ) {
    temp.codon[++temp.k] = '|';
    temp.bar[++temp.b]   = temp.k;
    temp.cell[++temp.c]  = temp.b;
  }
  // Darwin: temp.ival isn't initialised anywhere else and it was causing weird bugs so I added this line.
  memset(temp.ival, 0, (temp.bar[temp.cell[temp.num_cells]] + 1)*sizeof(Index));
  return( code_from_template( &temp ));
}

/********************************************************//**

*/
void free_code( Code *cd )
{
  if( cd != NULL ) {
    if( cd->codon != NULL ) {
      free( cd->codon ); cd->codon = NULL;
    }
    if( cd->bar != NULL ) {
      free( cd->bar ); cd->bar = NULL;
    }
    if( cd->cell  != NULL ) {
      free( cd->cell ); cd->cell  = NULL;
    }
    if( cd->ival != NULL ) {
      free( cd->ival ); cd->ival = NULL;
    }
    if( cd->fval != NULL ) {
      free( cd->fval ); cd->fval = NULL;
    }
    if( cd->param != NULL ) {
      free( cd->param ); cd->param = NULL;
    }
    free( cd );
  }
}

/********************************************************//**
   Return pointer to a new, empty channel
*/
Channel * new_channel()
{
  Channel *chan;
  chan = (Channel *)malloc(sizeof(Channel));
  check_null( chan );
  chan->max_message = 64; //  initial values
  chan->max_item  = 1024; // (will be expanded when necessary)
  chan->val=(Floating *)malloc(chan->max_item*sizeof(Floating));
  check_null( chan->val );
  chan->index  = (int *)malloc(chan->max_message*sizeof(int));
  check_null( chan->index );
  chan->length = (int *)malloc(chan->max_message*sizeof(int));
  check_null( chan->length );
  chan->domain = UNRESTRICTED;
  chan->min_length = -1;
  chan->max_length = 65536;
  chan->max_im = -1;   //     no limit
  chan->max_om = -1;   //     no limit
  chan->im = 0;        //  input message
  chan->om = 0;        // output message
  chan->si = 0;        //   scan-item
  chan->wi = 0;        //  write-item
  chan->length[0]= -1; // (initially) non-message
  chan->index[0] =  0;
  chan->index[1] =  0;
  clear_message(chan,1);

  return( chan );
}

/********************************************************//**
   Return TRUE,  if the messages of the shorter channel are identical
                 to the leading messages of the longer channel;
          FALSE, otherwise.
*/
Boolean matching_channel( Channel *chan0, Channel *chan1 )
{
  int om;
  if( chan0 == NULL || chan1 == NULL ) {
    return( chan0 == NULL && chan1 == NULL );
  }
  om = chan0->om;
  if( chan1->om < om ) {
    om = chan1->om;
  }
  return(   !memcmp(chan0->index, chan1->index, (om+2)*sizeof(int)) 
         && !memcmp(chan0->length,chan1->length,(om+1)*sizeof(int))
         && !memcmp(chan0->val,chan1->val,chan0->index[om+1]*sizeof(Floating)));
}

/********************************************************//**
   Free channel and all the space it occupies.
*/
void free_channel( Channel *chan )
{
  if( chan != NULL ) {
    if(chan->val    != NULL) { free(chan->val); }
    if(chan->index  != NULL) { free(chan->index); }
    if(chan->length != NULL) { free(chan->length); }
    free( chan );
  }
}

/********************************************************//**
   Fetch input and return TRUE or FALSE
*/
Boolean fetch_input( Channel *in )
{
  if( in == NULL ) {
    return( FALSE );
  }
  else {
    int max_im = in->om;
    if( in->max_im >= 0 && in->max_im < max_im ) {
      max_im = in->max_im;
    }
    if( in->im >= max_im ) {
      in->im = max_im;
      in->si = in->index[in->im+1]; // end of message
      return( FALSE );
    }
    else {
      in->im++;
      in->si = in->index[in->im]; // start of message
      return( in->length[in->im] >= 0 );
    }
  }
}

/********************************************************//**
   Return TRUE  if we have tried to input beyond the
                last available message,
          FALSE otherwise.
*/
Boolean input_exhausted( Channel *in )
{
  if( in == NULL ) {
    return TRUE;
  }
  else {
    int max_im = in->om;
    if( in->max_im >= 0 && in->max_im < max_im ) {
      max_im = in->max_im;
    }
    return(    in->im > max_im
           ||( in->im == max_im && in->si == in->index[in->im+1]));
  }
}

/********************************************************//**
   Set chan->om to be one less than the specified om_1,
   so the next message written will be om_1.
<pre>
              wi
               |
     +----+---+----+
 val |    |   |    |
     +----+---+----+
      |    |  /
     +-+--+-+-+-------+
index|0|..| | |       |
     +-+--+-+-+-------+
           | |
          om(om+1)
</pre>
*/
void clear_message( Channel *chan, int om_1 )
{
  if( chan != NULL && om_1 > 0 &&(   chan->max_om < 0
                                  || chan->max_om >= om_1 )) {
#ifdef DEBUG
    //printf("Clearing message %d for output.\n",om_1);
#endif
    // chan->max_message must be at least om_1+2
    accommodate_message( chan,om_1 );

    // if chan->om smaller than specified (om_1)-1,
    // mark the intervening messages as non-existent.
    while( chan->om < om_1 - 1 ) {
      chan->om++;
      chan->index[chan->om+1] = chan->index[chan->om];
      chan->length[chan->om+1]= -1;
    }
    chan->om = om_1 - 1; // next message written will be om_1
    chan->wi = chan->index[om_1];
    chan->length[om_1] = -1; // message initially empty
    // chan->index[om_1+1] = chan->wi; // IS THIS NEEDED ??????
    // cannot read beyond what has been written
    if( chan->im > chan->om ) {
      chan->im = chan->om;
      chan->si = chan->index[chan->im];
    }
  }
}

/********************************************************//**
   Set channel so that the next message input will be im_1,
   and any scan in the meantime will fail.
*/
void cue_message( Channel *chan, int im_1 )
{
  if( chan != NULL && im_1 <= chan->om+1 ) { // chan->max_om ?????
#ifdef DEBUG
    //printf("  Cueing message %d for  input.\n",im);
#endif
    chan->im = im_1 - 1;          // next message input will be im_1
    chan->si = chan->index[im_1]; // prior to input, scan will fail
  }
}

/********************************************************//**
   Flush output buffer and send its contents as a new message.
   Increment out->om, ready to start writing next message.
*/
void output_message( Channel *out )
{
  if( out != NULL ) {
    if( out->max_om < 0 || out->max_om > out->om ) {
      out->om++;
      accommodate_message( out,out->om );
      out->length[out->om]  = out->wi - out->index[out->om];
      out->index[out->om+1] = out->wi;
      out->length[out->om+1]= -1;
    }
    else { // output an empty message and clear output buffer
      /*
      out->om = out->max_om+1;
      out->length[out->om] = -1;
      out->wi = out->index[out->om];
      out->index[out->om+1] = out->wi;
      */
    }
  }
}

/********************************************************//**
   output a non-existant message.
   Increment out->om, ready to start writing next message.
*/
void output_non_message( Channel *out )
{
  if( out != NULL ) {
    if( out->max_om < 0 || out->max_om > out->om ) {
      out->om++;
      accommodate_message( out,out->om );
      out->length[out->om]  = -1;
      out->index[out->om+1] = out->index[out->om];
      out->length[out->om+1]= -1;
    }
  }
}

/********************************************************//**
   Delete messages from new_m up to old_m-1, shift values
   downwards and adjust im,om,is and iw accordingly.
*/
void shift_message( Channel *chan, int new_m, int old_m )
{
  if( chan != NULL && new_m < old_m && new_m <= chan->om ) {

    if( chan->om >= old_m ) { // there is content to be shifted
      int length = chan->wi - chan->index[old_m];
      int shift  = chan->index[new_m] - chan->index[old_m];
      int n;
      if( length > 0 && shift < 0 ) {
        memmove( &chan->val[chan->index[new_m]],
                 &chan->val[chan->index[old_m]],
                     length*sizeof(Floating));
      }
      memmove( &chan->length[new_m],&chan->length[old_m],
                 (chan->om+1 - old_m)*sizeof(int));

      for( n=0; n <= chan->om+1 - old_m; n++ ) {
        chan->index[new_m+n] = chan->index[old_m+n]+shift;
      }
      //printf("om %d im %d si %d wi %d |",chan->om,chan->im,
      //                                   chan->si,chan->wi);
      chan->om += ( new_m - old_m );
      chan->wi += shift;
      if( chan->im >= old_m ) { // already reached message old_m
        chan->im += ( new_m - old_m );
        chan->si += shift;
      }
      else if( new_m > 0 ) {    // cue message new_m for input
        chan->im = new_m-1;
        chan->si = chan->index[new_m];
      }
      else {                    // cue message 1 for input
        chan->im = 0;
        chan->si = 0;
      }
    }
    else { // no content to shift, clear message new_m
      chan->om = ( new_m > 0 ) ? new_m-1 : 0 ;
      chan->wi = chan->index[chan->om+1];
      chan->length[chan->om+1] = -1;
      if( chan->im >= new_m ) {
        chan->im = ( new_m > 0 ) ? new_m-1 : 0 ;
        chan->si = chan->index[chan->im+1];
      }
    }
  }
}

/********************************************************//**
   Return TRUE  if messages m0 in chan0 and m1 in chan1
   are identical, FALSE otherwise.
*/
int same_message( Channel *chan0, int m0,
                  Channel *chan1, int m1 )
{
  int i;
  if(   chan0 == NULL   ||   chan0->om < m0
     || chan1 == NULL   ||   chan1->om < m1
     || chan0->length[m0] != chan1->length[m1] ) {
    //printf("messages not the same length\n");
    return FALSE;
  }
  for( i=0; i < chan0->length[m0]; i++ ) {
    if(   chan0->val[chan0->index[m0]+i]
       != chan1->val[chan1->index[m1]+i] ) {
      //printf("messages not the same contents\n");
      return FALSE;
    }
  }
  //printf("messages the same\n");
  return TRUE;
}

/********************************************************//**
   increment copy->om and copy message m from chan
   to message copy->om of copy
*/
void copy_message( Channel *copy, Channel *chan, int m )
{
  if(   chan != NULL && m <= chan->om
     && copy != NULL &&(   copy->max_om < 0
                        || copy->max_om > copy->om )) {
    copy->om++;
    accommodate_message( copy,copy->om );
    copy->length[copy->om] = chan->length[m];
    copy->index[copy->om]  = chan->index[m];
    copy->index[copy->om+1]= chan->index[m+1];
    copy->wi = copy->index[copy->om+1];
    if( chan->length[m] > 0 ) {
      memcpy(&copy->val[copy->index[copy->om]],
             &chan->val[chan->index[m]],
              chan->length[m]*sizeof(Floating));
    }
  }
}

/********************************************************//**
   Scan item from input
*/
Boolean scan_value( Channel *chan, double *p_value )
{
  if( chan != NULL && chan->si < chan->index[chan->im+1] ) {
    *p_value = chan->val[chan->si++];
    return( TRUE );
  }
  else {
    return( FALSE );
  }
}

/********************************************************//**
   Write item to output channel
  (expanding the size of value array, if necessary).
*/
void write_value( Channel *chan, double value )
{
  if( chan != NULL &&( chan->max_om < 0 || chan->max_om > chan->om )) {
    if( chan->wi >= chan->max_item ) {
      chan->max_item *= 2;
      chan->val=(Floating *)realloc(chan->val,chan->max_item*sizeof(Floating));
      check_null( chan->val );
    }
    chan->val[chan->wi++] = value;
  }
}

/********************************************************//**
   If om+2 exceeds chan->max_message,
   increase chan->max_message to the next higher power of 2
   and reallocate memory for chan->index and chan->length
*/
void accommodate_message( Channel *chan, int om )
{
  if( chan->max_message < om+2 ) {
    while( chan->max_message < om+2 ) {
      chan->max_message *= 2;
    }
    chan->index =(int *)realloc(chan->index, chan->max_message*sizeof(int));
    check_null( chan->index );
    chan->length=(int *)realloc(chan->length,chan->max_message*sizeof(int));
    check_null( chan->length );
  }
}

/********************************************************//**
   Create a new, empty agent.
*/
Agent * new_agent()
{
  Agent *agt;
  agt = (Agent *)malloc(sizeof(Agent));
  check_null( agt );
  agt->cd  = NULL;

  agt->reg   = NULL;
  agt->stack = NULL;
  agt->mem   = NULL;

  agt->call_stack = NULL;
  agt->bstack = NULL;
  agt->bstack_size = 0;
  agt->call_stack_size = 0;
  agt->stack_items = 0;

  agt->running = FALSE;

  return( agt );
}

/********************************************************//**
   Reset the state of the agent.
*/
void reset_agent( Agent *agt )
{
  Code  *cd = agt->cd;
  if( agt->stack == NULL ) {
    agt->stack=(Floating *)malloc(cd->stack_size*sizeof(Floating));
    check_null( agt->stack );
  }
  memset(agt->stack,0,cd->stack_size*sizeof(Floating));
  if( agt->cd->num_reg > 0 ) {
    if( agt->reg   == NULL ) {
      agt->reg = (Floating *)malloc(cd->num_reg*sizeof(Floating));
      check_null( agt->reg );
    }
    memset(agt->reg,0,cd->num_reg*sizeof(Floating));
  }
  if( cd->mem_size > 0 ) {
    if ( agt->mem  == NULL ) {
      agt->mem = (Floating *)malloc(cd->mem_size*sizeof(Floating));
      check_null( agt->mem );
    }
    memset(agt->mem,0,cd->mem_size*sizeof(Floating));
  }
  if( agt->bstack == NULL ) {
    agt->bstack_size = 16;
    agt->bstack=(Boolean *)malloc(agt->bstack_size*sizeof(Boolean));
    check_null( agt->bstack );
  }
  if( agt->call_stack == NULL ) {
    agt->call_stack_size = 16;
    agt->call_stack =(Snode *)malloc(agt->call_stack_size*sizeof(Snode));
    check_null( agt->call_stack );
  }
  agt->fp = -1;

  // keep code, num_codons, *_size, *_power

  agt->stack_items = 0;
  agt->sp = cd->stack_size-1;
  agt->bp = -1;

  agt->running = TRUE;
  agt->step    = 0;

  push_stack_frame( agt, cd->num_cells-1 );
}

/********************************************************//**
   Free memory used by reg, stack, mem and call_stack.
*/
void compress_agent( Agent *agt )
{
  if( agt != NULL ) {
    if( agt->reg   != NULL ) {
      free( agt->reg );     agt->reg   = NULL;
    }
    if( agt->stack != NULL ) {
      free( agt->stack );   agt->stack = NULL;
    }
    if( agt->mem   != NULL ) {
      free( agt->mem );     agt->mem   = NULL;
    }
    if( agt->bstack != NULL ) {
      free( agt->bstack );
      agt->bstack_size = 0;
    }
    if( agt->call_stack != NULL ) {
      free( agt->call_stack );
      agt->call_stack_size = 0;
    }
    agt->call_stack = NULL;
  }
}

/********************************************************//**
*/
void next_boolean( Agent *agt )
{
  agt->bp++;
  if( agt->bp >= agt->bstack_size ) {
    agt->bstack_size *= 2;
    agt->bstack=(Boolean *)realloc(agt->bstack,agt->bstack_size*sizeof(Snode));
    check_null( agt->bstack );
  }
  agt->bstack[agt->bp] = TRUE; // default value
}

/********************************************************//**
   Adjust stack pointer up or down, wrapping if necessary.
*/
int stack_adjust( Agent *agt, int epsilon )
{
  int stack_size = agt->cd->stack_size;
  if( epsilon < 0 ) {
    return((agt->sp+stack_size+epsilon)%stack_size);
  }
  else {
    return((agt->sp + epsilon)% stack_size);
  }
}

/********************************************************//**
   Return value on top of stack.
*/
double top_of_stack( Agent *agt )
{
  return( agt->stack[agt->sp] );
}

/********************************************************//**
   If an INP has just failed, undo its effect
*/
void undo_input( Agent *agt, Channel *in )
{
  Snode *sf;                        // stack frame
  //printf("UNDO INPUT\n");
  if( agt->bstack[agt->bp] == TRUE ) {
    if( in != NULL && in->im > 0 ) {
      in->im--;
      in->si = in->index[in->im+1]; // end of message
    }
  }
  agt->bp--;
  sf = agt->call_stack + agt->fp;
  sf->k--;
  while( is_dot_or_digit( agt->cd->codon[sf->k] )) {
    sf->k--;
  }
}

/********************************************************//**
   undo_output
*/
void undo_output( Agent *agt, Channel *out )
{
  Snode *sf;                        // stack frame
  if( out != NULL ) {
    if( out->om > 0 ) {
      out->om--;
      out->length[out->om+1]= -1;
    }
  }
  sf = agt->call_stack + agt->fp;
  sf->k--;
  while( is_dot_or_digit( agt->cd->codon[sf->k] )) {
    sf->k--;
  }
}

/********************************************************//**
   Execute next instruction.
*/
int step( Agent *agt, Channel *in, Channel *out )
{
  Floating *p_reg;
  Code  *cd = agt->cd;   // code
  Snode *sf;             // stack frame
  int   mem_index=-1;    // sentinel
  int   sp1;
  Index op;              // operator
  Index d;               // cell, digit pointer

  sf = agt->call_stack + agt->fp;
  agt->step++;

  d = sf->k;            // start of instruction
  sf->k += cd->ival[d]; // end of instruction
#ifdef ASSERT
  assert( sf->k <= cd->last_codon+1 );
#endif
  op = cd->codon[sf->k++];

  switch( op ) {

  case ADD:
    if( agt->stack_items > 1 ) {
      sp1 = stack_adjust( agt, -1 );
      if( isfinite(agt->stack[sp1])) {
        agt->stack[sp1] += agt->stack[agt->sp];
      }
      agt->stack[agt->sp] = 0.0;
      agt->stack_items--;
      agt->sp = sp1;
    }
    break;

  case AND:
    if( agt->bp > sf->bp_reset+1 ) {
      agt->bp--;
      agt->bstack[agt->bp] &= agt->bstack[agt->bp+1];
    }
    break;

  case ASN:
    if( agt->stack_items > 0 ) {
      if( agt->stack[agt->sp] < -1.0 ) {
        agt->stack[agt->sp] = -1.0;
      }
      else if( agt->stack[agt->sp] > 1.0 ) {
        agt->stack[agt->sp] = 1.0;
      }
      agt->stack[agt->sp] = asin(agt->stack[agt->sp]);
    }
    break;

  case BLN:
    if( cd->codon[d] == '.' ) { // return
      return_from_jump( agt );
    }
    else if( sf->k > d+1 ) {    // halt
      while( agt->running ) {
        return_from_jump( agt );
      }
    }
    else {
      sf->b++;
      if( sf->b == cd->cell[sf->c+1] ) {
        return_from_jump( agt );
      }
    }
    break;

  case BRB:
    {
      int offset=0;
      if( sf->k > d+1 ) {
        offset = cd->ival[d+1];
      }
      if( agt->bstack[agt->bp] ) {
        sf->b -= offset;
        if( sf->b < cd->cell[sf->c] ) {
          sf->b = cd->cell[sf->c];
        }
        sf->k = cd->bar[sf->b]+1;
      }
      agt->bp = sf->bp_reset;
      agt->bstack[agt->bp] = TRUE;
    }
    break;

  case BRF:
    {
      Boolean branch_forward = agt->bstack[agt->bp];
      int offset=0;
      agt->bp = sf->bp_reset;
      agt->bstack[agt->bp] = TRUE;
      if( sf->k > d+1 ) {
        offset = cd->ival[d+1];
      }
      if( branch_forward ) {
        sf->b += offset+1;
        if( sf->b >= cd->cell[sf->c+1] ) {
          return_from_jump( agt );
        }
        else {
          sf->k = cd->bar[sf->b]+1;
        }
      }
    }
    break;

  case CPY:
    if( agt->stack_items > 0 ) {
      push_to_stack(agt,agt->stack[agt->sp]);
    }
    break;

  case DEC:
    p_reg = set_reg( agt,sf,d );
    *p_reg -= 1.0;
    break;

  case EQL:
    p_reg = set_reg( agt,sf,d );
    next_boolean( agt );
    agt->bstack[agt->bp] =(   agt->stack_items > 0
                           && agt->stack[agt->sp] == *p_reg );
    break;

  case GRT:
    p_reg = set_reg( agt,sf,d );
    next_boolean( agt );
    if( agt->stack_items > 0 ) {
      agt->bstack[agt->bp] = ( *p_reg > agt->stack[agt->sp]);
    }
    else {
      agt->bstack[agt->bp] = ( *p_reg > 0.0 );
    }
    break;

  case INC:
     p_reg = set_reg( agt,sf,d );
    *p_reg += 1.0;
    break;

  case EXP:
    if( agt->stack_items > 0 ) {
      // we rely on IEEE standard regarding infinity
      agt->stack[agt->sp] = exp(agt->stack[agt->sp]);
    }
    break;

  case GET:
    p_reg = set_reg( agt,sf,d );
    push_to_stack( agt,*p_reg );
    break;

  case INP:
    next_boolean( agt );
    agt->bstack[agt->bp] = fetch_input( in );
    break;

  case JSR:
    {
      Index c=0;
      if( sf->k > d+1 ) {
        c = cd->ival[d+1];
      }
      if( c < cd->num_cells ) {
        push_stack_frame( agt, c );
      }
    }
    break;

  case LOD:
    p_reg = set_reg( agt,sf,d );
    if( cd->mem_size > 0 && isfinite(*p_reg)) {
      mem_index = get_mem_index(*p_reg,cd->mem_size);
      push_to_stack(agt,agt->mem[mem_index]);
    }
    else { // zero memory, or infinite memory location
      push_to_stack( agt,0.0 ); // dummy value
    }
    break;

  case STO:
    next_boolean( agt );
    if( agt->stack_items == 0 ) {
      agt->bstack[agt->bp] = FALSE;
    }
    else {
      p_reg = set_reg( agt,sf,d );
      if( cd->mem_size > 0 && isfinite(*p_reg)) {
        mem_index = get_mem_index(*p_reg,cd->mem_size);
        agt->mem[mem_index] = agt->stack[agt->sp];
      }
      // if zero memory, or infinite memory location
      // just pop top item from stack and discard
      agt->stack[agt->sp] = 0.0;
      agt->sp = stack_adjust( agt, -1 );
      agt->stack_items--;
    }
    break;

  case LOG:
    if( agt->stack_items > 0 ) {
      if( agt->stack[agt->sp] <= 0.0 ) {
        agt->stack[agt->sp] = -INFINITY;
      }
      else {
        agt->stack[agt->sp] = log(agt->stack[agt->sp]);
      }
    }
    break;

  case MLT:
    if( agt->stack_items > 1 ) {
      double result;
      sp1 = stack_adjust( agt, -1 );
      result = agt->stack[sp1]*agt->stack[agt->sp];
      if( isnan(result)) {// infinity times zero is zero
        result = 0.0;
      }
      agt->stack[sp1] = result;
      agt->stack[agt->sp] = 0.0;
      agt->stack_items--;
      agt->sp = sp1;
    }
    break;

  case MOD: // ..(x)(m) -> ..(x/m)(x % m)
    if( agt->stack_items > 1 ) {
      sp1 = stack_adjust( agt, -1 );
      if( agt->stack[agt->sp] == 0.0 ) { 
        if( agt->stack[sp1] > 0.0 ) {
          agt->stack[sp1] =  INFINITY;
        }
        else if( agt->stack[sp1] < 0.0 ) {
          agt->stack[sp1] = -INFINITY;
        }
      }
      else {
        double quotient,remainder,dividend;
        double divisor = agt->stack[agt->sp];
        dividend = agt->stack[sp1];
        if( isfinite( divisor )) {
          if( isfinite( dividend )) {
            quotient  = floor(dividend/divisor);
            remainder = dividend - quotient*divisor;
          }
          else {
            quotient  =(divisor > 0.0) ? dividend : -dividend;
            remainder = 0.0;
          }
        }
        else { // divisor is +infty or -infty
          remainder = dividend;
          quotient = isfinite(dividend) ? 0.0 : dividend*divisor;
        }
        agt->stack[sp1] = quotient;
        agt->stack[agt->sp] = remainder;
      }
    }
    break;

  case NEG:
    if( sf->k > d+1 ) { // #- combination
      //if( cd->codon[sf->k-2] == PSH ) {      // 
      if( is_dot_or_digit(cd->codon[d])) { // value specified
        push_to_stack(agt,cd->fval[cd->ival[d+1]]);
      }
      else {
        push_to_stack(agt,-0.0);           // negative zero
      }
    }
    else if( agt->stack_items > 0 ) { // stand-alone negation
      agt->stack[agt->sp] = -agt->stack[agt->sp];
    }
    break;

  case NOT:
    agt->bstack[agt->bp] = !agt->bstack[agt->bp];
    break;

  case OR:
    if( agt->bp > sf->bp_reset+1 ) {
      agt->bp--;
      agt->bstack[agt->bp] |= agt->bstack[agt->bp+1];
    }
    break;

  case OUT:
    output_message( out );
    break;

  case PAY:
    // NOT IMPLEMENTED
    break;

  case PLR: // ..(y)(x) -> ..(theta)(r)
    if( agt->stack_items > 1 ) {
      double x,y;
      sp1 = stack_adjust( agt, -1 );
      x = agt->stack[agt->sp];
      y = agt->stack[sp1];
      agt->stack[sp1]      = atan2(y,x);
      agt->stack[agt->sp] = sqrt(x*x + y*y);
    }
    break;

  case POP:
    if( agt->stack_items > 0 ) {
      agt->stack[agt->sp] = 0.0;
      agt->sp = stack_adjust( agt, -1 );
      agt->stack_items--;
    }
    break;

  case PSH:
    if( sf->k > d+1 ) {
      push_to_stack(agt,cd->fval[cd->ival[d+1]]);
    }
    else {
      push_to_stack(agt,0.0);
    }
    break;

  case PUT:
    next_boolean( agt );
    if( agt->stack_items == 0 ) {
      agt->bstack[agt->bp] = FALSE;
    }
    else {
      p_reg = set_reg( agt,sf,d );
     *p_reg = agt->stack[agt->sp];
      agt->stack[agt->sp] = 0.0;
      agt->sp = stack_adjust( agt, -1 );
      agt->stack_items--;
    }
    break;

  case RAN:
    push_to_stack( agt, random_uniform());
    //printf("? %1.4f\n",agt->stack[agt->sp]);
    break;

  case RCP:
    if( agt->stack_items > 0 ) {
      agt->stack[agt->sp] = 1.0 / agt->stack[agt->sp];
    }
    break;

  case RND:
    if( agt->stack_items > 0 ) {
      agt->stack[agt->sp] = floor(0.5+agt->stack[agt->sp]);
    }
    break;

  case ROT:
    if( agt->stack_items > 2 ) {
      int sp2;
      double tmp;
      sp1 = stack_adjust( agt, -1 );
      sp2 = stack_adjust( agt, -2 );
      tmp = agt->stack[agt->sp];
      agt->stack[agt->sp] = agt->stack[sp1];
      agt->stack[sp1]     = agt->stack[sp2];
      agt->stack[sp2]     = tmp;
    }
    break;

  case SCN: // move one value from input channel to stack
    next_boolean( agt );
    /*
    if( in->si == 5 ) {
      printf("SI\n");
    }
    if( in->im > 1 ) {
      printf("si=%d,im=%d,index=%d\n",in->si,in->im,in->index[in->im+1]);
    }
    if( in->index[in->im+1] == 17 ) {
      printf("INDEX\n");
    }
    */
    if( in != NULL && in->si < in->index[in->im+1] ) {
      push_to_stack(agt,in->val[in->si++]);
      //agt->sp = stack_adjust( agt, +1 );
      //agt->stack[agt->sp] = in->val[in->si++];
    }
    else {
      agt->bstack[agt->bp] = FALSE;
    }
    break;

  case SQT:
    if( agt->stack_items > 0 ) {
      if( agt->stack[agt->sp] < 0.0 ) {
        agt->stack[agt->sp] = 0.0;
      }
      else {
        agt->stack[agt->sp] = sqrt(agt->stack[agt->sp]);
      }
    }
    break;

  case SWP:
    if( agt->stack_items > 1 ) {
      double tmp = agt->stack[agt->sp];
      sp1 = stack_adjust( agt, -1 );
      agt->stack[agt->sp] = agt->stack[sp1];
      agt->stack[sp1]     = tmp;
    }
    break;

  case TNH:
    if( agt->stack_items > 0 ) {
      agt->stack[agt->sp] = tanh(agt->stack[agt->sp]);
    //agt->stack[agt->sp] = 1.0/(1.0 + exp(-agt->stack[agt->sp]));
    }
    break;

  case TRG: // ..(theta)(r) -> (y)(x)
    if( agt->stack_items > 1 ) {
      double radius,theta,x=0.0,y=0.0;
      sp1 = stack_adjust( agt, -1 );

      theta  = agt->stack[sp1];
      radius = agt->stack[agt->sp];

      if( isfinite( radius )) {
        if( isfinite( theta )) {
          x = radius * cos(theta);
          y = radius * sin(theta);
        }
        else {
          x = 0.0;
          y = 0.0;
        }
      }
      else { // infinite radius
        if( isfinite(theta)) {
          switch( mysign(cos(theta))) {
           case  1: x =  radius; break;
           case  0: x =       0; break;
           case -1: x = -radius; break;
          }
          switch( mysign(sin(theta))) {
           case  1: y =  radius; break;
           case  0: y =       0; break;
           case -1: y = -radius; break;
          }
        }
        else {
          x = radius;
          y = radius;
        }
      }
      agt->stack[sp1]     = y;
      agt->stack[agt->sp] = x;
    }
    break;

  case WRT:
    next_boolean( agt );
    if( agt->stack_items == 0 ) {
      agt->bstack[agt->bp] = FALSE;
    }
    else {
      write_value( out,agt->stack[agt->sp] );
      agt->stack[agt->sp] = 0.0;
      agt->sp = stack_adjust( agt, -1 );
      agt->stack_items--;
    }
    break;
  }

  if( globals.learning_rate > 0.0 ) {
    grad_step( op,agt,sf,out,d,mem_index );
  }

  return( op );
}

/********************************************************//**
   Execute next instruction.
*/
int get_mem_index(
                  Floating raw_index,
                  int mem_size
                 )
{
  int mem_index =(int)floor(raw_index + 0.5);
  mem_index = mem_index % mem_size;
  while( mem_index < 0 ) {
    mem_index += mem_size;
  }
  return( mem_index );
}

/********************************************************//**
   Push a new stack frame onto call stack.
*/
void push_stack_frame( Agent *agt, int c )
{
  Snode *sf;

  agt->fp++;
  if( agt->fp >= agt->call_stack_size ) {
    if( agt->fp >= MAX_CALL_STACK ) {
      agt->fp--;
      while( agt->running ) {
        return_from_jump( agt );
      }
      return;
    }
    agt->call_stack_size *= 2;
    agt->call_stack =(Snode *)realloc(agt->call_stack,
                      agt->call_stack_size*sizeof(Snode));
    check_null( agt->call_stack );
  }
  sf = agt->call_stack + agt->fp;

  sf->c = c;
  sf->b = agt->cd->cell[c];
  sf->k = agt->cd->bar[sf->b]+1;
  sf->r = -1;                    // local register
  sf->local_reg = 0.0;
  next_boolean( agt );
  //agt->bstack[agt->bp] = TRUE;
  sf->bp_reset = agt->bp;
}

/********************************************************//**
   Push value to stack.
*/
void push_to_stack( Agent *agt, double value )
{
  agt->sp = stack_adjust( agt, +1 );
  agt->stack[agt->sp] = value;
  if( agt->stack_items < agt->cd->stack_size ) {
    agt->stack_items++;
  }
}

/********************************************************//**
   Return from jump by popping stack frame.
*/
void return_from_jump( Agent *agt )
{
  int return_boolean = agt->bstack[agt->bp];
  agt->bp = agt->call_stack[agt->fp].bp_reset;
  agt->fp--;
  if( agt->fp < 0 ) {
    push_stack_frame( agt, agt->cd->num_cells-1 );
    agt->running = FALSE;
  }
  agt->bstack[agt->bp] = return_boolean;
}

/********************************************************//**
   Update the current register, if appropriate,
   and return pointer to current register
*/
Floating *set_reg( Agent *agt, Snode *sf, int d )
{
  Floating *p_reg;
  if( agt->cd->num_reg == 0 || agt->cd->codon[d] == '.' ) {
    sf->r = -1;                    // local register
  }
  else if( sf->k > d+1 ) {
    sf->r = agt->cd->ival[d+1];    // global register
  }
  if( sf->r < 0 ) {
    p_reg = &sf->local_reg;        // local
  }
  else {
    p_reg = &agt->reg[sf->r];      // global
  }
  return( p_reg );
}
