/** \file point.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "step.h"
#include "scan_print.h"
#include "point.h"

#define  MAX_INT   2147483647

#define  STRETCH_OFFSETS    1

#define  CAN_INSERT_BRANCH  1

#define STACK_REG    0   /*     [x](r)[g=][x]  */
#define REG_VAL      1   /* (val).#(r)[g=][!]  */
#define STACK_VAL    2   /* (val).# .>[g=]     */
#define STACK_STACK  3   /*         .>[g=][<]  */        
#define REG_REG      4   /*    (r)<(r)[g=][!]  */

double pseudo_cauchy(double median);
int    random_offset(int offset,int max_offset);
void   focus_weights(int w[],int c0,int Nc,int num_cells);
void adjust_fraction(Fraction *fr,double delta_num,int force_change);

void   negate_push_value(Template *tp);
Fraction random_fraction(Template *tp,int force_integer,int force_dot);
int      random_register(Template *tp);
void     insert_register(Template *tp,int r);
void        insert_query(Template *tp,int query);

void        fix_register(Template *tp,Code *cd,Clip *cp,int must_change);
void          fix_offset(Template *tp,Code *cd,Clip *cp,int must_change);
void      fix_jump_index(Template *tp,Code *cd,Clip *cp,int must_change);
int        align_bar(int b1_min,int b1_max,int b0_min,int b0_max,int b0);

#ifdef DEBUG
char ds[64*MAX_CODE];
int  di = 0;
int  print_debug = FALSE;
#endif

/********************************************************//**
   Return random variable between 1 and max (inclusive)
   from a geometric distribution with the specified mean
*/
int random_geometric( double mean, int max )
{
  double x;
  int n=0;
  if( mean > 0.0 ) {
    do {
      x = 1.0 - random_uniform();
      n = 1+(int)floor(log(x)/log(1.0-1.0/mean));
    } while( n > max );
  }
  return( n );
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
   Return a random number between 1 and max_offset, with
   a preference for values close (but not equal) to offset.
*/
int random_offset(
                  int offset,
                  int max_offset
                 )
{
  int longer_side;
  int j,k;
  longer_side = offset;
  if( max_offset - offset > offset ) {
    longer_side = max_offset - offset;
  }
  if( offset > 0 &&(random()% max_offset < offset)) {
    // choose a lower offset
    j = my_random()% offset;
    k = offset-1 - my_random()%(longer_side+1);
    return(( j > k ) ? j : k );
  }
  else if( max_offset > offset ) {
    // choose a higher offset
    j = offset+1 + my_random()%(max_offset-offset);
    k = offset+1 + my_random()%(longer_side+1);
    return(( j < k ) ? j : k );
  }
  else {
    return( offset );
  }
}

/********************************************************//**
   Generate a random jump index appropriate to cell b
*/
int choose_jump_index( int c, int num_cells )
{
  int w[MAX_CELL];
  int sum_weight;
  int c0,m;
  // more likely to jump to cells within same block
  focus_weights( w,c,1,num_cells );
  // increased weight to lower cells
  for( c0=0; c0 < c; c0++ ) {
    w[c0] *= 2;
  }
  // reduced weight for the cell to call itself
  if( w[c] > 1 ) {
    w[c] /= 2;
  }
  sum_weight = 0;
  for( c0=0; c0 < num_cells; c0++ ) {
    sum_weight += w[c0];
  }
  m = my_random()% sum_weight;
  c0 = -1;
  while( m >= 0 ) {
    m -= w[++c0];
  }
  return( c0 );
}

/********************************************************//**
   Mutate fraction fr.
   If POINT or TUNE, with probability 1/2 generate
      a new value (similar size and precision to old value)
   Otherwise, mutate the value (possibly change precision)
*/
void mutate_fraction(
                     Fraction *fr,
                     int    level_m,
                     double scale,
                     int    singular
                    )
{
  double delta_num,scale_num,shift_num=0.0,roundoff;

  //if( singular &&( frac.denom > 1 || level_m != TUNE )) {

  if( fr->denom > 1 || level_m != TUNE ) {

    switch( my_random()% 64 ) { // allow lengthen or shorten

    case 0: // chop off final digit, or toggle dot/non-dot
      if(    fr->denom > 1
         &&( fr->denom > globals.base || level_m != TUNE )) {
        fr->denom /= globals.base; // chop off the last digit
        roundoff = (fr->num % globals.base)/(double)globals.base;
        fr->num /= globals.base;
        if( roundoff > 0.0 && random_uniform() < roundoff ) {
          fr->num++;
        }
      }
      else if ( fr->denom == 1 && level_m != TUNE ) {
        fr->has_dot = !fr->has_dot;
      }
      break;
#ifndef INTEGERS_ONLY
    case 1: // add an extra digit
      if( fr->denom < MAX_INT/globals.base && level_m != TRIM
                 &&( fr->denom > 1 || level_m != TUNE )) {
        // TUNE mutations cannot switch between integer & non-integer
        fr->has_dot = TRUE;
        fr->denom  *= globals.base;
        fr->num    *= globals.base;
      }
      break;
#endif
    }
  }
#ifdef DEBUG
  //di += sprintf( &ds[di],"pc=%d, fd=%d n=%ld,",
  //             precision_change,frac_digits,n);
#endif
  if( fr->whole == 0 && fr->num == 0 ) {
    delta_num = 0.5*scale*pow((double)fr->denom,0.25)* random_gaussian();
  }
  else if( level_m != TUNE && random_bit()) { // generate new value
    delta_num = fr->num *( -0.25 + random_gaussian());
  }
  else { // just tune the previous value by a small amount
    double numerator = ((fr->whole*(double)fr->denom))+(double)fr->num;
    double value = fr->whole + fr->num/(double)fr->denom;
    scale_num = 0.5*scale*pow(numerator*(double)fr->denom,0.25);
    if( value > 1.0 ) {
      shift_num = 0.1*(value + 1.0)/value;
    }
    delta_num = scale_num *(-shift_num + random_gaussian());
    //delta_num = (1.0 + log(1.0+fr->whole*fr->denom+fr->num))*random_gaussian();
    //delta_num = sqrt((double)fr->denom)*random_gaussian();
    //delta_num = 0.5 * globals.base * random_gaussian();
    //delta_num = pow((double)fr->denom,0.25)*random_gaussian();
  }

  adjust_fraction( fr,delta_num,singular );

#ifdef DEBUG
    //di += sprintf( &ds[di]," fd=%d n=%ld\n",frac_digits,n);
#endif
}

/********************************************************//**
   Assign weights w[c] to all cells from 1 to num_cells-1,
   according to the size of the smallest block common to
   c,c0 (considering only blocks of size Nc or larger)
*/
void focus_weights(int w[],int c0,int Nc,int num_cells)
{
  int block_size = Nc;
  int c;
  for( c=0; c < num_cells; c++ ) {
    w[c] = 1;
  }
  while( block_size < num_cells ) {
    c0 = ( c0 / block_size ) * block_size;
    for( c = c0; c < c0 + block_size; c++ ) {
      w[c] = w[c] + w[c];
    }
    block_size = block_size + block_size;
  }
}

/********************************************************//**
   Assume fr->whole and fr->num are non-negative.
   Add delta_num to fr->num (rounding probabilistically).
   Ensure fracion remains proper and non-negative
  (flipping fr->is_negative if necessary)
*/
void adjust_fraction(
                     Fraction *fr,
                     double delta_num,
                     int    force_change
                    )
{
  Long   delta_num_long,delta_whole;
  double delta_num_frac;

  //delta_num_long = floor( delta_num - 0.001*num ); ELASTIC??
  delta_num_long = floor( delta_num );
  delta_num_frac = delta_num - delta_num_long;
  if( delta_num_frac > 0.0 ) { // round probabilistically
    if( random_uniform() < delta_num_frac ) {
      delta_num_long++;
    }
  }
  if( delta_num_long == 0 && force_change ) {
    delta_num_long = random_bit() ? 1 : -1 ;
  }
  if( delta_num_long != 0 ) {
    // tp->value_changed = TRUE;
    fr->num += delta_num_long;
  }
  // adjust whole number, sign if necessary
  if( fr->num >= fr->denom ) { // ensure proper fraction
    if( fr->whole > MAX_INT - fr->num/fr->denom ) {
      fr->whole  = MAX_INT;
      fr->num    = fr->denom - 1;
    }
    else {
      fr->whole += fr->num / fr->denom;
      fr->num    = fr->num % fr->denom;
    }
  }
  else if( fr->num < 0 ) {
    if( (-fr->num-1)/fr->denom >= fr->whole ){// change sign
      fr->whole =  -fr->num/fr->denom - fr->whole;
      fr->num   = (-fr->num) % fr->denom;
      fr->is_negative = !fr->is_negative;
    }
    else {
      delta_whole = (fr->denom-1 - fr->num)/fr->denom;
      fr->whole -= delta_whole;
      fr->num   += delta_whole * fr->denom;
    }
  }
}

/********************************************************//**
   set Lb,Rb ready for copying (remainder of???) this cell
*/
void reset_cell_bars( Template *tp,Code *cd,Clip *cp )
{
  if( cp->c < cd->num_cells && tp->c < tp->num_cells) {
    cp->Lb            = cd->cell[cp->c];
    cp->Rb            = cd->cell[cp->c+1];
    tp->cell[tp->c+1] = tp->b + cp->Rb - cp->b;
    tp->Lb            = tp->cell[tp->c];
    tp->Rb            = tp->cell[tp->c+1];
  }
}

/********************************************************//**
   Store the digits encoding integer n to tp->codon[],
   starting at tp->codon[tp->k+1]
   Update tp->k to be the last index of the digits.
   Also update tp->ival[d,d+1] as appropriate
*/
void insert_integer( Template *tp, int n )
{
  int digit;
  int size=1;
  tp->d = tp->k+1;
  tp->ival[tp->d+1] = n;
  if( n == 0 ) {
    //if( random_bit()) { // ?????????????????????
    transcribe_codon(tp,'0');
  }
  else {
    while( size <= n / globals.base ) {
      size *= globals.base;
    }
    while( size > 0 ) {
      digit = n / size;
      transcribe_codon(tp,codon_from_digit(digit));
      n -= digit * size;
      size /= globals.base;
    }
  }
  tp->ival[tp->d] = tp->k+1 - tp->d;
}

/********************************************************//**
   Increment tp->k, Store a DOT '.' into tp->codon[tp->k],
   and set tp->ival[tp->k] equal to 1
*/
void insert_dot( Template *tp )
{
  transcribe_codon(tp,DOT);
  tp->ival[tp->k] = 1;
}

/********************************************************//**
   Increment tp->k and Store the specified operator to
   tp->codon[tp->k]. If it is not preceded by dot/digits,
   set tp->ival[tp->k] to zero.
*/
void insert_operator( Template *tp, int op )
{
  if( op == NEG ) {
    if(      tp->codon[tp->k] == PSH
       ||(   tp->codon[tp->k] == NEG
          && tp->codon[tp->k-1] == PSH )) {
      negate_push_value( tp );
    }
    else {
      transcribe_codon( tp,NEG );
      tp->ival[tp->k] = 0; // negation instruction
    }
  }
  else {
    transcribe_codon( tp,op );
    if( !is_dot_or_digit(tp->codon[tp->k-1])) {
      tp->ival[tp->k] = 0; // single-codon instruction
    }
  }
}

/********************************************************//**
    Combine #  with - to produce #-
 or combine #- with - to produce #
*/
void negate_push_value( Template *tp )
{
  int d = tp->k;
  if( tp->codon[tp->k] == PSH ) { // #  - to #-
    transcribe_codon( tp,NEG );
  }
  else {                          // #- - to #
    tp->k--;
    d--;
  }
  while( is_dot_or_digit(tp->codon[d-1])) {
    d--;
  }
  tp->ival[d] = tp->k - d;
  if( is_dot_or_digit(tp->codon[d])) {
    tp->fval[tp->ival[d+1]] *= -1.0;
  }
}

/********************************************************//**
   Store the digits encoding fraction whole + num/denom
   in globals.base, to tp->codon[], starting at tp->k+1.
   Set tp->d to be (current) tp->k+1,
       tp->k to be the index of the operator (PSH or NEG)
   Also update ival[],fval[],param[], tp->v,p
*/
void insert_push( Template *tp, Fraction frac )
{
  Floating value = val_from_fraction(&frac);
  Long size=1;
  int digit;
  tp->d = tp->k+1;
  while( size <= frac.whole / globals.base ) {
    size *= globals.base;
  }
  if( frac.whole > 0 ) {
    while( size > 0 ) {
      digit = frac.whole / size;
      transcribe_codon(tp,codon_from_digit(digit));
      frac.whole -= digit * size;
      size /= globals.base;
    }
  }
  if( frac.denom > 1 || frac.has_dot ) {
    transcribe_codon( tp,DOT );
    while( frac.denom > 1 ) {
      // assume denom is a power of base
      frac.denom /= globals.base;
      digit = frac.num / frac.denom;
      transcribe_codon(tp,codon_from_digit(digit));
      frac.num -= digit * frac.denom;
    }
  }
  transcribe_codon( tp, PSH );
  if( frac.is_negative ) {
    transcribe_codon( tp, NEG );
  }
  tp->ival[tp->d] = tp->k - tp->d;
  if( is_dot_or_digit(tp->codon[tp->d])) {
    tp->ival[tp->d+1] = tp->v;
    tp->fval[tp->v]   = value;
    next_value( tp );
  }
  if( frac.has_dot ) {
    tp->param[tp->p] = tp->d;
    next_param( tp );
  }
}

/********************************************************//**
   Insert a comparison to preceed a branch instruction
*/
void insert_comparison( Template *tp, int op )
{
  Fraction frac;
  int pair,query;
  int force_integer;
  if( op == BRF ) {
    pair  = random()% 5;
    query = random()% 6;
    if( pair == REG_VAL || pair == STACK_VAL ) {
      force_integer = ( query < 2 );
      frac = random_fraction( tp,force_integer,TRUE );
      insert_push( tp,frac );
    }
    switch( pair ) {
    case STACK_REG:   //  [x](r)[g=][x]
      if(random_bit()) insert_operator(tp,SWP);
      insert_register(tp,random_register(tp));
      insert_query(tp,query);
      if(random_bit()) insert_operator(tp,SWP);
      break;
    case REG_VAL:     //     (r)[g=][!]
      insert_register(tp,random_register(tp));
      insert_query(tp,query);
      if(random_bit()) insert_operator(tp,POP);
      break;
    case STACK_VAL:   //      .>[g=]
      insert_dot(tp);
      insert_operator(tp,PUT);
      insert_query(tp,query);
      break;
    case STACK_STACK: //      .>[g=][<]
      insert_dot(tp);
      insert_operator(tp,PUT);
      insert_query(tp,query);
      if(random_bit()) insert_operator(tp,GET);
      break;
    case REG_REG:     // (r)<(r)[g=][!]
      insert_register(tp,random_register(tp));
      insert_operator(tp,GET);
      insert_register(tp,random_register(tp));
      insert_query(tp,query);
      if(random_bit()) insert_operator(tp,POP);
    }
    if( random_bit()) {
      if( random_bit()) {
        insert_operator( tp, AND );
      }
      else {
        insert_operator( tp, OR );
      }
    }
  }
}

/********************************************************//**
   Return a random number between -2 and tp->num_reg-1,
   with a preference for smaller values. ???????????
   -1 indicates local register (.)
   -2 indicates an unspecified register
*/
int random_register( Template *tp )
{
  int r;
  if( tp->num_reg == 0 || random_bit()) {
    r = -2 + random_bit();     // local, or unspecified
  }
  else {
    r = random()% tp->num_reg; // all equally likely
  }
  return( r );
}

/********************************************************//**
   Store the digits encoding register r, in globals.base,
   to tp->codon[], starting at tp->codon[tp->k+1]
   Update tp->k to be the last index of the digits.
*/
void insert_register( Template *tp, int r )
{
  // if( r == -2 ) leave blank (unspecified register)
  if( r == -1 ) {      // local register
    insert_dot( tp );
  }
  else if( r >= 0 ) {
    insert_integer( tp, r );
  }
}

/********************************************************//**
   Insert a query
*/
void insert_query( Template *tp, int query )
{
  switch( query ) {
   case 0: //                 equal
     insert_operator( tp, EQL );
     break;
   case 1: //             not equal
     insert_operator( tp, EQL );
     insert_operator( tp, NOT );
     break;
   case 2: // greater than
     insert_operator( tp, GRT );
     break;
   case 3: // greater than or equal
     insert_operator( tp, GRT );
     insert_operator( tp, EQL );
     insert_operator( tp, OR  );
     break;
   case 4: //    less than or equal
     insert_operator( tp, GRT );
     insert_operator( tp, NOT );
     break;
   case 5: //    less than
     insert_operator( tp, GRT );
     insert_operator( tp, EQL );
     insert_operator( tp, OR  );
     insert_operator( tp, NOT );
     break;
  }
}

/********************************************************//**
   Insert a randomly generated instruction to
   tp->codon[], starting at tp->k+1.
   Update tp->k to be the index of the last
   codon in the instruction.
*/
void insert_instruction( Template *tp )
{
  Fraction frac;
  int offset;
  int op;
  int c,r;

  do {
    op = phen[random()% NUM_CODONS];
  } while(   is_dot_or_digit(op) || !is_allowed[op]
#ifndef CAN_INSERT_BRANCH
          || op == BRF || op == BRB
#endif
          || op == BLN || op == JSR );
  // ||(tp->num_cells == 1 && op == JSR));

  if( op_uses_digits[op] ) {

    switch( op ) {

    case GET: case PUT: case LOD: case STO:
    case INC: case DEC: case EQL: case GRT:
      r = random_register( tp );
      insert_register( tp,r );
      insert_operator( tp,op );
      break;

    case JSR:
      c = choose_jump_index( tp->c,tp->num_cells );
      insert_integer( tp, c );
      insert_operator( tp,JSR );
      break;

    case BRB: case BRF:
      if( op == BRF && random_bit()) {
        insert_comparison( tp,op );
      }
      if( random_bit()) {
        if( op == BRB ) { // ????????? -1 ????????????
          offset=random_offset(-1,tp->b - tp->cell[tp->c]);
        }
        else {
          offset=random_offset(-1,tp->cell[tp->c+1]-1 - tp->b);
        }
        if( offset > 0 ) {
          insert_integer( tp, offset );
        }
      }
      insert_operator( tp,op );
      break;

    case PSH:
      frac = random_fraction( tp,FALSE,TRUE );
      insert_push( tp,frac );
      break;
    }
  }
  else {
    insert_operator( tp,op );
  }
}



/********************************************************//**
   Generate a random fraction (push value) and return it.
*/
Fraction random_fraction(
                         Template *tp,
                         int force_integer,
                         int force_dot // ALWAYS TRUE ?????????
                        )
{
  Fraction frac;
  frac.num = 0;
  frac.denom = 1;
  frac.is_negative = ( random()% 3 == 0 );
#ifndef INTEGERS_ONLY
  if( !force_integer ) {
    while( frac.denom < MAX_INT/globals.base && random_bit()) {
      frac.denom *= globals.base;
    }
  }
#endif
  if( frac.denom == 1 ) { // integer
    frac.whole =   frac.is_negative
                + (Long)floor(pseudo_cauchy(globals.base));
    frac.has_dot = force_dot || random_bit();
  }
  else {
    frac.num = frac.is_negative + (Long)floor(pseudo_cauchy(
                         globals.base*sqrt((double)frac.denom)));
    frac.whole = frac.num / frac.denom;
    frac.num   = frac.num % frac.denom;
    frac.has_dot = TRUE;
  }
  return( frac );
}

/********************************************************//**
   Interpolate cd0->param[p] with cd1->param[p] and insert
   the result to tp. Also update tp->p,v,param[],value[]
*/
void interpolate_fraction(
                          Template *tp,
                          Code  *cd0,
                          Code  *cd1,
                          int    p,
                          double alpha
                         )
{
  Fraction frac0,frac1;
  double   delta_num,roundoff=0.0;

  scan_fraction( &frac0,cd0->codon,cd0->param[p] );
  scan_fraction( &frac1,cd1->codon,cd1->param[p] );

  if( frac1.is_negative != frac0.is_negative ) {
    frac1.is_negative = !frac1.is_negative;
    frac1.whole = -frac1.whole;
    frac1.num   = -frac1.num;
  }
  if( frac1.denom > frac0.denom ) {    // increase precision
    frac0.num  *= ( frac1.denom/frac0.denom );
    frac0.denom =   frac1.denom;
  }
  else if( frac1.denom < frac0.denom ){// decrease precision
    Long ratio = frac0.denom / frac1.denom;
    roundoff = ( frac0.num % ratio )/(double)(ratio);
    frac0.num   /= ratio;
    frac0.denom /= ratio;
  }

  delta_num = (frac1.whole - frac0.whole)*(double)frac0.denom
             + frac1.num  - (frac0.num + roundoff);

  adjust_fraction( &frac0,roundoff+alpha*delta_num,FALSE );

  insert_push( tp,frac0 );
}

/********************************************************//**
*/
void fix_bar_codons(
                    Template *tp,
                    char d_codon,
                    int  must_change
                   )
{
  int final_bar  = ( tp->b+1 == tp->cell[tp->c+1] );
  int final_cell = ( tp->c+1 == tp->num_cells );
  if(       final_bar
     &&(                  d_codon == DOT
        ||( final_cell && d_codon != BLN ))) {
    must_change = TRUE;
  }
  if( must_change ) {
    if( random_bit()) {   // insert dot or digit
      if( random_bit()) { // insert dot
        if( !final_bar ) {
          insert_dot( tp );
        }
      }
      else {              // insert digit
        if( !( final_bar && final_cell )) {
          insert_integer( tp,8 );
        }
      }
    }          //     otherwise, no dot or digit
  }
  else {       // keep any existing dot or digit
    if( d_codon == DOT ) {
      insert_dot( tp );
    }
    else if( d_codon != BLN ) {
      insert_integer( tp,8 );
    }
  }
  insert_operator( tp, BLN );
  next_bar( tp );
  tp->bar[tp->b] = tp->k;
}

/********************************************************//**
   Assume cp->d is the index of the first codon in the
   instruction, cp->k that of the operator.
   Fix the digits, if necessary, and copy them to
   tp->codon[], starting at tp->k+1. Update tp->k, cp->k
   to be the index of the last codon in the instruction.
*/
void fix_digits(
                Template *tp,
                Code *cd,
                Clip *cp,
                int must_change
               )
{
  Fraction frac;

  switch( cd->codon[cp->k] ) {

  case GET: case PUT: case LOD: case STO:
  case INC: case DEC: case EQL: case GRT:
    fix_register( tp,cd,cp,must_change );
    break;

  case BRF: case BRB:
    fix_offset( tp,cd,cp,must_change );
    break;

  case JSR:
    fix_jump_index(tp,cd,cp,must_change);
    break;

  case PSH: case NEG:
    scan_fraction( &frac,cd->codon,cp->d );
    if( must_change ) {
      mutate_fraction( &frac,tp->level.m,1.0,TRUE );
    }
    insert_push( tp,frac );
    break;

  case BLN:
    fix_bar_codons(tp,cd->codon[cp->d],must_change);
    cp->b++;   // move to next bar
    if( cp->b == cd->cell[cp->c+1] ) { // end of cell
      cp->c++;
    }
    if( tp->b == tp->cell[tp->c+1] ) { // end of cell
#ifdef DEBUG
      print_completed_cell( tp,tp->c );
#endif
      tp->c++;
      reset_cell_bars( tp,cd,cp ); // reset Lb,Rb,cell[c+1]
    }
    break;
  }
}

/********************************************************//**
   Assume cd->codon[cp->k] is either BRB or BRF.
   If the branch offset needs to be fixed, fix it.
   If must_change is TRUE, and offset has not already
   changed, mutate it.

   Lb and Rb are the left and right edge,
   respectively, of preserved section(s) of code.

   If Rb <= Lb, it is assumed that
   [0..Rb] and [Lb..end] are preserved,
   while the code in the middle is modified, i.e.
<pre>
   [0..gp->Rb][old middle][gp->Lb..end]
    |  | |  |               \   \  \  \
   [0..cp->Rb][new middle ][cp->Lb..end]
</pre>
   If Rb > Lb, it is assumed that
   [Lb..Rb] will be preserved, while the code
   to to the left and/or right is modified, i.e.
<pre>
   [old_beginning][cp->Lb..cp->Rb][old_ending]
                    \   \   \   \
   [new_beginning ][gp->Lb..gp->Rb][new_ending]
</pre>
*/
void fix_offset(
                Template *tp,
                Code *cd,
                Clip *cp,
                int must_change
               )
{
  int op = cd->codon[cp->k]; //  BRB or BRF
  int offset,max_offset; // old_offset
  int Bb;
  offset = int_from_codons( cd->codon,cp->d );
  //old_offset = offset;
  if( op == BRB ) {          //  branch back
    Bb = cp->b - offset;     //  branching to start of Bb
    if( cp->Lb < cp->Rb ) {  //  middle section is preserved
      if( Bb <= cp->Lb ) {   //       .....+-----+.....
                             //            L       R
        offset = tp->b - align_bar(tp->cell[tp->c],tp->Lb,
                                   cd->cell[cp->c],cp->Lb,Bb );
      }
    }              // otherwise, middle section is replaced                 
    else {                   //      +-----+.....+-----+
      if( Bb <= cp->Rb ) {   //            R     L
        offset += (tp->b - tp->Rb) - (cp->b - cp->Rb);
      }
      else if( Bb <= cp->Lb && cp->b >= cp->Lb ) {
                             //  branch into middle section
        offset = tp->b - align_bar(tp->Rb,tp->Lb,cp->Rb,cp->Lb,Bb);
      }
    }
  }
  else if( op == BRF ) {     //  branch forward
    Bb = cp->b+1 + offset;   //  branching to start of Bb
    if( cp->Lb < cp->Rb ) {  //  middle section is preserved
      if( Bb >= cp->Rb ) {   //       .....+-----+.....
                             //            L       R
#ifdef STRETCH_OFFSETS
        offset=align_bar(tp->Rb-1,tp->cell[tp->c+1],
                         cp->Rb-1,cd->cell[cp->c+1],Bb)-(tp->b+1);
#else
        offset=align_bar(tp->Rb,tp->cell[tp->c+1],
                         cp->Rb,cd->cell[cp->c+1],Bb)-(tp->b+1);
#endif
      }
    }              // otherwise, middle section is replaced
    else {                   //      +-----+.....+-----+
      if( Bb > cp->Lb ) {    //            R     L
        offset += (tp->Lb - tp->b) - (cp->Lb - cp->b);
      }
      else if( Bb > cp->Rb && cp->b <= cp->Rb ) {
                             // branch into middle section
        offset=align_bar(tp->Rb,tp->Lb,cp->Rb,cp->Lb,Bb)-(tp->b+1);
      }
    }
  }
#ifdef DEBUG
  //if( offset != old_offset ) {// && tp->level.m == CELL ) {
    //print_debug = TRUE;
  //}
#endif
  if( must_change ) {
    if( op == BRB ) {
      max_offset = tp->b - tp->cell[tp->c];
    }
    else { // ( op == BRF ) {
      max_offset = tp->cell[tp->c+1] - ( tp->b+1 );
    }
    // randomly mutate branch offset
    offset = random_offset( offset, max_offset );
  }
  if( offset > 0 ) {
    insert_integer( tp, offset );
  }
  insert_operator( tp,op );
}

/********************************************************//**
   Change to a new register with probability 1/2,
   or if the register index exceeds tp->num_reg
*/
void fix_register(
                  Template *tp,
                  Code *cd,
                  Clip *cp,
                  int must_change
                 )
{
  int r;
  if( !must_change ) {
    if( codon_is_digit( cd->codon[cp->d])) {
      r = int_from_codons( cd->codon,cp->d );
      if( r < tp->num_reg ) {
        insert_register( tp, r );
      }
      else {
        must_change = TRUE;
      }
    }
    else if( cd->codon[cp->d] == '.' ) {
      insert_dot( tp );
    }
  }
  if( must_change ) {
    r = random_register( tp );
    insert_register( tp, r );
  }
  insert_operator( tp,cd->codon[cp->k] );
}

/********************************************************//**
   Return TRUE if the jump index has changed,
   FALSE otherwise.
*/
void fix_jump_index(
                    Template *tp,
                    Code *cd,
                    Clip *cp,
                    int must_change
                   )
{
  int c_relative;
  int c=0,r;
  if( !must_change ) {
    //if( codon_is_digit( cd->codon[cp->d])) {
    c = int_from_codons( cd->codon,cp->d ); // absolute address
    if( tp->c != cp->c ) {
      r=0;
      // find largest "virtual block"
      while(( tp->c - cp->c )% power2[r+1] == 0 ) {
        r++;
      }
      // preserve relativity of calls within virtual block
      if( c/power2[r] == cp->c/power2[r] ) {
        c_relative = c + tp->c - cp->c; // relative address
        if( c_relative < tp->num_cells ) {
          c = c_relative;
        }
      }
      else if( c >= tp->c ) {
        must_change = TRUE;
      }
    }
    if( c >= tp->num_cells ) {
      must_change = TRUE;
    }
  }
  if( must_change ) { // choose cell randomly
    c = choose_jump_index( tp->c,tp->num_cells );
  }
  //if( c != tp->c || random_bit()) {
  //if( c != 0 ) {
  insert_integer( tp,c );
    //}
  insert_operator( tp,JSR );
}

/********************************************************//**
*/
int align_bar( int b1_min, int b1_max,
               int b0_min, int b0_max, int b0 )
{
  int b;
  if( b1_max < b1_min ) {
    b      = b1_min;
    b1_min = b1_max;
    b1_max = b;
  }
  b = b1_min + (int) floor(
       ( ( b0 - b0_min + random_uniform())
        *(b1_max+1 - b1_min)/(double)(b0_max+1 - b0_min)));
  return( b );
}

#ifdef DEBUG
/********************************************************//**

*/
void print_spaces( int num_spaces )
{
  int j;
  for( j=0; j < num_spaces; j++ ) {
    di += sprintf(&ds[di]," ");
  }
}

/********************************************************//**
*/
void print_line()
{
  di += sprintf( &ds[di],
 "____________________________________________\n");
  di += sprintf( &ds[di],
 "                             ");
}

/********************************************************//**
*/
void print_parent_cell( Code *cd, int c )
{
  int k;
  if( c >= cd->Fc && c < cd->Fc + cd->Nc ) {
    di += sprintf(&ds[di],"f");
  }
  else {
    di += sprintf(&ds[di]," ");
  }
  di += sprintf(&ds[di],"  [");
  k = cd->bar[cd->cell[c]]+1;
  while( k < cd->bar[cd->cell[c+1]] ) {
    di += sprintf(&ds[di],"%c",cd->codon[k++]);
  }
  di += sprintf(&ds[di],"]\n");
}

/********************************************************//**
*/
void print_clip( Code *cd, Clip *cp )
{
  int k;
  for( k = cp->k+1; k <= cp->k_stop; k++ ) {
    di += sprintf(&ds[di],"%c",cd->codon[k]);
  }
  di += sprintf(&ds[di],"\n");
}

/********************************************************//**
  Print the current cell of the PRIMARY parent (cd), up to the
  mutation point, then print the code that is being mutated.
*/
int print_clip_template(
                        Code *cd,
                        Clip *cp,
                        Template *tp
                       )
{
  int gap=0;
  int ch;
  int k,n;

  if( cp->c >= cd->Fc && cp->c < cd->Fc + cd->Nc ) {
    di += sprintf(&ds[di],"f"); // focus cell
  }
  else {
    di += sprintf(&ds[di]," ");
  }
  ch = '[';
  for( n=0; n < cd->num_focus_bars; n++ ) {
    if( cd->cell[cp->c] == cd->Fb[n].b ) {
      ch = 'F';
    }
  }
  di += sprintf(&ds[di],"  %c",ch);
  if( tp->level.m < BAR || cp->k_stop < cp->k ) {
    // align cp->k_stop with tp->k
    gap = (cp->k_stop - cd->bar[cd->cell[cp->c]])
         -(tp->k      - tp->bar[tp->cell[tp->c]]);
    print_spaces( -gap );
  }
  k = cd->bar[cd->cell[cp->c]]+1; // first codon in cell
  if( k == cp->k+1 ) {
    di += sprintf(&ds[di]," ");
  }
  if( k == cp->k_stop+1 && tp->level.m >= BAR ) {
    di += sprintf(&ds[di]," ");
  }
  while( k < cd->bar[cd->cell[cp->c+1]] ) {
    di += sprintf(&ds[di],"%c",cd->codon[k++]);
    if( k == cp->k+1 ) {
      di += sprintf(&ds[di]," ");
    }
    if( k == cp->k_stop+1 && tp->level.m >= BAR ) {
      di += sprintf(&ds[di]," ");
    }
  }
  di += sprintf(&ds[di],"]\n");
  di += sprintf(&ds[di],"   [");
  if( tp->level.m < BAR || cp->k_stop < cp->k ) {
    print_spaces( gap );
  }
  for( k=tp->bar[tp->cell[tp->c]]+1; k <= tp->k; k++ ) {
    di += sprintf(&ds[di],"%c",tp->codon[k]);
  }
  if( gap > 0 ) {
    return( 4+gap );
  }
  else {
    return( 4 );
  }
}

/********************************************************//**
*/
void print_template( Template *tp )
{
  int k;
  for( k=tp->bar[tp->cell[tp->c]]+1; k <= tp->k; k++ ) {
    di += sprintf(&ds[di],"%c",tp->codon[k]);
  }
  di += sprintf(&ds[di],"\n");
}

/********************************************************//**
*/
void print_completed_cell( Template *tp, int c )
{
  int ch='[';
  int k,n;
  if( c >= tp->Fc && c < tp->Fc + tp->Nc ) {
    di += sprintf(&ds[di],"f");
  }
  else {
    di += sprintf(&ds[di]," ");
  }
  for( n=0; n < tp->num_focus_bars; n++ ) {
    if( tp->cell[c] == tp->Fb[n].b ) {
      ch = 'F';
    }
  }
  di += sprintf(&ds[di],"%2d%c",c,ch);
  if( tp->overflow ) {
    di += sprintf(&ds[di]," OVERFLOW ]\n");
  }
  else {
    k = tp->bar[tp->cell[c]]+1;
    while( k < tp->bar[tp->cell[c+1]] ) {
      ch = tp->codon[k];
      if( ch == '|' ) {
        for( n=0; n < tp->num_focus_bars; n++ ) {
          if( tp->bar[tp->Fb[n].b] == k ){
            ch = 'F';
          }
        }
      }
      di += sprintf(&ds[di],"%c",ch);
      k++;
    }
    di += sprintf(&ds[di],"]\n");
  }
}

/********************************************************//**
*/
void flush_debug()
{
  print_debug = TRUE;
  if( print_debug ) {
    di += sprintf(&ds[di],"%c",'\0');
    printf("%s",ds);
    print_debug = FALSE;
  }
  di=0;
}
#endif


/********************************************************//**
   Assume min is either 0 or 1.
   Return random variable from a Cauchy distribution
   (starting from min)

double random_cauchy( int min )
{
  double x = random_uniform();
  if( min == 1 ) {
    x = 0.5 *( x + 1.0 );
  }
  return( tan( x * 3.14159265359/2.0));
}
*/

