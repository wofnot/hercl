/** \file cross_mutate.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include "step.h"
#include "scan_print.h"
#include "point.h"
#include "cross_mutate.h"

#define  PRIMARY   0
#define  SECONDARY 1

#define  MAX_MUTATIONS  20
#define  ALIGNED_CODON   1

void try_cross_mutate(Template *tp,Code *cd0,Code *cd1);
void       copy_focus(Template *tp,Code *cd);
void        copy_code(Template *tp,Code *cd0);
void    interpolation(Template *tp,Code *cd0,Code *cd1);
void    tune_mutation(Template *tp,Code *cd);
void    trim_mutation(Template *tp,Code *cd);
void   point_mutation(Template *tp,Code *cd);
void    bar_crossover(Template *tp,Code *cd0,Code *cd1);
void branch_crossover(Template *tp,Code *cd0,Code *cd1);
void   cell_crossover(Template *tp,Code *cd0,Code *cd1);
//void   jump_crossover(Template *tp,Code *cd0,Code *cd1);
void  block_crossover(Template *tp,Code *cd0,Code *cd1,int block_size);

void      cross_cells(Template *tp,Code *cd0,Clip *cp0,int c,Code *cd1);
void    insert_branch(Template *tp,Code *cd0,Clip *cp0,Code *cd1);
void      insert_jump(Template *tp,Code *cd0,Clip *cp0);
void bar_line_removal(Template *tp,Code *cd);

int    check_template(Template *tp);

void        clip_cell(Template *tp,Code *cd,Clip *cp,int c);
void    copy_verbatim(Template *tp,Code *cd,Clip *cp);
void  copy_fix_digits(Template *tp,Code *cd,Clip *cp);

void    add_focus_bar(Template *tp,int b);

int          adjust_focus(Code *cd,int Fc[],int Nc[]);
int          choose_codon(Code *cd,int Fc[],int Nc[],int r_max,int branch_or_jump);
int choose_block_by_focus(Code *cd,int block_size);
int   choose_bar_by_focus(Code *cd);
int      choose_mutations(Code *cd,int level_m,int locus[20],int type[20]);

Boolean     code_is_empty(Code *cd);
Boolean  identical_blocks(Code *cd0,Code *cd1,int c,int n);
Boolean    identical_bars(Code *cd0,Code *cd1,int b0,int c);
Boolean   identical_clips(Code *cd0,Clip *cp0,Code *cd1,Clip *cp1,
                          int bar_or_cell,int portion);

void       random_clip(Code *cd,Clip *cp,int bar_or_cell,int portion,int parent);
void exclude_final_bar(Code *cd,Clip *cp,int bar_or_cell,int portion,int parent);
void     matching_clip(Code *cd0,Clip *cp0,Code *cd1,Clip *cp1,
                       int bar_or_cell,int portion);
void          cue_clip(Code *cd,Clip *cp,int ki,int kf,
                       int bar_or_cell,int portion,int parent);
int      aligned_codon(Code *cd0,int b0,int k0,Code *cd1,int b1);
int end_of_instruction(Code *cd,int k);
int  last_codon_in_bar(Code *cd,int b);
int last_codon_in_cell(Code *cd,int c);
int       cell_for_bar(Code *cd,int b);
int    get_param_index(Code *cd,int k);

int  choose_index(int max);
int similar_block(int block0,int num_blocks);


extern long long ncomp;


/********************************************************//**
   Produce new code, using cd0 as primary parent,
   cd1 as secondary parent and mlevel as mutation level.
*/
Code * cross_mutate(
                    Code *cd0,
                    Code *cd1,
                    Level level
                   )
{
  Template temp;

  temp.level      = level;
  temp.num_cells  = cd0->num_cells;
  temp.num_reg    = cd0->num_reg;
  temp.stack_size = cd0->stack_size;
  temp.mem_size   = cd0->mem_size;

  if( cd1 == NULL ) {
    cd1 = cd0;
  }
  if( code_is_empty(cd0) && code_is_empty(cd1)) {
    temp.level.m = POINT; // both empty code
    temp.level.g = NON_ALIGNED;
  }
  else if(        temp.level.g == ALIGNED
          &&(   cd0->num_cells != cd1->num_cells
             || strcmp(cd0->codon,cd1->codon) == 0 )) {
    temp.level.g = NON_ALIGNED;
  }
  do {
    try_cross_mutate( &temp,cd0,cd1 );
  } while( temp.overflow );

  if( level.m == COPY ) { // || mlevel == STRIP || TRIM ) {
    // treat new code as replacement for cd0
    temp.level = cd0->level;
  }

  return( code_from_template(&temp));
}

/********************************************************//**
   Try to produce new code, using cd0 as primary parent,
   cd1 as secondary parent and tp->level.m as mutation level.
   If new code exceeds certain limits, set tp->overflow to TRUE.
*/
void try_cross_mutate(
                      Template *tp,
                      Code     *cd0,
                      Code     *cd1
                     )
{
  int k;
  reset_template( tp );
  copy_focus( tp,cd0 );

  switch( tp->level.m ) {
   case COPY:         copy_code(tp,cd0);       break;
     //case GRAD:  gradient_descent(tp,cd0);   break;
   case INTERP:   interpolation( tp,cd0,cd1);  break;
   case TUNE:     tune_mutation( tp,cd0);      break;
   case TRIM:     trim_mutation( tp,cd0);      break;
   case POINT:   point_mutation( tp,cd0);      break;
   case BAR:       bar_crossover(tp,cd0,cd1);  break;
   case BRANCH: branch_crossover(tp,cd0,cd1);  break;
   case CELL:     cell_crossover(tp,cd0,cd1);  break;
     // case JUMP:     jump_crossover(tp,cd0,cd1);  break;
   default:      block_crossover(tp,cd0,cd1,
                  power2[tp->level.m-BLOCK]);  break;
  }
#ifdef DEBUG
  //if( strcmp(cd0->codon,tp->codon) == 0 ) {
  //  printf("%d\n",tp->level.m);
  //}
#endif
  k = check_template( tp );
  if( k > 0 ) {
    printf("mlevel: %d\n",tp->level.m);
    tp->codon[k+1] = '\0';
    printf("[%s]",tp->codon);
    print_code( cd0,stdout );
    print_code( cd1,stdout );
    printf("Something is wrong with the evolved code.\n");
    exit(1);
  }
}

/********************************************************//**
   Let code inherit focus cell(s) and bar(s) from template.
*/
void copy_focus( Template *tp, Code *cd )
{
  int m,n;
  tp->Fc = cd->Fc;
  tp->Nc = cd->Nc;
  tp->num_focus_bars = 0;
  for( n=0; n < cd->num_focus_bars; n++ ) {
    m=0;
    while( m < tp->num_focus_bars && tp->Fb[m].b != cd->Fb[m].b ) {
      m++;
    }
    if( m == tp->num_focus_bars ) {
      tp->Fb[m] = cd->Fb[n];
      tp->num_focus_bars++;
    }
  }
}

/********************************************************//**
   Generate random code with specified number of cells
*/
Code *random_code( Code *cd0, int num_cells )
{
  Template temp;
  Boolean is_allowed_BRB;
  int num_bars,num_instructions;
  int b,c,i;

  temp.num_cells = num_cells;
  if( cd0 == NULL ) {
    temp.level.m = POINT;
    temp.level.g = NON_ALIGNED;
    temp.num_reg = 10;
    temp.stack_size = 256;
    temp.mem_size = 1024;
  }
  else {
    temp.level      = cd0->level;
    temp.num_reg    = cd0->num_reg;
    temp.stack_size = cd0->stack_size;
    temp.mem_size   = cd0->mem_size;
  }

  reset_template(&temp);

  is_allowed_BRB = is_allowed[(int)BRB];
  // don't allow loops in random code
  is_allowed[(int)BRB] = FALSE;
  for( c=0; c < num_cells; c++ ) {
    num_bars = random_geometric(2,6);
    temp.cell[c+1] = temp.cell[c] + num_bars;
    for( b=0; b < num_bars; b++ ) {
      num_instructions = -1 + random_geometric(6,16)
                            + random_geometric(6,16);
      for( i=0; i < num_instructions; i++ ) {
        insert_instruction( &temp );
      }
      insert_operator( &temp,BLN );
      next_bar( &temp );
      temp.bar[temp.b] = temp.k;
    }
  }
  is_allowed[(int)BRB] = is_allowed_BRB;
  return( code_from_template(&temp));
}

/********************************************************//**
   Make an exact copy of the specified code.
*/
void copy_code( Template *tp, Code *cd0 )
{
  Clip clip;
#ifdef DEBUG
  print_line();
  di += sprintf( &ds[di],"     COPY CODE\n");
#endif
  //copy_focus( tp,cd0 );
  clip_cell( tp,cd0,&clip,0 );
  copy_verbatim( tp,cd0,&clip );
#ifdef DEBUG
  flush_debug();
#endif
}

/********************************************************//**
   Interpolate parameters between cd0 and cd1.
*/
void interpolation(
                   Template *tp,
                   Code *cd0,
                   Code *cd1
                  )
{
  Clip clip0;
  int p;
  clip_cell(tp,cd0,&clip0,0);
  for( p=0; p < cd0->num_param; p++ ) {
    clip0.k_stop = cd0->param[p]-1;
    copy_verbatim(tp,cd0,&clip0);
    interpolate_fraction(tp,cd0,cd1,p,globals.interp_rate);
    clip0.d = clip0.k+1;
    clip0.k = clip0.d + cd0->ival[clip0.d];
  }
  clip0.k_stop = cd0->last_codon;
  copy_verbatim(tp,cd0,&clip0);
#ifdef DEBUG
  flush_debug();
#endif
}

/********************************************************//**
   Mutate one or more push values.
*/
void tune_mutation( Template *tp, Code *cd )
{
  Fraction frac;
  Clip clip;
  int param[MAX_PUSH]; // which params to mutate
  int Fc[MAX_LEVEL];   // start cell (telescoping focus blocks)
  int Nc[MAX_LEVEL];   // number of cells
  int Fcp[MAX_LEVEL];  // start param in focus cell(s)
  int Ncp[MAX_LEVEL];  //  num params in focus cell(s)
  int Fbp[5];          // start param in focus bar(s)
  int Nbp[5];          //  num params in focus bar(s)
  double scale;        // size of changes to params
  double weight=1.0;   // weight of focus cell(s)
  double Np_total=0.0; // weighted total number of params
  int Np_bar=0;        // number of params in focus bar(s)
  int singular;        // whether there is only one mutation
  int num_mutations;
  int num_focus_bars=0;
  int r_min=0,r_max;
  int mutation=0;
  int i,j=0,m,r,p=0;
#ifdef DEBUG
  int gap=0;
  print_line();
  di += sprintf( &ds[di],"  TUNE MUTATION\n");
#endif
  clip_cell( tp,cd,&clip,0 );
  if( cd->num_param > 0 ) {
    r_max = adjust_focus( cd,Fc,Nc );
    for( r = 0; r <= r_max; r++ ) {
      Fcp[r]=get_param_index(cd,cd->bar[cd->cell[Fc[r]]]);
      Ncp[r]=get_param_index(cd,cd->bar[cd->cell[Fc[r]+Nc[r]]])-Fcp[r];
    }
    while( Ncp[r_min] == 0 ) {
      r_min++; // smallest index r for which Ncp[r] > 0
    }
    if( cd->num_focus_bars > 0 ) {
      // find the focus bars which actually contain params,
      // and count the number of params in each
      for( i=0; i < cd->num_focus_bars; i++ ) {
        Fbp[j]=get_param_index(cd,cd->bar[cd->Fb[i].b]);
        Nbp[j]=get_param_index(cd,cd->bar[cd->Fb[i].b+1])-Fbp[j];
        if( Nbp[j] > 0 ) {
          Np_bar += Nbp[j];
          j++;
        }
      }
      num_focus_bars = j;
    }
    if( Np_bar > 0 ) {
      weight  *= 0.5;
      Np_total = weight * Np_bar;
    }
    for( r = r_min; r < r_max; r++ ) {
      weight *= 0.5;
      Np_total += weight * Ncp[r];
    }
    Np_total += weight * Ncp[r_max];

    num_mutations = random_geometric(sqrt(Np_total),cd->num_param);

    // choose which params to mutate, and store in param[]
    while( mutation < num_mutations ) {
      // choose param p, from focus bar(s) or cell(s)
      if( num_focus_bars > 0 && random_bit()) {
        j = random()% num_focus_bars;
        p = Fbp[j] + random()% Nbp[j];
      }
      else {
        r = r_min + choose_index( r_max - r_min );
        p = Fcp[r] + random()% Ncp[r];
      }
      m = mutation;
      while( m > 0 && param[m-1] > p ) {// keep them sorted
        m--;
      }
      if( m == 0 || param[m-1] != p ) { // exclude duplicates
        m = mutation;
        while( m > 0 && param[m-1] > p ) {
          param[m] = param[m-1];
          m--;
        }
        param[m] = p;
        mutation++;
      }
    }
    scale = pow((double)num_mutations,-0.25);
    singular = ( num_mutations == 1 );
    for( m=0; m < num_mutations; m++ ) {
      clip.k_stop = cd->param[param[m]]-1;
      copy_verbatim(tp,cd,&clip);
#ifdef DEBUG
      gap = print_clip_template( cd,&clip,tp );
      di += sprintf(&ds[di],"\n");
#endif
      tp->d  = tp->k+1; // not needed ???????
      clip.d = clip.k+1;
      scan_fraction( &frac,cd->codon,clip.d );
      clip.k = clip.d + cd->ival[clip.d];
      mutate_fraction(&frac,TUNE,scale,singular);
      insert_push( tp,frac );
#ifdef DEBUG
      print_spaces( gap );
      print_template( tp );
#endif
    }
  }
  clip.k_stop = cd->last_codon;
  copy_verbatim(tp,cd,&clip);
  if( tp->overflow ) {
    reset_template( tp ); // just copy with no changes
    copy_focus( tp,cd );
    copy_code( tp,cd );
  }
#ifdef DEBUG
  //print_debug = TRUE;
  if(!strcmp(cd->codon,tp->codon)) {
    //print_debug = TRUE;
  }
  flush_debug();
#endif
}

/********************************************************//**
   Trim the code with a (restricted) point mutation,
   or removal of a random sub-bar.
*/
void trim_mutation( Template *tp, Code *cd )
{
  Clip clip;
  int portion; // FRONT, MIDDLE, BACK or WHOLE
#ifdef DEBUG
  int gap;
#endif
  int c,b,k,m,n;
  if( random()% 4 != 0 ) {
    //if( random()% 4 >= 0 ) {
    point_mutation( tp,cd ); // cannot insert
    return;
  }
#ifdef DEBUG
  print_line();
  di += sprintf( &ds[di]," SUB-BAR REMOVAL\n");
#endif
  // randomly choose bar b (and cell c) from cd
  //b = random()%(cd->cell[cd->num_cells]);
  b = choose_bar_by_focus( cd );
  c = cell_for_bar( cd,b );

  clip_cell(tp,cd,&clip,0);

  if( b+1 < cd->cell[c+1] ) { // not final bar
    if( end_of_instruction(cd,cd->bar[b]+1) == cd->bar[b+1]) {
      portion = WHOLE;
    }
    else {
      portion =   random()% 4; // WHOLE, FRONT, MIDDLE or BACK
    }
  }
  else {
    portion = 1 + random()% 3; // can't remove WHOLE of final bar
  }
  if( portion == WHOLE ) { // remove WHOLE bar
    // copy from beginning to start of cell
    clip.k_stop = cd->bar[cd->cell[c]];
    copy_verbatim(tp,cd,&clip);
    clip.Rb = b;       //    ---|-----|---
    clip.Lb = b;       //        Rb=Lb
    tp->Lb  = b-1;     //     -----|-----
    tp->Rb  = b;       //       Lb   Rb
    tp->cell[c+1]--;
    // copy from start of cell to start of bar
    clip.k_stop = cd->bar[b];
    copy_fix_digits(tp,cd,&clip);

    // remove deleted bar and decrement later focus bar(s)
    m=0;
    for( n=0; n < tp->num_focus_bars; n++ ) {
      if( tp->Fb[n].b < tp->b ) {
        tp->Fb[m++] = tp->Fb[n];
#ifdef DEBUG
        printf("bar %d kept\n",tp->Fb[n].b);
#endif
      }
      else if( tp->Fb[n].b > b ) {
        tp->Fb[m] = tp->Fb[n];
#ifdef DEBUG
        printf("bar %d -> %d\n",tp->Fb[m].b,tp->Fb[m].b-1);
#endif
        tp->Fb[m++].b--;
      }
      // otherwise, focus bar is removed
#ifdef DEBUG
      else {
        printf("bar %d removed\n",tp->Fb[n].b);
      }
#endif
    }
    tp->num_focus_bars = m;

    clip.b++;          // skip to next bar
    clip.k = cd->bar[clip.b];
#ifdef DEBUG
    gap = print_clip_template( cd,&clip,tp );
    di += sprintf(&ds[di],"\n");
#endif
    clip.k_stop = cd->bar[cd->cell[c+1]];
    copy_fix_digits(tp,cd,&clip); // copy to end of cell
  }
  else { // FRONT, MIDDLE or BACK
    // copy from beginning to start of bar
    clip.k_stop = cd->bar[b];
    copy_verbatim(tp,cd,&clip);
    // choose random subbar from cd
    random_clip(cd,&clip,BAR,portion,0);
    exclude_final_bar(cd,&clip,BAR,portion,0);
    // copy from start of bar to start of subbar
    k      = clip.k;          // store for later use
    clip.k = cd->bar[clip.b]; // start of bar
    copy_verbatim(tp,cd,&clip);
    clip.k = k;               // prepare to copy from end of subbar
#ifdef DEBUG
    gap = print_clip_template( cd,&clip,tp );
    di += sprintf(&ds[di],"\n");
    print_spaces( gap );
    print_template( tp );
#endif
  }
  clip.k_stop = cd->last_codon; // copy to end of code
  copy_verbatim(tp,cd,&clip);
#ifdef DEBUG
  if(!strcmp(cd->codon,tp->codon)) {
    //print_debug = TRUE;
  }
  flush_debug();
#endif
}

/********************************************************//**
   Choose one or more points to remove, insert or replace an
   instruction, or modify the dot/digits of an instruction.
*/
void point_mutation( Template *tp, Code *cd )
{
  Clip clip;
  int locus[MAX_MUTATIONS];
  int  type[MAX_MUTATIONS];
  int num_mutations;
  int m;
#ifdef DEBUG
  int gap=0;
#endif
  if(   cd->cell[cd->num_cells] > cd->num_cells
     && cd->level.m >= CELL && random()% 8 == 0 ) {
    bar_line_removal( tp,cd );
    return;
  }
#ifdef DEBUG
  print_line();
  if( tp->level.m == TRIM ) {
    di += sprintf( &ds[di],"  TRIM MUTATION\n");
  }
  else {
    di += sprintf( &ds[di]," POINT MUTATION\n");
  }
#endif
  num_mutations = choose_mutations(cd,tp->level.m,locus,type);

  clip_cell(tp,cd,&clip,0); // prepare to start copying

  for( m=0; m < num_mutations; m++ ) {
    // copy up to end of previous instruction
    clip.k_stop = locus[m]-1;
    while( is_dot_or_digit(cd->codon[clip.k_stop])) {
      clip.k_stop--;
    }
    copy_verbatim(tp,cd,&clip);
#ifdef DEBUG
    if( type[m] != NONE ) {
      gap = print_clip_template( cd,&clip,tp );
      if( type[m] == INSERT ) {
        di += sprintf(&ds[di],"       <insert>\n");
      }
      else if( type[m] == DIGIT ) {
        di += sprintf(&ds[di],"       <digit>\n");
      }
      else if( type[m] == EXCHANGE ) {
        di += sprintf(&ds[di],"       <replace>\n");
      }
      else {
        di += sprintf(&ds[di],"       <delete>\n");
      }
    }
#endif
    if( type[m] == INSERT || type[m] == EXCHANGE ) {
      insert_instruction( tp );
    }
    if( type[m] == REMOVE || type[m] == EXCHANGE ) {
      // skip to end of (next) instruction
      clip.k = end_of_instruction( cd,clip.k+1 ); // ival[] ????
    }
    else if( type[m] == DIGIT ) {
      //mutate_digits( tp,cd,&clip );
      clip.d = clip.k+1;
      clip.k = clip.d + cd->ival[clip.d];
      fix_digits( tp,cd,&clip,TRUE );
#ifdef DEBUG
  /*
  print_clip_template( cd,cp,tp );
  print_debug = TRUE;
  if( print_debug ) {
    printf("%s",ds);
    print_debug = FALSE;
  }
  di = 0;
  */
#endif
    }
#ifdef DEBUG
    if( type[m] != NONE ) {
      print_spaces( gap );
      print_template( tp );
    }
#endif
  }
  clip.k_stop = cd->last_codon;
  copy_verbatim(tp,cd,&clip); // copy to end of code
#ifdef DEBUG
  print_template( tp );
  if(!strcmp(cd->codon,tp->codon)) {
    //print_debug = TRUE;
  }
  flush_debug();
#endif
}

/********************************************************//**
   Replace the front, middle, back, fringe or whole of a
   randomly chosen bar in cd0 with the front, middle, back,
   fringe or whole of a bar from cd1.
*/
void bar_crossover(
                   Template *tp,
                   Code *cd0,
                   Code *cd1
                  )
{
  Clip clip0,clip1;
  int portion; // FRONT, BACK, MIDDLE, FRINGE or WHOLE
  int b0,c0;
  int k;
#ifdef DEBUG
  int gap;
#endif
#ifdef DEBUG
  print_line();
  di += sprintf( &ds[di],"  BAR CROSSOVER\n");
#endif

  b0 = choose_bar_by_focus( cd0 );
  c0 = cell_for_bar( cd0,b0 );

  if( tp->level.g == ALIGNED ) {
    while( identical_bars(cd0,cd1,b0,c0)) {
      b0 = choose_bar_by_focus( cd0 );
      c0 = cell_for_bar( cd0,b0 );
    }
    if(   cd1->cell[c0+1] - cd1->cell[c0]
       != cd0->cell[c0+1] - cd0->cell[c0] ) {
      tp->level.g = NON_ALIGNED;
    }
  }

  // c0 will become the focus cell for tp
  tp->Fc = c0;
  tp->Nc = 1;
  add_focus_bar( tp,b0 );

  // copy from cd0, up to start of bar
  clip_cell(tp,cd0,&clip0,0);
  clip0.k_stop = cd0->bar[b0];
  copy_verbatim(tp,cd0,&clip0);

  if( tp->level.g == ALIGNED ) {
    do {
      portion = random()% 5; // WHOLE, FRONT, MIDDLE, BACK or FRINGE
      random_clip(cd0,&clip0,BAR,portion,0);
      clip1.c = c0;
      matching_clip(cd0,&clip0,cd1,&clip1,BAR,portion);
    } while(identical_clips(cd0,&clip0,cd1,&clip1,BAR,portion));
  }
  else {
    portion = random()% 5; // WHOLE, FRONT, MIDDLE, BACK or FRINGE
    random_clip(cd0,&clip0,BAR,portion,0);
    clip1.b = random()% cd1->cell[cd1->num_cells];
    clip1.c = cell_for_bar( cd1,clip1.b );
    random_clip( cd1,&clip1,BAR,portion,1 );
  }

  if( portion == FRINGE ) { // replace front and back of bar
#ifdef DEBUG
    print_clip_template( cd0,&clip0,tp );
    di += sprintf(&ds[di],"\n");
    print_spaces(3+cd0->bar[clip0.b]-cd0->bar[cd0->cell[clip0.c]]);
    for( k = cd1->bar[clip1.b]; k <= cd1->bar[clip1.b+1]; k++ ) {
      di += sprintf(&ds[di],"%c",cd1->codon[k]);
      if( k == clip1.k ) {
        di += sprintf(&ds[di]," ");
      }
      if( k == clip1.k_stop ) {
        di += sprintf(&ds[di]," ");
      }
    }
    di += sprintf(&ds[di],"\n");
#endif
    // copy the front of the bar from cd1
    k       = clip1.k;           // store for later use
    clip1.k = cd1->bar[clip1.b]; // start of bar
    tp->Lb  = tp->b;             // preserve branch offsets only
    tp->Rb  = tp->b+1;           // within the bar itself
    copy_fix_digits(tp,cd1,&clip1);
    clip1.k = k;

    // copy preserved middle portion of bar from cd0
    tp->Lb  = tp->cell[tp->c];   // preserve all branch offsets
    tp->Rb  = tp->cell[tp->c+1];
    copy_fix_digits(tp,cd0,&clip0);

    // copy the back of the bar from cd1
    clip1.k_stop = cd1->bar[clip1.b+1]-1;
    if(is_dot_or_digit(cd1->codon[clip1.k_stop])) {
      clip1.k_stop--;            // bar ends with .| or 8|
    }
    tp->Lb  = tp->b;             // preserve branch offsets only
    tp->Rb  = tp->b+1;           // within the bar itself
    copy_fix_digits(tp,cd1,&clip1);
    fix_bar_codons(tp,cd1->codon[clip1.k_stop+1],FALSE);
    // don't need Lb,Rb because we will copy verbatim ??????????
    clip0.b++; // move to next bar
    clip0.k = cd0->bar[clip0.b];
    if( clip0.b == cd0->cell[clip0.c+1] ) { // end of cell
      clip0.c++;
#ifdef DEBUG
      print_completed_cell( tp,tp->c );
#endif
      tp->c++;
      reset_cell_bars(tp,cd0,&clip0);// is this needed ??????
    }
  }
  else { // non-FRINGE
    // copy from start of bar to start of subbar
    k       = clip0.k; // store for later use
    clip0.k = cd0->bar[clip0.b]; // start of bar
    copy_verbatim(tp,cd0,&clip0);
    clip0.k = k;       // prepare to copy from end of subbar
#ifdef DEBUG
    gap = print_clip_template( cd0,&clip0,tp );
    di += sprintf(&ds[di],"\n");
    print_spaces(gap+clip0.k_stop - cd0->bar[cd0->cell[clip0.c]]);
    print_clip( cd1,&clip1 );
#endif
    // copy subbar from cd1
    tp->Lb  = tp->b;             // preserve only branch offsets
    tp->Rb  = tp->b+1;           // within the bar itself
    copy_fix_digits(tp,cd1,&clip1);
#ifdef DEBUG
    print_spaces( gap );
    print_template( tp );
#endif
  }
  // copy from end of subbar to end of code
  clip0.k_stop = cd0->last_codon;
  if( portion == WHOLE || portion == BACK ) {
    clip0.b++;   // skip to next bar
    if( clip0.b >= cd0->cell[clip0.c+1] ) {
      clip0.c++; // skip to next cell
    }
  }
  // don't need Lb,Rb because we will copy verbatim
  copy_verbatim(tp,cd0,&clip0);
#ifdef DEBUG
  //if(!strcmp(cd0->codon,tp->codon)) {
  if( tp->level.g == NON_ALIGNED ) {
    //print_debug = TRUE;
  }
  flush_debug();
#endif
}

/********************************************************//**
   Randomly insert a jump instruction, or a branch sequence
*/
void branch_crossover(Template *tp,Code *cd0,Code *cd1)
{
  Clip clip0;
#ifdef DEBUG
  print_line();
  di += sprintf( &ds[di],"BRANCH CROSSOVER\n");
#endif
  clip_cell(tp,cd0,&clip0,0);

  if( cd0->num_cells > 1 && random()% 4 == 0 ) {
    insert_jump( tp,cd0,&clip0 );
  }
  if( clip0.k == 0 ) { // insert_jump unsuccessful, or not attempted
    insert_branch( tp,cd0,&clip0,cd1 );
  }

  // copy from end of cell to end of code
  clip0.k_stop = cd0->last_codon;
  copy_verbatim(tp,cd0,&clip0);
#ifdef DEBUG
  if(!strcmp(cd0->codon,tp->codon)) {
    // print_debug = TRUE;
  }
  flush_debug();
#endif
}

/********************************************************//**
   Replace the front, back, middle or fringe of a cell in cd0
   with the front, back, middle or fringe of a cell from cd1.
   If align_cell is specified, choose same cell for c0,c1
   and try to align the crossover point(s) in the code;
   otherwise choose the cell randomly (by focus) and choose
   the crossover point(s) randomly.
*/
void cell_crossover(
                    Template *tp,
                    Code *cd0,
                    Code *cd1
                   )
{
  Clip clip0;
  int c0;
#ifdef DEBUG
  print_line();
  di += sprintf( &ds[di]," CELL CROSSOVER\n");
#endif

  c0 = choose_block_by_focus( cd0,1 );

  if( tp->level.g == ALIGNED ) {
    while( identical_blocks(cd0,cd1,c0,1)) {
      c0 = choose_block_by_focus( cd0,1 );
    }
    if(   cd1->cell[c0+1] - cd1->cell[c0]
       != cd0->cell[c0+1] - cd0->cell[c0] ) {
      tp->level.g = NON_ALIGNED;
    }
  }

  clip_cell( tp,cd0,&clip0,0 );

  cross_cells( tp,cd0,&clip0,c0,cd1 );

  if( tp->level.g != ALIGNED && c0 < cd0->num_cells-1
     && power2[cd0->level.m-BLOCK] > cd0->num_cells
                             && random()% 4 == 0 ) {
    if( random_bit()) {
      //c0 = c+1 + random()%( cd0->num_cells - (c+1));
      insert_jump( tp,cd0,&clip0 );
    }
    else {
      insert_branch( tp,cd0,&clip0,cd1 );
    }
  }

  // copy from end of cell to end of code
  clip0.k_stop = cd0->last_codon; // IS THIS NEEDED ????
  copy_verbatim(tp,cd0,&clip0);
#ifdef DEBUG
  if(!strcmp(cd0->codon,tp->codon)) {
    // print_debug = TRUE;
  }
  flush_debug();
#endif
}

/********************************************************//**
   Replace a randomly chosen block of cells from cd0
   with a randomly chosen block of the same size from cd1.
*/
void block_crossover(
                     Template *tp,
                     Code *cd0,
                     Code *cd1,
                     int block_size
                    )
{
  Clip clip0,clip1;
  int c0,c1=0,m,n;
#ifdef DEBUG
  int c;
  print_line();
  di += sprintf( &ds[di],"BLOCK CROSSOVER\n");
#endif
  c0 = choose_block_by_focus( cd0,block_size );

  if( tp->level.g == ALIGNED ) {
    int num_choices = cd0->num_cells / block_size;
    if( !identical_blocks(cd0,cd1,0,num_choices*block_size)) {
      while( identical_blocks(cd0,cd1,c0,block_size)) {
        c0 = choose_block_by_focus( cd0,block_size );
      }
    }
    c1 = c0;
  }
  if( tp->level.g == NON_ALIGNED ) {
    c1 = block_size * similar_block(
                c0/block_size, cd1->num_cells/block_size );
    if( cd1->num_cells == cd0->num_cells && c1 == c0 ) {
      tp->level.g = ALIGNED;
    }
  }
  else if( tp->level.g == TRANSGENIC ) {
    c1 = block_size * similar_block(
                c0/block_size, cd1->num_cells/block_size );
  //c1 = block_size * random()%(cd1->num_cells/block_size);
  }

  // [c0..c0+block_size-1] will be focus block for tp
  tp->Fc = c0;
  tp->Nc = block_size;
  //copy_focus(tp,cd0);

  // adjust focus bar(s), as appropriate
  m=0;
  for( n=0; n < tp->num_focus_bars; n++ ) {
    if( tp->Fb[n].b < cd0->cell[c0] ) {
      tp->Fb[m++] = tp->Fb[n];
#ifdef DEBUG
      printf("bar %d kept\n",tp->Fb[n].b);
#endif
    }
    else if( tp->Fb[n].b >= cd0->cell[c0+block_size] ) {
      tp->Fb[m] = tp->Fb[n];
#ifdef DEBUG
      printf("bar %d",tp->Fb[m].b);
#endif
      tp->Fb[m++].b += (cd1->cell[c1+block_size] - cd1->cell[c1])
                     - (cd0->cell[c0+block_size] - cd0->cell[c0]);
#ifdef DEBUG
        printf(" -> %d\n",tp->Fb[m-1].b);
#endif
    }
    // otherwise, focus bar is removed
#ifdef DEBUG
    else {
      printf("bar %d removed\n",tp->Fb[n].b);
    }
#endif
  }
  tp->num_focus_bars = m;

  // copy cells [0..c0-1] from cd0
  clip_cell(tp,cd0,&clip0,0);
  clip0.k_stop = cd0->bar[cd0->cell[c0]];
  copy_verbatim(tp,cd0,&clip0);

#ifdef DEBUG
  for( c=c1; c < c1 + block_size; c++ ) {
    print_parent_cell( cd1,c );
  }
#endif
  // copy cells [c1..c1+block_size-1] from cd1
  clip_cell(tp,cd1,&clip1,c1);
  clip1.k_stop = cd1->bar[cd1->cell[c1+block_size]];
  copy_fix_digits(tp,cd1,&clip1);

  // prepare to copy cells from cd0, starting from end of block
  clip_cell(tp,cd0,&clip0,c0+block_size);

  if( tp->level.g != ALIGNED && c0+block_size < cd0->num_cells-1
                && power2[cd0->level.m-BLOCK] > cd0->num_cells
                && random_bit()) {
      //     && random() % 4 == 0 ) {
    if( random_bit()) {
      //c0 = c+1 + random()%( cd0->num_cells - (c+1));
      insert_jump( tp,cd0,&clip0 );
    }
    else {
      insert_branch( tp,cd0,&clip0,cd1 );
    }
  }
  /*
  // optionally, insert a jump instruction to jump to head of new block
  // choose jump insertion point
  k0 = cd0->bar[cd0->cell[c0+block_size]];
  k1 = cd0->last_codon;
  // check it is CHAMP(?)
  if( k1 > k0 && tp->level.g != ALIGNED ) {
    k = k0 + random()%( 2*(cd0->bar[cd0->cell[c0]] + k1 - k0 ));
    if( k < k1 ) {
      int op;
      do { // choose the last codon in an instruction
        k = clip0.k + random()%( k1 - clip0.k );
        op = cd0->codon[k];
      } while(    is_dot_or_digit(op)
              ||( op == '#' && cd0->codon[k+1] == '-' ));

      //clip0.k_stop = end_of_instruction(cd0,k);
      clip0.k_stop = k;
      copy_verbatim(tp,cd0,&clip0);

      //add_focus_bar( tp,tp->b );

      insert_integer( tp,c0+block_size-1 );
      insert_operator(tp,'j');
    }
  }
  */
  // copy from end of cell to end of code
  clip0.k_stop = cd0->last_codon;
  copy_verbatim(tp,cd0,&clip0);
#ifdef DEBUG
  print_debug = TRUE;
  flush_debug();
#endif
}

/********************************************************//**
   Replace the front, back, middle or fringe of cell c in cd0
   with the front, back, middle or fringe of a cell from cd1.
   If level.g is ALIGNED, choose same-numbered cell and
   align the crossover points within the cell;
   otherwise, choose cell from cd1 randomly (by focus)
   and choose crossover points (uniformly) randomly.
*/
void cross_cells(
                 Template *tp,
                 Code *cd0,
                 Clip *cp0,
                 int   c,
                 Code *cd1
                )
{
  Clip clip1;
  int portion;
  int focus_bar;
  int b0,b1,k,m,n;
#ifdef DEBUG
  int gap;
#endif

  // c will become the focus cell for tp
  tp->Fc = c;
  tp->Nc = 1;
  //tp->num_focus_bars is preserved from parent

  // copy from cd0, up to start of cell c
  cp0->k_stop = cd0->bar[cd0->cell[c]];
  copy_verbatim( tp,cd0,cp0 );

  if( tp->level.g == ALIGNED ) {
    do {
      if( tp->level.g == TRANSGENIC ) { //|| tp->level.m == JUMP){
        portion =     random()% 5; // allow WHOLE as well
      }
      else {
        portion = 1 + random()% 4; // FRONT, MIDDLE, BACK or FRINGE
      }
      random_clip( cd0,cp0,CELL,portion,0 );
      clip1.c = c; // find matching clip from same same-numbered cell
      matching_clip( cd0,cp0,cd1,&clip1,CELL,portion );
    } while(identical_clips(cd0,cp0,cd1,&clip1,CELL,portion));
  }
  else {
    if( tp->level.g == TRANSGENIC ) { //|| tp->level.m == JUMP){
      portion =     random()% 5; // allow WHOLE as well
    }
    else {
      portion = 1 + random()% 4; // FRONT, MIDDLE, BACK or FRINGE
    }
    random_clip( cd0,cp0,CELL,portion,0 );
    clip1.c = similar_block( c,cd1->num_cells );
    random_clip( cd1,&clip1,CELL,portion,1 );
  }

  if( portion == FRINGE ) { // replace front and back of cell
    // adjust the number of bars in the resulting cell
    b0 = tp->cell[c] + (clip1.Rb - cd1->cell[clip1.c]);
    b1 = b0 + (cp0->Rb-1 - cp0->Lb);
    tp->cell[c+1] = b1 + (cd1->cell[clip1.c+1] - clip1.Lb);

    // adjust focus bar(s), as appropriate
    m=0;
    for( n=0; n < tp->num_focus_bars; n++ ) {
      focus_bar = tp->Fb[n].b - tp->cell[c];
      if( focus_bar < 0 ) { // <= 0 ???
        tp->Fb[m++] = tp->Fb[n];
#ifdef DEBUG
        printf("bar %d kept\n",tp->Fb[n].b);
#endif
      }
      else if( focus_bar >= cd0->cell[c+1] - cd0->cell[c] ) {
        tp->Fb[m] = tp->Fb[n];
#ifdef DEBUG
        printf("bar %d\n",tp->Fb[m].b);
#endif
        tp->Fb[m++].b += (  tp->cell[c+1] -  tp->cell[c] )
                       - ( cd0->cell[c+1] - cd0->cell[c] );
#ifdef DEBUG
        printf(" -> %d\n",tp->Fb[m-1].b);
#endif
      }
      else if(   focus_bar >= cp0->Lb - cd0->cell[c]
              && focus_bar <  cp0->Rb - cd0->cell[c] ) {
        tp->Fb[m] = tp->Fb[n];
#ifdef DEBUG
        printf("bar %d\n",tp->Fb[m].b);
#endif
        tp->Fb[m++].b += ( clip1.Rb - cd1->cell[clip1.c])
                       - (  cp0->Lb - cd0->cell[cp0->c] );
#ifdef DEBUG
        printf(" -> %d\n",tp->Fb[m-1].b);
#endif
      }
      // otherwise, focus bar is removed
#ifdef DEBUG
      else {
        printf("bar %d removed\n",tp->Fb[n].b);
      }
#endif
    }
    tp->num_focus_bars = m;

    add_focus_bar( tp,b0 );
    add_focus_bar( tp,b1 );
#ifdef DEBUG
    print_clip_template( cd0,cp0,tp );
    di += sprintf(&ds[di],"\n");
    //print_spaces( gap-(clip1.k_stop - cd1->bar[clip1.b]));
    print_spaces(3+cd0->bar[cp0->b]-cd0->bar[cd0->cell[cp0->c]]);
    for( k  = cd1->bar[cd1->cell[clip1.c]];
         k <= cd1->bar[cd1->cell[clip1.c+1]]; k++ ) {
      di += sprintf(&ds[di],"%c",cd1->codon[k]);
      if( k == clip1.k ) {
        di += sprintf(&ds[di]," ");
      }
      if( k == clip1.k_stop ) {
        di += sprintf(&ds[di]," ");
      }
    }
    di += sprintf(&ds[di],"\n");
    //print_spaces(gap+cp0->k_stop - cd0->bar[cd0->cell[cp0->c]]);
#endif
    // copy the front of the cell from cd1
    k       = clip1.k; // keep for later use
    clip1.k = cd1->bar[cd1->cell[clip1.c]];
    tp->Rb  = b0;
    tp->Lb  = b1;
    copy_fix_digits(tp,cd1,&clip1);
    clip1.k = k; // prepare to copy end of cell from cd1
    clip1.b = clip1.Lb;

    // copy preserved middle section of cell from cd0
    tp->Lb = b0;
    tp->Rb = b1+1;
    copy_fix_digits(tp,cd0,cp0);

    // copy the back of the cell from cd1
    clip1.k_stop = cd1->bar[cd1->cell[clip1.c+1]]-1;
    if(is_dot_or_digit(cd1->codon[clip1.k_stop])) {
      clip1.k_stop--;           // bar ends with .| or 8|
    }
    tp->Rb = b0;
    tp->Lb = b1;
    copy_fix_digits(tp,cd1,&clip1);
    fix_bar_codons(tp,cd1->codon[clip1.k_stop+1],FALSE);
    cp0->c++; // move to next cell
    cp0->b = cd0->cell[cp0->c];
    cp0->k = cd0->bar[cp0->b];
#ifdef DEBUG
    print_completed_cell( tp,tp->c );
#endif
    tp->c++;
    reset_cell_bars(tp,cd0,cp0); // reset Lb,Rb,cell[c+1]
  }
  else { // not FRINGE
    // adjust the number of bars in the resulting cell
    b0 = cp0->Rb + tp->cell[c] - cd0->cell[c];
    b1 = b0 + ( clip1.Rb-1 - clip1.Lb );
    tp->cell[c+1] = b1 + cd0->cell[c+1] - cp0->Lb;
    //tp->cell[tp->c+1] +=  ( clip1.Rb-1 - clip1.Lb )
    //                    - (  cp0->Lb   -  cp0->Rb );
    m=0;
    for( n=0; n < tp->num_focus_bars; n++ ) {
      focus_bar = tp->Fb[n].b - tp->cell[c];
      if( focus_bar <= cp0->Rb - cd0->cell[c] ) {
        tp->Fb[m++] = tp->Fb[n];
#ifdef DEBUG
        printf("bar %d kept\n",tp->Fb[n].b);
#endif
      }
      else if( focus_bar >= cp0->Lb - cd0->cell[c] ) {
        tp->Fb[m] = tp->Fb[n];
#ifdef DEBUG
        printf("bar %d",tp->Fb[n].b);
#endif
        tp->Fb[m++].b += ( clip1.Rb-1 - clip1.Lb )
                       - (  cp0->Lb   -  cp0->Rb );
#ifdef DEBUG
        printf("-> %d\n",tp->Fb[m-1].b);
#endif
      }
      // otherwise, focus bar is removed
#ifdef DEBUG
      else {
        printf("bar %d removed\n",tp->Fb[n].b);
      }
#endif
    }
    tp->num_focus_bars = m;

    if( !( portion == FRONT || portion == WHOLE )) {
      add_focus_bar( tp,cp0->Rb + tp->cell[c] - cd0->cell[c] );
    }
    if( !( portion == WHOLE || portion == BACK )) {
      add_focus_bar( tp,cp0->Rb + clip1.Rb-1 - clip1.Lb
                                + tp->cell[c] - cd0->cell[c] );
    }
    // copy from start of cell to start of subcell
    k      = cp0->k; // save for later use
    cp0->b = cd0->cell[cp0->c];
    cp0->k = cd0->bar[cp0->b];
    tp->Rb = b0;
    tp->Lb = b1;
    copy_fix_digits(tp,cd0,cp0);
    cp0->k = k; // prepare to copy from end of subcell
#ifdef DEBUG
    gap = print_clip_template( cd0,cp0,tp );
    di += sprintf(&ds[di],"\n");
    print_spaces(gap+cp0->k_stop - cd0->bar[cd0->cell[cp0->c]]);
    print_clip( cd1,&clip1 );
#endif
    // copy subcell from cd1
    tp->Lb = tp->b;
    tp->Rb = tp->Lb + ( clip1.Rb - clip1.Lb );
    copy_fix_digits(tp,cd1,&clip1);
#ifdef DEBUG
    print_spaces( gap );
    print_template( tp );
#endif
    // copy from end of subcell to end of cell
    tp->Rb = b0;
    tp->Lb = b1;
    cp0->b = cp0->Lb;
    cp0->k_stop = cd0->bar[cd0->cell[cp0->c+1]];
    copy_fix_digits(tp,cd0,cp0);
  }
  //if( portion == WHOLE || portion == BACK ) {
  //  cp0->c++;
  //}
  // prepare to start copying from next cell in cd0
  clip_cell(tp,cd0,cp0,c+1);
}

/********************************************************//**
   If CELL or BLOCK crossover, insert a jump instruction
   from a cell >= tp->c, to jump to cell tp->c
  (i.e. mutated cell, or last cell in transplanted block).
   If BRANCH crossover, try to insert a jump instruction
   in a location chosen by focus.
*/
void insert_jump(
                 Template *tp,
                 Code *cd0,
                 Clip *cp0
                )
{
  int Fc[MAX_LEVEL];   // start cell (telescoping focus blocks)
  int Nc[MAX_LEVEL];   // number of cells
  int locus; // POINT, BAR or CELL
  int r_max;
  int op;
  int b0,c0=0,c,k0,k1,m,n;
  /*
  switch( random()% 6 ) {
    case 0: case 1: case 2:
      locus = POINT; break;
    case 3: case 4:
      locus = BAR;   break;
    case 5:
      locus = CELL;  break;
  }
  */
  locus = POINT;

  if( tp->level.m == BRANCH ) {
    if( locus == CELL ) {     // jump replaces subcell
      c0 = choose_block_by_focus( cd0,1 );
    }
    else if( locus == BAR ) { // jump replaces subbar
      b0 = choose_bar_by_focus( cd0 );
    }
    else { // locus == POINT,    jump is simply inserted
      // choose insertion point k0 (=k1)
      r_max = adjust_focus( cd0,Fc,Nc );
      k0 = k1 = choose_codon( cd0,Fc,Nc,r_max,TRUE );
      while( cd0->bar[cd0->cell[c0+1]] <= k0 ) {
        c0++;   // find cell c0 containing codon k0
      }
      b0 = cd0->cell[c0];
      while( cd0->bar[b0+1] <= k0 ) {
        b0++;   // find bar  b containing codon k0
      }
    }
    if( c0 == 0 ) { // no opportunity to jump down
      return;
    }
  }
  else { // CELL or BLOCK
    if( locus == CELL ) { // choose c0 > c uniformly randomly
      c0 = tp->c + random()%( cd0->num_cells - tp->c );
    }
    else { // POINT or BAR,  choose b0 > cp0->b randomly
      int b = cd0->cell[tp->c];
      b0 = b + random()%( cd0->cell[cd0->num_cells] - b );
      c0 = cell_for_bar( cd0,b0 );
    }
  }

  // choose which cell to jump to
  if( tp->level.m == BRANCH ) {
    //c = choose_jump_index( c0,cd0->num_cells );
    c = random()% c0; // choose c < c0 randomly
  }
  else {
    c = tp->c-1;      // jump to most recently copied cell
  }
  // copy from cd0, up to start of cell c0
  cp0->k_stop = cd0->bar[cd0->cell[c0]];
  copy_verbatim(tp,cd0,cp0);

  if( locus == CELL ) {
    random_clip(cd0,cp0,CELL,MIDDLE,0);
    k0 = cp0->k_stop; // start of subcell
    k1 = cp0->k;      //  end  of subcell
    tp->Rb = tp->b + cp0->Rb - cp0->b; //cd0->cell[c];
    tp->Lb = tp->Rb;
    tp->cell[c0+1] -= cp0->Lb - cp0->Rb;
    cp0->k = cd0->bar[cd0->cell[c0]];

    m=0;    // adjust focus bars
    for( n=0; n < tp->num_focus_bars; n++ ) {
      if( tp->Fb[n].b <= tp->Rb ) {
        tp->Fb[m++] = tp->Fb[n];
#ifdef DEBUG
        printf("bar %d kept\n",tp->Fb[n].b);
#endif
      }
      else if( tp->Fb[n].b >= tp->Rb + cp0->Lb - cp0->Rb ) {
        tp->Fb[m] = tp->Fb[n];
#ifdef DEBUG
        printf("bar %d",tp->Fb[n].b);
#endif
        tp->Fb[m++].b -= cp0->Lb - cp0->Rb;
#ifdef DEBUG
        printf(" -> %d\n",tp->Fb[m-1].b);
#endif
      }
      // otherwise, focus bar is removed
#ifdef DEBUG
      else {
        printf("bar %d removed\n",tp->Fb[n].b);
      }
#endif
    }
    tp->num_focus_bars = m;
  }
  else { // locus is POINT or BAR
    // copy from start of cell to start of bar
    cp0->k_stop = cd0->bar[b0];
    copy_verbatim(tp,cd0,cp0);
    if( locus == BAR ) {
      // choose random subbar from cd
      random_clip(cd0,cp0,BAR,MIDDLE,0);
      //exclude_final_bar(cd0,cp0,BAR,MIDDLE,0); ???
      k0 = cp0->k_stop; // start of subbar
      k1 = cp0->k;      //  end  of subbar
    }
    else { // ( locus == POINT )
      // choose k0=k1 randomly within b0
      k0 = cd0->bar[b0];
      do { // choose the last codon in an instruction
        k1 = k0 + random()%( cd0->bar[b0+1] - k0 );
        op = cd0->codon[k1];
      } while(    is_dot_or_digit(op)
              ||( op == '#' && cd0->codon[k1+1] == '-' ));
      k0 = k1;
    }
  }
  // copy up to insertion point
  cp0->k_stop = k0;
  copy_fix_digits( tp,cd0,cp0 );

  add_focus_bar( tp,tp->b );

  insert_integer( tp,c );
  insert_operator(tp,'j');

  // copy to end of cell
  if( locus == CELL ) {
    cp0->b = cp0->Lb;
  }
  cp0->k = k1;
  cp0->k_stop = cd0->bar[cd0->cell[c0+1]];
  copy_fix_digits(tp,cd0,cp0);

  /*
  k0 = cd0->bar[cd0->cell[c0+block_size]];
  k1 = cd0->last_codon;
  // check it is CHAMP(?)
  if( k1 > k0 && tp->level.g != ALIGNED ) {
    k = k0 + random()%( 2*(cd0->bar[cd0->cell[c0]] + k1 - k0 ));
    if( k < k1 ) {
      int op;
      do { // choose the last codon in an instruction
        k = clip0.k + random()%( k1 - clip0.k );
        op = cd0->codon[k];
      } while(    is_dot_or_digit(op)
              ||( op == '#' && cd0->codon[k+1] == '-' ));

      //clip0.k_stop = end_of_instruction(cd0,k);
      clip0.k_stop = k;
      copy_verbatim(tp,cd0,cp0);

      //add_focus_bar( tp,tp->b );

      insert_integer( tp,c0+block_size );
      insert_operator(tp,'j');
    }
  }
  */
}

/********************************************************//**
   Insert a branch forward sequence at cell c0 and
  (optionally) include a jump to cell c.
   If level is BRANCH, choose c0 by focus and c < c0 uniformly.
   If level is CELL or BLOCK, c = mutated cell or top cell
   in transplanted block, and c0 > c is chosen uniformly.

   insert:

            +--------- ------+
            /       /   \    \
          +---------g:..|------+

   skip:

            +--------- ----- ----+
            /       /  |   | |  |
          +---------g: -----|----+

   insert + skip:

          +---------   ----- ----+
                   /    \   \ \  \
        +---------g:..1:|-----|----+
*/
void insert_branch(
                   Template *tp,
                   Code *cd0,
                   Clip *cp0,
                   Code *cd1
                  )
{
  int Fc[MAX_LEVEL];
  int Nc[MAX_LEVEL];
  Clip clip1;
  Boolean insert_code,skip_code;
  Boolean include_jump = ( tp->level.m > BRANCH );
  int portion_insert=0;
  int portion_skip;
  int new_final_bar,num_new_bars;
  int r_max;
  int op;
  int b,b0,c=0,c0=0,k0,k1,n;
#ifdef DEBUG
  int gap;
#endif

  // must either insert something or skip code (or both)
  skip_code   =  random_bit(); // random()% 5 < 3 ???
  insert_code = !include_jump &&( !skip_code || random_bit());

  // plan to insert code in the first instance, may later
  // decide to insert jump instruction instead of code
  if( insert_code ) {
    // choose subbar uniformly randomly from cd1
    clip1.b = random()%(cd1->cell[cd1->num_cells]);
    clip1.c = cell_for_bar( cd1,clip1.b );
    //portion_insert = random_bit() ? MIDDLE : BACK ;
    portion_insert = random()% 4;
    random_clip(cd1,&clip1,BAR,portion_insert,1);
    //exclude_final_bar(cd1,&clip1,BAR,portion,1);
    if( clip1.k == clip1.k_stop ) { // no code to insert
      insert_code = FALSE;
      skip_code   = TRUE;
    }
  }

  b = cd0->cell[tp->c]; // first bar to be copied from c0

  // choose cell c0, bar b0 and subbar [k0+1..k1]
  if( skip_code ) {
    if( tp->level.m == BRANCH ) {
      b0 = choose_bar_by_focus( cd0 );
    }
    else { // level is CELL or BLOCK
      b0 = b + random()%( cd0->cell[cd0->num_cells] - b );
    }
    c0 = cell_for_bar( cd0,b0 );
    cp0->b = b0;
    cp0->c = c0;
    //portion_skip = 1 + random()% 3; // FRONT, MIDDLE or BACK
    portion_skip = random()% 4;
    random_clip(cd0,cp0,BAR,portion_skip,0);
    exclude_final_bar(cd0,cp0,BAR,portion_skip,0);
    k0 = cp0->k_stop; // start of subbar
    k1 = cp0->k;      //  end  of subbar
    if( k0 == k1 ) {  //  no code to skip
      skip_code = FALSE;
    }
  }
  else if( tp->level.m == BRANCH ) {
    // choose insertion point k0 (=k1) by focus
    r_max = adjust_focus( cd0,Fc,Nc );
    k0 = k1 = choose_codon( cd0,Fc,Nc,r_max,TRUE );
    while( cd0->bar[cd0->cell[c0+1]] <= k0 ) {
      c0++;   // find cell c0 containing codon k0
    }
    b0 = cd0->cell[c0];
    while( cd0->bar[b0+1] <= k0 ) {
      b0++;   // find bar  b0 containing codon k0
    }
  }
  else {
    // choose bar uniformly randomly
    b0 = b + random()%(cd0->cell[cd0->num_cells] - b);
    c0 = cell_for_bar( cd0,b0 );
    // choose k0=k1 randomly within b0
    k0 = cd0->bar[b0];
    do { // choose the last codon in an instruction
      k1 = k0 + random()%( cd0->bar[b0+1] - k0 );
      op = cd0->codon[k1];
    } while(    is_dot_or_digit(op)
            ||( op == '#' && cd0->codon[k1+1] == '-' ));
    k0 = k1;
  }

  // choose whether to include_jump
  if( c0 > 0 && !include_jump &&( random()% 3 == 0 )) {
    insert_code = FALSE; // can't do both
    include_jump = TRUE;
  }

  new_final_bar = (end_of_instruction(cd0,k1+1) < cd0->bar[b0+1]);
  num_new_bars = ((insert_code || include_jump)&& skip_code)+ new_final_bar;

  /*   exclude final bar of the inserted code if:
   (a) there is already a bar following it, or
   (b) it is not preceded by . or 8 (in this case the bar will be
       inserted later, and . or 8 may be included, randomly)
  */
  if(       insert_code
     &&( portion_insert == WHOLE || portion_insert == BACK )
     &&(   !new_final_bar
        || !is_dot_or_digit(cd1->codon[clip1.k_stop-1]))) {
    exclude_final_bar(cd1,&clip1,BAR,portion_insert,1);
    //clip1.k_stop--;
  }

  if( include_jump ) { // choose which cell to jump to
    if( tp->level.m == BRANCH ) {
      do {
        c = choose_jump_index( c0,cd0->num_cells );
      } while( c == c0 );
    }
    else {         // jump to mutated cell,
      c = tp->c-1; // or highest cell in transplanted block
    }
  }

  b = tp->b + b0 - b; // modified bar in new code

  // adjust focus bars
  for( n=0; n < tp->num_focus_bars; n++ ) {
    if( tp->Fb[n].b > b ) {
#ifdef DEBUG
      printf("bar %d -> %d\n",tp->Fb[n].b,tp->Fb[n].b+num_new_bars);
#endif
      tp->Fb[n].b += num_new_bars;
    }
#ifdef DEBUG
    else {
      printf("bar %d kept\n",tp->Fb[n].b);
    }
#endif
  }
  add_focus_bar( tp,b ); // + new_initial_bar
  if( num_new_bars > 0 ) {
    add_focus_bar( tp,b+1 ); // + new_initial_bar
  }

  clip_cell(tp,cd0,cp0,tp->c);
  /*
  if( include_jump && tp->level.m == BRANCH ) {
    do {
      c = choose_jump_index( c0,cd0->num_cells );
    } while( c == c0 );
    if( c < c0 ) {
      cross_cells( tp,cd0,&clip0,c,cd1 );
    }
  }
  //else {
  //  tp->Fc = c;   NOT NEEDED ????
  //  tp->Nc = 1;
  //}
  */
  // copy from cd0, up to start of cell c0
  cp0->k_stop = cd0->bar[cd0->cell[c0]];
  copy_verbatim(tp,cd0,cp0);

  // copy from start of cell to start of subbar
  cp0->k_stop = k0;          //      k0 k1
  cp0->Rb = b0;              // ----+.....+-----
  cp0->Lb = b0;              //      R = L
  tp->Rb  = b; // + tp->b - clip0.b;
  tp->Lb  = tp->Rb + num_new_bars;
  tp->cell[c0+1] += num_new_bars;
  copy_fix_digits( tp,cd0,cp0 );

  //if( new_initial_bar ) {
  //  fix_bar_codons( tp,BLN,FALSE );
  //}
#ifdef DEBUG
  cp0->k = k1;
  gap = print_clip_template( cd0,cp0,tp );
  cp0->k = k0;
  di += sprintf(&ds[di],"\n");
#endif

  // insert comparison and branch forward
  insert_comparison( tp,BRF );
  insert_operator( tp,BRF );

  if( insert_code || include_jump ) {
#ifdef DEBUG
    print_spaces(gap+1+cp0->k_stop - cd0->bar[cd0->cell[cp0->c]]);
    if( insert_code ) print_clip( cd1,&clip1 );
#endif
    if( include_jump ) {
      insert_integer( tp,c );
      insert_operator(tp,'j');
    }
    else { // insert_code
      // copy selected subbar from cd1
      tp->Lb = tp->b;
      tp->Rb = tp->b+1;
      copy_fix_digits( tp,cd1,&clip1 );

      // prepare for further copying from cd0
      tp->Rb = b;
      tp->Lb = b + num_new_bars;
    }

    if( skip_code && tp->codon[tp->k] != '|' ) { // if .. else
      //switch( random()% 3 ) {
      //case 0: insert_dot( tp );       break;
      //case 1: insert_integer( tp,8 ); break;
      //case 2:
      insert_integer( tp,1 );
      insert_operator(tp,BRF);
        //  break;
        //}
      fix_bar_codons( tp,BLN,FALSE );
    }
  }

  //if( skip_code ) {
  // copy subbar from cd0 (Lb,Rb remain as before)
  cp0->k_stop = k1; // if( k0==k1 ), nothing will happen
  copy_fix_digits( tp,cd0,cp0 );
  //}

  //if( new_final_bar ) {
  if( new_final_bar && tp->codon[tp->k] != '|' ) {
    if( random_bit()) {
      if( random_bit())
        insert_dot( tp );
      else
        insert_integer(tp,8);
    }
    fix_bar_codons( tp,BLN,FALSE );
  }
  // copy to end of cell
  cp0->k_stop = cd0->bar[cd0->cell[c0+1]];
  copy_fix_digits( tp,cd0,cp0 );
#ifdef DEBUG
  print_spaces( gap );
  print_template( tp );
#endif
  //if( include_jump && c0 > c ) {
  //  cross_cells( tp,cd0,cp0,c0,cd1 );
  //}
  // copy to end of code
  cp0->k_stop = cd0->last_codon;
  copy_verbatim( tp,cd0,cp0 );
  if( tp->overflow ) {
    reset_template( tp ); // just copy with no changes
    copy_focus( tp,cd0 );
    copy_code( tp,cd0 );
  }
#ifdef DEBUG
  if(strcmp(cd0->codon,tp->codon)) {
    //print_debug = TRUE;
  }
  flush_debug();
#endif
}

/********************************************************//**
   Randomly choose an internal bar line and remove it.
*/
void bar_line_removal( Template *tp, Code *cd )
{
  Clip clip;
  int Fc[MAX_LEVEL]; // telescoping focus blocks
  int Nc[MAX_LEVEL];
  int r_max;
  int b,c,n,r;
#ifdef DEBUG
  int gap;
  print_line();
  di += sprintf( &ds[di],"BAR LINE REMOVAL\n");
#endif
  clip_cell(tp,cd,&clip,0); // prepare to start copying

  if(cd->cell[cd->num_cells] > cd->num_cells){// more bars than cells

    r_max = adjust_focus( cd,Fc,Nc ); // telescoping focus blocks
    // choose a block with at least one internal bar
    //r = random_geometric(2.0,r_max);
    r = 0;
    while( r < r_max && random_bit()) {
      r++;
    }
    while(cd->cell[Fc[r]+Nc[r]] == cd->cell[Fc[r]]+Nc[r]){
      r++;
    }
    // choose internal bar uniformly randomly, within the block
    b =            cd->cell[Fc[r]]+1
       + random()%(cd->cell[Fc[r]+Nc[r]] - (cd->cell[Fc[r]]+Nc[r]));
    c = Fc[r];
    while( cd->cell[c+1] <= b ) {
      c++;
      b++;
    }
    // copy from beginning to start of cell
    clip.k_stop = cd->bar[cd->cell[c]];
    copy_verbatim(tp,cd,&clip);
    clip.Rb = b-1;     //    -----|-----
    clip.Lb = b;       //      Rb   Lb
    tp->Rb  = clip.Rb; //    -----------
    tp->Lb  = tp->Rb;  //      Rb = Lb
    tp->cell[c+1]--;
    // copy from start of cell to (last instruction before) bar
    clip.k_stop = cd->bar[b]-1;
    while(is_dot_or_digit(cd->codon[clip.k_stop])) {
      clip.k_stop--;
    }
    copy_fix_digits(tp,cd,&clip);
    // adjust focus bar(s)
    for( n=0; n < tp->num_focus_bars; n++ ) {
      if( tp->Fb[n].b >= b ) {
#ifdef DEBUG
        printf("bar %d -> %d\n",tp->Fb[n].b,tp->Fb[n].b-1);
#endif
        tp->Fb[n].b--;
      }
#ifdef DEBUG
      else {
        printf("bar %d kept\n",tp->Fb[n].b);
      }
#endif
    }
    clip.b++;          // skip to next bar
    clip.k = cd->bar[b];
#ifdef DEBUG
    gap = print_clip_template( cd,&clip,tp );
    di += sprintf(&ds[di],"%c",'\n');
#endif
    clip.k_stop = cd->bar[cd->cell[c+1]];
    copy_fix_digits(tp,cd,&clip); // copy to end of cell
  }
  clip.k_stop = cd->last_codon;
  copy_verbatim(tp,cd,&clip);     // copy to end of code
#ifdef DEBUG
  if(!strcmp(cd->codon,tp->codon)) {
    //print_debug = TRUE;
  }
  flush_debug();
#endif
}

/********************************************************//**
   Check the template for syntax error (debugging).
*/
int check_template( Template *tp )
{
  Clip clip;
  int is_fraction;
  int value;
  int op;
  int d,n,v=0,p=0;
  clip.c = 0; clip.b = 0; clip.k = 1;
  if( tp->overflow ) {
    return( 0 );
  }
  for( n=0; n < tp->num_focus_bars; n++ ) {
    if( tp->Fb[n].b >= tp->b ) {
      printf("B\n");
      return( tp->k );
    }
  }
  while( clip.k < tp->k ) {
    op = tp->codon[clip.k];
    if( op_uses_digits[op] ) {
      d = clip.k;
      if( op == PSH && tp->codon[clip.k+1] == NEG ) {
        clip.k++;
      }
      is_fraction = FALSE;
      while( is_dot_or_digit(tp->codon[d-1])) {
        d--;
        is_fraction |= ( tp->codon[d] == DOT );
      }
      if( tp->ival[d] != clip.k - d ) {
        printf("V\n");
        return( tp->k );
      }
      if(  ( op == PSH || op == NEG )
         &&( is_dot_or_digit(tp->codon[d]))) {
        if( clip.k > d ) {
          if( v >= tp->v || tp->ival[d+1] != v ) {
            printf("C\n");
            return( clip.k );
          }
          v++;
        }
        if( is_fraction ) {
          if( p >= tp->p || d != tp->param[p] ) {
            printf("D\n");
            return( clip.k );
          }
          p++;
        }
      }
      value = 0;
      while(codon_is_digit(tp->codon[d])) {
        value = globals.base*value+digit_from_codon(tp->codon[d++]);
      }
      if( op == BRB && value > clip.b - tp->cell[clip.c] ) {
        printf("E\n");
        return( clip.k );
      }
      if(op == BRF && clip.b + value >= tp->cell[clip.c+1]){
        printf("F\n");
        return( clip.k );
      }
      if( op == BLN ) {
        clip.b++;
        if( tp->bar[clip.b] != clip.k ) {
          printf("G\n");
          return( clip.k );
        }
        if( tp->cell[clip.c+1] == clip.b ) {
          clip.c++;
        }
      }
    }
    else if( !is_dot_or_digit(op) && tp->ival[clip.k] != 0 ) {
      printf("Z\n");
      return( clip.k );
    }
    clip.k++;
  }
  if( clip.k != tp->bar[tp->cell[tp->num_cells]] ) {
    printf("H\n");
    return( clip.k );
  }

  return( 0 );
}

/********************************************************//**
   Initialize clip, to copy from beginning of the specified
   cell in cd, to the beginning of a cell in tp,
   preserving all offsets and jump indices.
   Update cp->c,n,k,Lb,Rb,cell[c+1] (but not cp->p).
*/
void clip_cell(
               Template *tp,
               Code *cd,
               Clip *cp,
               int c
              )
{
  cp->c = c;                   // specified cell
  cp->b = cd->cell[cp->c];     // first bar in cell
  cp->k = cd->bar[cp->b];      // first codon in bar
  cp->k_stop = cd->last_codon; // end of code
  reset_cell_bars(tp,cd,cp);   // set Lb,Rb,cell[c+1]
}

/********************************************************//**
   Copy code from cd->codon[(cp->k+1)..(cp->k_stop)]
               to tp->codon[(tp->k+1)..(..........)].
   Update k,n,c,p,bar[],cell[],param[], but not Fb[]
*/
void copy_verbatim(
                   Template *tp,
                   Code *cd,
                   Clip *cp
                  )
{
  int k_offset;
  int p;
#ifdef DEBUG
  int c=tp->c;
#endif
  while(   cp->k < cp->k_stop
        && cd->codon[cp->k+1] == NEG ) {
    insert_operator( tp,NEG );
    cp->k++;
  }
  if( cp->k < cp->k_stop ) {
    k_offset = tp->k - cp->k;
    if( cd->bar[cd->cell[cp->c+1]] <= cp->k_stop ) {
      // copy cell boundaries, with appropriate offset
      while(   cp->c < cd->num_cells
            && cd->bar[cd->cell[cp->c+1]] <= cp->k_stop ) {
        tp->cell[++tp->c] = cd->cell[++cp->c]+(tp->b - cp->b);
      }
      reset_cell_bars( tp,cd,cp ); // set Lb,Rb,cell[c+1]
    }
    // copy bar boundaries, with appropriate offset
    while(   cp->b < cd->cell[cd->num_cells]
          && cd->bar[cp->b+1] <= cp->k_stop ) {
      next_bar( tp );
      tp->bar[tp->b] = cd->bar[++cp->b] + k_offset;
    }
    // copy param indices, with appropriate offset
    p = get_param_index( cd,cp->k );
    while(             p  <  cd->num_param
          && cd->param[p] <= cp->k_stop ) {
      tp->param[tp->p] = cd->param[p++] + k_offset;
      next_param( tp );
    }
    // copy codons (if space permits)
    if( cp->k_stop + k_offset < MAX_CODE ) {
      memcpy( &tp->codon[tp->k+1],
              &cd->codon[cp->k+1],(cp->k_stop - cp->k));
      memcpy( &tp->ival[tp->k+1],
              &cd->ival[cp->k+1],(cp->k_stop - cp->k)*sizeof(Index));
      while( cp->k < cp->k_stop ) {
        cp->d = cp->k+1;                 // first codon in instruction
        cp->k = cp->d + cd->ival[cp->d]; //  last codon in instruction
        if(  ( cd->codon[cp->k] == PSH || cd->codon[cp->k] == NEG )
           &&  is_dot_or_digit(cd->codon[cp->d])) {
          tp->ival[cp->d+1 + k_offset] = tp->v;
#ifdef ASSERT
          assert( cd->ival[cp->d+1] < cd->num_value );
#endif
          tp->fval[tp->v] = cd->fval[cd->ival[cp->d+1]];
          next_value( tp );
        }
      }
      tp->k = cp->k_stop + k_offset;
    }
    else {
      tp->overflow = TRUE;
    }
#ifdef DEBUG
    while( c < tp->c ) {
      print_completed_cell( tp,c++ );
    }
#endif
    //cp->k = cp->k_stop;
  }
}

/********************************************************//**
   Copy code from cd->codon[(cp->k+1)..(cp->k_stop)]
               to tp->codon[(tp->k+1)..(..........)].
   Check digits along the way, fixing them if necessary.
   Update k,n,c,p,bar[],cell[],param[].
   Within the first (copied) cell, Lb and Rb specify the left
   and right edge of the preserved section(s) of that cell.
   [see fix_branch_offset() for details].
   If copying reaches the end of one cell and starts on a
   new cell (c), Lb and Rb are reset to cell[c],cell[c+1]
  [this is done by fix_digits()].
*/
void copy_fix_digits(
                     Template *tp,
                     Code *cd,
                     Clip *cp
                    )
{
  while( cp->k < cp->k_stop ) {
    cp->d  = cp->k+1;
    cp->k += cd->ival[cp->d]+1;
    if( cp->k > cp->d || cd->codon[cp->k] == BLN ) {
      fix_digits(tp,cd,cp,FALSE);
    }
    else {
      insert_operator( tp,cd->codon[cp->k] );
    }
  }
}

/********************************************************//**
   If b is not already a focus bar for tp, add it.
*/
void add_focus_bar( Template *tp, int b )
{
  int n=0;
  while( n < tp->num_focus_bars && tp->Fb[n].b != b ) {
    n++;
  }
  if( n == tp->num_focus_bars ) {
    tp->Fb[n].m = tp->level.m;
    tp->Fb[n].b = b;
    tp->num_focus_bars++;
#ifdef DEBUG
    printf("bar %d added\n",tp->Fb[n].b);
#endif
  }
}

/********************************************************//**
   Set Fc[r],Nc[r] to define a series of
   telescoping focus blocks expanding outwards from cd->Fc,Nc.
*/
int adjust_focus( Code *cd, int Fc[], int Nc[] )
{
  int r=0;
  Fc[0] = cd->Fc;
  Nc[0] = cd->Nc;
  while( cd->Nc * power2[r] < cd->num_cells ) {
    r++;
    Nc[r] = cd->Nc * power2[r];
    Fc[r] = Nc[r]*( cd->Fc / Nc[r] );
    if( Fc[r] + Nc[r] > cd->num_cells ) {
      Nc[r] = cd->num_cells - Fc[r];
    }
  }
  return r;
}

/********************************************************//**
   Choose codon randomly, giving preference to focus bar(s)
   and telescoping focus block(s).
   If branch_or_jump, choose the last codon in an instruction
  (new code will be inserted after this instruction).
   Otherwise, choose either the last codon or the last
   dot/digit in an instruction (any new code will be inserted
   after the previous instruction).
*/
int choose_codon(
                 Code *cd,
                 int   Fc[],
                 int   Nc[],
                 int   r_max,
                 int   branch_or_jump
                )
{
  int op;
  int k0,k1;
  int k,n,r;
  if( cd->num_focus_bars > 0 && random_bit()) {
    n = random()% cd->num_focus_bars;
    k0 = cd->bar[cd->Fb[n].b]; // choose from focus bar
    k1 = cd->bar[cd->Fb[n].b+1];
  }
  else {
    r = 0;
    while( r < r_max && random_bit()) {
      r++;
    }
    k0 = cd->bar[cd->cell[Fc[r]]];
    k1 = cd->bar[cd->cell[Fc[r]+Nc[r]]];
  }
  if( branch_or_jump ) {
    do { // choose the last codon in an instruction
      k = k0 + random()%( k1 - k0 );   // could be initial bar
      op = cd->codon[k];
    } while(    is_dot_or_digit(op)
            ||( op == '#' && cd->codon[k+1] == '-' ));
  }
  else {
    do { // choose last codon or last dot/digit in an instruction
      k = k0+1 + random()%( k1 - k0 ); // could be final bar
      op = cd->codon[k];
    } while((is_dot_or_digit(op) && is_dot_or_digit(cd->codon[k+1]))
            ||( op == '-' && cd->codon[k-1] == '#' )
            ||( op == '|' && random_bit())); // make bar less likely,
                     //  because cannot REMOVE or EXCHANGE a bar line
  }
  return( k );
}

/********************************************************//**
   Choose a block of the specified size randomly,
   weighted by proximity to focus block.
*/
int choose_block_by_focus( Code *cd, int block_size )
{
  int Fc[MAX_LEVEL];
  int Nc[MAX_LEVEL];
  int r_max;
  int block,r;
  r_max = adjust_focus( cd,Fc,Nc );
  r = choose_index( r_max );
  //block = random()%((cd->Nc * power2[r])/ block_size );
  block = random()%( Nc[r] / block_size );
  return( Fc[r] + block * block_size );
}

/********************************************************//**
   Choose bar randomly, weighted by proximity to focus block.
*/
int choose_bar_by_focus( Code *cd )
{
  int Fc[MAX_LEVEL];
  int Nc[MAX_LEVEL];
  int r_max;
  int b,r;
  r_max = adjust_focus( cd,Fc,Nc );
  r = choose_index( r_max );
  b =            cd->cell[Fc[r]]
     + random()%(cd->cell[Fc[r]+Nc[r]] - cd->cell[Fc[r]]);
  return( b );
}

/********************************************************//**
   Choose number of mutations as well as locus and type of
   each mutation; sort them and ensure there are no conflicts.
*/
int choose_mutations(
                     Code *cd,
                     int   level_m,
                     int   locus[20],
                     int   type[20]
                    )
{
  int Fc[MAX_LEVEL];
  int Nc[MAX_LEVEL];
  int is_trim = ( level_m == TRIM );
  int num_mutations,mutation=0;
  int type_m;
  int r_max;
  int k,m;

  r_max = adjust_focus( cd,Fc,Nc );

  num_mutations = random_geometric(1.0,MAX_MUTATIONS);
  /*
  num_mutations = random_cauchy(1);
  if( num_mutations > MAX_MUTATIONS ) {
    num_mutations = MAX_MUTATIONS;
  }
  */
  // choose locus and type of each mutation
  while( mutation < num_mutations ) {
    k = choose_codon( cd,Fc,Nc,r_max,FALSE );
    if( is_dot_or_digit( cd->codon[k] )) {
      type_m = DIGIT;
      // scan forward to find operator
      while( is_dot_or_digit( cd->codon[k] )) {
        k++;
      }
    }
    else if( cd->codon[k] == BLN ){
      // cannot REMOVE or EXCHANGE a bar line
      if(      is_trim
         ||(  !is_dot_or_digit(cd->codon[k-1])
            && random_bit())) {
        type_m = DIGIT;
      }
      else {
        type_m = INSERT;
      }
    }
    else { // codon is an operator
      int option[4] = { DIGIT,REMOVE,EXCHANGE,INSERT };
      int m0,m1; // mutation type
      // DIGIT mutation allowed only if this operator
      // could take dot/digits, but currently does not
      m0 = (    op_uses_digits[(int)cd->codon[k]]
            && !is_dot_or_digit(cd->codon[k-1])) ? 0 : 1 ;
      m1 = is_trim ? 3 : 4 ; // if TRIM, cannot INSERT
      type_m = option[m0 + random()%(m1 - m0)];
    }
    // keep mutations sorted by locus (location)
    m = mutation;
    while( m > 0 && locus[m-1] > k ) {
      m--;
    }
    // at each location, we allow any number of INSERT mutations
    // followed by at most one DIGIT, REMOVE or EXCHANGE mutation
    if( m > 0 && locus[m-1] == k && type[m-1] != INSERT ) {
      if( type_m == INSERT ) {
        type_m  = type[m-1]; // swap order of mutations
        type[m-1] = INSERT;
      }
      else {
        type_m = NONE; // scrap this mutation
      }
    }
    if( type_m == NONE ) {
      num_mutations--;
    }
    else {
      m = mutation;
      while( m > 0 && locus[m-1] > k ) {
        locus[m] = locus[m-1];
        type[m]  = type[m-1];
        m--;
      }
      locus[m] = k;
      type[m] = type_m;
      mutation++;
    }
  }
  return( num_mutations );
}

/********************************************************//**
   Return TRUE if code contains no instructions; FALSE otherwise.
*/
Boolean code_is_empty( Code *cd )
{
  return( cd->last_codon == cd->num_cells );
}

/********************************************************//**
*/
Boolean identical_blocks(
                         Code *cd0,
                         Code *cd1,
                         int c,
                         int n
                        )
{
  int k0 = cd0->bar[cd0->cell[c]];
  int k1 = cd1->bar[cd1->cell[c]];
  int length = cd0->bar[cd0->cell[c+n]] - k0;
  return(   length == cd1->bar[cd1->cell[c+n]] - k1
         && strncmp(&cd0->codon[k0+1],&cd1->codon[k1+1],length) == 0);
}

/********************************************************//**
   Assume cd0,cd1 have the same number of cells.
   Return TRUE if cell c of cd0,cd1 have the same number of bars
   and the code for the corresponding bars is identical;
   FALSE otherwise
*/
Boolean identical_bars(
                       Code *cd0,
                       Code *cd1,
                       int b0,
                       int c
                      )
{
  if(   cd1->cell[c+1] - cd1->cell[c]
     != cd0->cell[c+1] - cd0->cell[c] ) {
    return FALSE;
  }
  else {
    int b1 = b0 + cd1->cell[c] - cd0->cell[c];
    int k0 = cd0->bar[b0];
    int k1 = cd1->bar[b1];
    int length = cd0->bar[b0+1] - k0;
    return(    length == cd1->bar[b1+1] - k1
           && strncmp(&cd0->codon[k0+1],&cd1->codon[k1+1],length) == 0);
    /*
    int length = last_codon_in_bar(cd0,b0) - k0;
    //printf("b0=%d,b1=%d,k0=%d,k1=%d,length=%d\n",b0,b1,k0,k1,length);
    return(   length == last_codon_in_bar(cd1,b1) - k1
           && strncmp(&cd0->codon[k0+1],&cd1->codon[k1+1],length) == 0);
    */
  }
}

/********************************************************//**
   Assume cp0 is parent 0, cp1 is parent 1.
*/
Boolean identical_clips(
                        Code *cd0,
                        Clip *cp0,
                        Code *cd1,
                        Clip *cp1,
                        int bar_or_cell,
                        int portion
                       )
{
  Boolean eq;
  int k0,k1;
  if( portion == FRINGE ) {
    if( bar_or_cell == BAR ) {
      k0 = cd0->bar[cp0->b];                // start of bar
      k1 = cd1->bar[cp1->b];
    }
    else {
      k0 = cd0->bar[cd0->cell[cp0->c]];     // start of cell
      k1 = cd1->bar[cd1->cell[cp1->c]];
    }
    eq =(   cp0->k - k0 == cp1->k_stop - k1
         && strncmp(&cd0->codon[k0+1],
                    &cd1->codon[k1+1],cp0->k-k0)==0);
    if( eq ) {
      if( bar_or_cell == BAR ) {
        k0 = cd0->bar[cp0->b+1];            // end of bar
        k1 = cd1->bar[cp1->b+1];
      }
      else {
        k0 = cd0->bar[cd0->cell[cp0->c+1]]; // end of cell
        k1 = cd1->bar[cd1->cell[cp1->c+1]];
      }
      eq &=(   k0 - cp0->k_stop == k1 - cp1->k
            && strncmp(&cd0->codon[cp0->k+1],
                       &cd1->codon[cp1->k+1],k1-cp1->k)==0);
    }
  }
  else {
    eq =(   cp0->k - cp0->k_stop == cp1->k_stop - cp1->k
         && strncmp(&cd0->codon[cp0->k+1],
                    &cd1->codon[cp1->k+1],cp0->k-cp0->k_stop)==0);
  }
  return( eq );
}

/********************************************************//**
   Choose either a subbar  from cp->b (BAR)
              or a subcell from cp->c (CELL).
   The clip will be delineated by cp->k and cp->k_stop.
*/
void random_clip(
                 Code *cd,
                 Clip *cp,
                 int bar_or_cell,
                 int portion,
                 int parent  // 0=PRIMARY, 1=SECONDARY
                )
{
  int k_min,k_max,ki,kf;
  if( bar_or_cell == BAR ) {
    k_min = cd->bar[cp->b];            // start of bar
    k_max = last_codon_in_bar(cd,cp->b);
  }
  else {
    k_min = cd->bar[cd->cell[cp->c]];  // start of cell
    k_max = last_codon_in_cell(cd,cp->c);
  }
  if( portion == FRONT || portion == WHOLE ) {
    ki = k_min;                        // start of bar or cell
  }
  else { // portion is BACK, MIDDLE or FRINGE
    ki = k_min + random()%( k_max+1 - k_min );
    ki = end_of_instruction(cd,ki);
  }
  if( portion == WHOLE || portion == BACK ) {
    //kf = k_max;          // end of bar or cell
    if( bar_or_cell == BAR ) {
      kf = cd->bar[cp->b+1];          // start of next bar
    }
    else {
      kf = cd->bar[cd->cell[cp->c+1]];// start of next cell
    }
  }
  else { // portion is FRONT, MIDDLE or FRINGE
    kf = k_min + random()%( k_max+1 - k_min );
    kf = end_of_instruction(cd,kf);
  }
  cue_clip( cd,cp,ki,kf,bar_or_cell,portion,parent );
}

/********************************************************//**
*/
void exclude_final_bar(
                       Code *cd,
                       Clip *cp,
                       int bar_or_cell,
                       int portion,
                       int parent
                      )
{
  if( portion == WHOLE || portion == BACK ) {
    if( parent == 1 ) {
      if( bar_or_cell == BAR ?
          (cp->k_stop == cd->bar[cp->b+1])
         :(cp->k_stop == cd->bar[cd->cell[cp->c+1]])) {
        cp->k_stop--;
        if( is_dot_or_digit(cd->codon[cp->k_stop])) {
          cp->k_stop--;
        }
      }
    }
    else { // ( parent == 0 )
      if( bar_or_cell == BAR ?
               (cp->k == cd->bar[cp->b+1])
             : (cp->k == cd->bar[cd->cell[cp->c+1]])) {
        cp->k--;
        if( is_dot_or_digit(cd->codon[cp->k])) {
          cp->k--;
        }
      }
    }
  }
}

/********************************************************//**
   Assume cell cp1->c has same number of bars in cd0
       as cell cp1->c in cd1.
   Choose subbar or subcell of cd1 that most closely matches
   the specified subbar or subcell in cd0.
*/
void matching_clip(
                   Code *cd0,
                   Clip *cp0,
                   Code *cd1,
                   Clip *cp1,
                   int bar_or_cell,
                   int portion
                  )
{
  int k0i,k0f,k1i,k1f;
  // k0i,k0f are cp0->k,k_stop in increasing order
  if( cp0->k_stop < cp0->k ) {
    k0i = cp0->k_stop;
    k0f = cp0->k;
  }
  else {
    k0i = cp0->k;
    k0f = cp0->k_stop;
  }
  //cp1->c = cp0->c;
  if( bar_or_cell == BAR ) {
    cp1->b = cp0->b + cd1->cell[cp1->c] - cd0->cell[cp0->c];
    if( portion == FRONT || portion == WHOLE ) {
      k1i = cd1->bar[cp1->b]; // start of bar
    }
    else { // portion is BACK, MIDDLE or FRINGE
#ifdef ALIGNED_CODON
      k1i = aligned_codon( cd0,cp0->b,k0i,cd1,cp1->b );
#else
      k1i = approximate_codon(cd0,cp0->b,k0i,cd1,cp1->b);
#endif
    }
    if( portion == WHOLE || portion == BACK ) {
      //if( k0f == cd0->bar[cp0->b+1] ) {
      k1f = cd1->bar[cp1->b+1]; // start of next bar
      //else {
      //k1f = last_codon_in_bar( cd1,cp1->b );
    }
    else {
#ifdef ALIGNED_CODON
      k1f = aligned_codon( cd0,cp0->b,k0f,cd1,cp1->b );
#else
      k1f = approximate_codon(cd0,cp0->b,k0f,cd1,cp1->b);
#endif
    }
  }
  else { // CELL
    int b0i,b0f,b1i,b1f;
    b0i = cd0->cell[cp0->c];
    if( !( portion == FRONT || portion == WHOLE )) {
      while( cd0->bar[b0i+1] <= k0i ) {
        b0i++;
      }
    }
    if( portion == WHOLE || portion == BACK ) {
      b0f = cd0->cell[cp0->c+1] - 1;
    }
    else {
      b0f = b0i;
      while( cd0->bar[b0f+1] <= k0f ) {
        b0f++;
      }
    }
    b1i = b0i + cd1->cell[cp1->c] - cd0->cell[cp0->c];
    b1f = b0f + cd1->cell[cp1->c] - cd0->cell[cp0->c];
    if( portion == FRONT || portion == WHOLE ) {
      k1i = cd1->bar[cd1->cell[cp1->c]]; // start of cell
    }
    else { // portion is BACK, MIDDLE or FRINGE
      //k1i = random_codon_in_bar( cd1,b1i );
#ifdef ALIGNED_CODON
      k1i = aligned_codon(cd0,b0i,k0i,cd1,b1i);
#else
      k1i = approximate_codon(cd0,b0i,k0i,cd1,b1i);
#endif
    }
    if( portion == WHOLE || portion == BACK ) {
      //if( k0f == cd0->bar[cd0->cell[cp0->c+1]] ) {
      k1f = cd1->bar[cd1->cell[cp1->c+1]]; // start of next cell
      //else {
      //k1f = last_codon_in_cell( cd1,cp1->c );
    }
    else {
      //k1f = random_codon_in_bar( cd1,b1f );
#ifdef ALIGNED_CODON
      k1f = aligned_codon(cd0,b0f,k0f,cd1,b1f);
#else
      k1f = approximate_codon(cd0,b0f,k0f,cd1,b1f);
#endif
    }
  }
  cue_clip( cd1,cp1,k1i,k1f,bar_or_cell,portion,1 );
}

/********************************************************//**
   Set cp->k,k_stop,Lb,Rb, based on ki,kf.
   If BAR, assume cp->b has already been set.
   If CELL,  set  cp->b to be the bar where copying begins.
   Set cp->Lb,Rb as appropriate.
   The portion (FRONT, MIDDLE, BACK, FRINGE or WHOLE) indicates
   the section of code copied from the PRIMARY parent.
   The complement will be copied from the SECONDARY parent
  (FRONT <-> BACK, MIDDLE <-> FRINGE, WHOLE <-> none).
   For all portions except FRINGE, stop_then_start is TRUE
   for the PRIMARY parent and FALSE for the SECONDARY parent,
   because code is copied from the PRIMARY, then SECONDARY,
   then PRIMARY again. When portion is FRINGE, stop_then_start
   is FALSE for the PRIMARY and TRUE for the SECONDARY parent.
   The clip will be:
   [cp->k_stop+1..cp->k] when stop_then_start is TRUE,
   [cp->k+1..cp->k_stop] when stop_then_start is FALSE.
*/
void cue_clip( Code *cd,Clip *cp,int ki,int kf,
               int bar_or_cell,int portion,int parent )
{
  int stop_then_start;
  if( ki > kf ) {        // ensure ki <= kf
    int k = ki;
    ki    = kf;
    kf    = k;
  }
  if( parent == PRIMARY ) {
    stop_then_start = ( portion != FRINGE );
  }
  else {                 // secondary parent
    stop_then_start = ( portion == FRINGE );
  }
  if( stop_then_start ) {//        ki     kf
    cp->k_stop = ki;     //   ------.......------
    cp->k      = kf;     //        k_stop k
  }
  else {                 //        ki     kf
    cp->k      = ki;     //   ......-------......
    cp->k_stop = kf;     //        k      k_stop
  }
  if( bar_or_cell == BAR ) {
    if( parent == SECONDARY ) {
      cp->Lb = cp->b;    // copying only (part of) one bar
      cp->Rb = cp->b+1;
    }
    else {               // preserve branch offsets
      cp->Lb = cd->cell[cp->c];
      cp->Rb = cd->cell[cp->c+1];
    }
  }
  else {       // CELL
    int b0,b1; // find bars b0,b1 containing ki,kf
    b0 = cd->cell[cp->c];
    while( cd->bar[b0+1] <= ki ) {
      b0++;
    }
    b1 = b0;
    while( b1 < cd->cell[cd->num_cells] && cd->bar[b1+1] <= kf ) {
      b1++;
    }
    if( stop_then_start ){//         ki   kf
      cp->Rb = b0;        //   -----+.......+-----
      cp->Lb = b1;        //         R     L   
    }
    else {                //         ki   kf
      cp->Lb = b0;        //   .....+-------+.....
      cp->Rb = b1+1;      //         L        R
    }
    cp->b = stop_then_start ? cd->cell[cp->c] : b0 ;
  }
}

/********************************************************//**
   Find the codon in bar b1 of cd1 most nearly aligned
   with codon k0 in bar b0 of cd0.
*/
int aligned_codon( Code *cd0, int b0, int k0,
                   Code *cd1, int b1 )
{
  int k0_min,k0_max;
  int k1_min,k1_max;
  int k1,k;

  // find start and end of bar b0,b1 in cd0,cd1
  k0_min = cd0->bar[b0];
  k0_max = last_codon_in_bar( cd0,b0 );
  k1_min = cd1->bar[b1];
  k1_max = last_codon_in_bar( cd1,b1 );
  // find k1 in cd1 aligned with k0 in cd0
  k = k0;
  while( k > k0_min && is_dot_or_digit(cd0->codon[k-1])) {
    k--;
  }
  k1 = k1_min + (int)floor((random_uniform()+0.5*(k+k0)-k0_min)
                   *(k1_max+1-k1_min)/(double)(k0_max+1-k0_min));
  //printf("[%d|%d|%d] -> [%d|%d|%d]\n",k0_min,k0,k0_max,k1_min,k1,k1_max);
  k1 = end_of_instruction( cd1,k1 );
  return( k1 );
}

/********************************************************//**
   Return index of the last codon in this instruction.
*/
int end_of_instruction( Code *cd, int k )
{
  while( is_dot_or_digit(cd->codon[k])) {
    k++;
  }
  // treat #- as a single instruction
  if(cd->codon[k] == PSH && cd->codon[k+1] == NEG) {
    k++;
  }
  return( k );
}

/********************************************************//**
   Return last codon in the specified bar (excluding . or 8)
*/
int last_codon_in_bar( Code *cd, int b )
{
  int k = cd->bar[b+1]-1;
  if( is_dot_or_digit( cd->codon[k] )) {
    k--;
  }
  return( k );
}

/********************************************************//**
   Return last codon in the specified cell (excluding 8)
*/
int last_codon_in_cell( Code *cd, int c )
{
  int k = cd->bar[cd->cell[c+1]] - 1;
  if( is_dot_or_digit( cd->codon[k] )) {
    k--;
  }
  return( k );
}

/********************************************************//**
   Find cell c containing bar b.
*/
int cell_for_bar( Code *cd, int b )
{
  int lo=0, hi = cd->num_cells-1;
  int mid;

  while( lo != hi ) {
    mid = ( lo + hi )/ 2;
    if( cd->cell[mid+1] <= b ) {
      lo = mid + 1;
    }
    else {
      hi = mid;
    }
  }

  return( lo );
}

/********************************************************//**
   Return the index of the first param after cd->codon[k].
*/
int get_param_index( Code *cd, int k )
{
  int lo=0, hi = cd->num_param;
  int mid;
  while( lo != hi ) {
    mid = ( lo + hi )/ 2;
    if( cd->param[mid] <= k ) {
      lo = mid + 1;
    }
    else {
      hi = mid;
    }
  }
  return( lo );
}

/********************************************************//**
   Choose r randomly between 0 and max, weighting max double
*/
int choose_index( int max )
{
  return( max - (random()%(max+2))%(max+1));
}

/********************************************************//**
   Choose a block randomly between 0 and num_blocks-1,
   weighted by similarity to block0.
*/
int similar_block( int block0, int num_blocks )
{
  int modulus,num_choices;
  int block,r,r_max=0;
  while( power2[r_max+1] <= num_blocks ) {
    r_max++;
  }
  // choose r between 0 and r_max, weighting 0 double
  r = (random()%( r_max+2 ))%( r_max+1 );
  modulus = power2[r];
  // how many numbers between 0 and num_blocks-1
  // are equal to block0 modulo power2[r]
  num_choices = 1+(num_blocks-1-(block0%modulus))/modulus;
  block =(random()% num_choices)*modulus + block0%modulus;
  //printf("%d %d\n",c0,c);
  return( block );
}

#ifndef ALIGNED_CODON
/********************************************************//**
   Return an integer within the range 0..max
   chosen by proximity to mid
*/
int approximate_int( double mid, int max )
{
  double n0,n1;
  n0 = random_uniform()* max;
  if( n0 > mid ) {
    n1 = mid + random_uniform()*( max - mid );
    if( n1 < n0 ) {
      n0 = n1;
    }
  }
  else { // n0 <= mid
    n1 = random_uniform()* mid;
    if( n1 > n0 ) {
      n0 = n1;
    }
  }
  return((int)floor(n0));
}

/********************************************************//**
   Find a bar in cell c1 of cd1 that is in close proximity
    to bar b0 in cell c0 of cd0.
*/
int approximate_bar( Code *cd0, int c0, int b0,
                     Code *cd1, int c1 )
{
  double b_mid;
  int num_bars = cd1->cell[c1+1] - cd1->cell[c1];

  b_mid = (0.5 + b0 - cd0->cell[c0])* num_bars/
                (double)(cd0->cell[c0+1] - cd0->cell[c0]);

  return( cd1->cell[c1] + approximate_int(b_mid,num_bars));
}

/********************************************************//**
   Find a codon in bar b1 of cd1 in close proximity
   to codon k0 in bar b0 of cd0.
*/
int approximate_codon( Code *cd0, int b0, int k0,
                       Code *cd1, int b1 )
{
  double k_mid;
  int k0_min,k0_max;
  int k1_min,k1_max,codon_range;
  int k,k1;

  // find start and end of bar b0,b1 in cd0,cd1
  k0_min = cd0->bar[b0];
  k0_max = last_codon_in_bar( cd0,b0 );
  k1_min = cd1->bar[b1];
  k1_max = last_codon_in_bar( cd1,b1 );
  codon_range = k1_max+1 - k1_min;
  // find k1 in cd1 aligned with k0 in cd0
  k = k0;
  while( k > k0_min && is_dot_or_digit(cd0->codon[k-1])) {
    k--;
  }
  k_mid = (0.5+0.5*(k+k0)-k0_min)*codon_range/
                         (double)(k0_max+1-k0_min);
  k1 = k1_min + approximate_int(k_mid,codon_range);
  k1 = end_of_instruction( cd1,k1 );

  return( k1 );
}

/********************************************************//**
   Choose subcell from cell cp->c of cd that is most nearly
   aligned with the specified subcell of cell cp0->c in cd0.
*/
void approximate_clip(
                      Code *cd0,
                      Clip *cp0,
                      Code *cd1,
                      Clip *cp1,
                      int bar_or_cell,
                      int portion
                     )
{
  int k0i,k0f,k1i,k1f;
  // k0i,k0f are cp0->k,k_stop in increasing order
  if( cp0->k_stop < cp0->k ) {
    k0i = cp0->k_stop;
    k0f = cp0->k;
  }
  else {
    k0i = cp0->k;
    k0f = cp0->k_stop;
  }
  if( bar_or_cell == BAR ) {
    cp1->b = approximate_bar(cd0,cp0->c,cp0->b,cd1,cp1->c);
    if( portion == FRONT || portion == WHOLE ) {
      k1i = cd1->bar[cp1->b]; // start of bar
    }
    else { // portion is BACK, MIDDLE or FRINGE
      k1i = approximate_codon( cd0,cp0->b,k0i,cd1,cp1->b );
    }
    if( portion == WHOLE || portion == BACK ) {
      //k1f = last_codon_in_bar( cd1,cp1->b );
      k1f = cd1->bar[cp1->b+1]; // start of next bar
    }
    else {
      k1f = approximate_codon( cd0,cp0->b,k0f,cd1,cp1->b );
    }
  }
  else { // CELL
    int b0i,b0f=0,b1i,b1f;
    b0i = cd0->cell[cp0->c];
    if( portion == FRONT || portion == WHOLE ) {
      b1i = cd1->cell[cp1->c]; // first bar in cell
    }
    else {
      // find bar b0i containing k0i
      while( cd0->bar[b0i+1] <= k0i ) {
        b0i++;
      }
      b1i = approximate_bar(cd0,cp0->c,b0i,cd1,cp1->c);
    }
    if( portion == WHOLE || portion == BACK ) {
      b1f = cd1->cell[cp1->c+1]-1;// last bar in cell
    }
    else {
      // find bar b0f containing k0f
      b0f = b0i;
      while( cd0->bar[b0f+1] <= k0f ) {
        b0f++;
      }
      b1f = approximate_bar(cd0,cp0->c,b0f,cd1,cp1->c);
    }
    if( b1i > b1f ) { // ensure b1i <= b1f
      int b = b1i;
      b1i   = b1f;
      b1f   = b;
    }
    if( portion == FRONT || portion == WHOLE ) {
      k1i = cd1->bar[cd1->cell[cp1->c]]; // start of cell
    }
    else { // portion is BACK, MIDDLE or FRINGE
      k1i = approximate_codon( cd0,b0i,k0i,cd1,b1i );
    }
    if( portion == WHOLE || portion == BACK ) {
      //k1f = last_codon_in_cell( cd1,cp1->c );
      k1f = cd1->bar[cd1->cell[cp1->c+1]]; // start of next cell
    }
    else {
      k1f = approximate_codon( cd0,b0f,k0f,cd1,b1f );
    }
  }
  cue_clip( cd1,cp1,k1i,k1f,bar_or_cell,portion,1 );
}
#endif

/********************************************************//**
   Choose a block randomly between 0 and num_blocks-1,
   weighted by similarity to block0, but not equal to block0.

int similar_not_equal_block( int block0, int num_blocks )
{
  int block;
  if( num_blocks == 1 ) {
    block = 0;
  }
  else {
    do {
      block = similar_block( block0, num_blocks );
    } while( block == block0 );
  }
  return( block );
}
*/

/********************************************************//**
   Replace cell c0 of cd0 with cell c1 from cd1 and,
   optionally, insert an instruction to jump to cell c0.

void jump_crossover(
                    Template *tp,
                    Code *cd0,
                    Code *cd1
                   )
{
  int Fc[MAX_LEVEL];
  int Nc[MAX_LEVEL];
  Clip clip0;
  Boolean chop_cell;
  int r_max;
  int b,c=0,c0,k0=0,k1=0;
#ifdef DEBUG
  print_line();
  di += sprintf( &ds[di]," JUMP CROSSOVER\n");
#endif
  if( random_bit()) {
    insert_branch( tp,cd0,&clip0,cd1 ); // combine with branch
    return;
  }

  chop_cell = random_bit();

  if( chop_cell ) {
             // choose cell c and subcell [k0+1..k1]
    c = choose_block_by_focus( cd0,1 );
  }
  else {     // choose insertion point k0 (=k1)
    r_max = adjust_focus( cd0,Fc,Nc );
    k0 = k1 = choose_codon( cd0,Fc,Nc,r_max,TRUE );
    while( cd0->bar[cd0->cell[c+1]] <= k0 ) {
      c++;   // find cell c containing codon k0
    }
    b = cd0->cell[c];
    while( cd0->bar[b+1] <= k0 ) {
      b++;   // find bar  b containing codon k0
    }
  }

  clip_cell(tp,cd0,&clip0,0);

  do {
    c0 = choose_jump_index( c,cd0->num_cells );
  } while( c0 == c );
  if( c0 < c ) {
    cross_cells( tp,cd0,&clip0,c0,cd1 );
  }

  // copy from cd0, up to start of cell c
  clip0.k_stop = cd0->bar[cd0->cell[c]];
  copy_verbatim(tp,cd0,&clip0);

  if( chop_cell ) {
    random_clip(cd0,&clip0,CELL,MIDDLE,0);
    k0 = clip0.k_stop; // start of subcell
    k1 = clip0.k;      //  end  of subcell
    tp->Rb = tp->b + clip0.Rb - clip0.b; //cd0->cell[c];
    tp->Lb = tp->Rb;
    tp->cell[c+1] -= clip0.Lb - clip0.Rb;
    clip0.k = cd0->bar[cd0->cell[c]];
  }

  // copy from start of cell to start of subcell
  //clip0.b = cd0->cell[c];
  //clip0.k = cd0->bar[clip0.b];
  clip0.k_stop = k0;
  copy_fix_digits( tp,cd0,&clip0 );

  add_focus_bar( tp,tp->b );

  insert_integer( tp,c0 );
  insert_operator(tp,'j');

  // copy from end of subcell to end of cell
  if( chop_cell ) {
    clip0.b = clip0.Lb;
    clip0.k = k1;
  }
  clip0.k_stop = cd0->bar[cd0->cell[c+1]];
  copy_fix_digits(tp,cd0,&clip0);

  // prepare to start copying from next cell
  //clip_cell(tp,cd0,&clip0,c+1);

  if( c0 > c ) {
    cross_cells( tp,cd0,&clip0,c0,cd1 );
  }

  // copy to end of code
  clip0.k_stop = cd0->last_codon;
  copy_verbatim(tp,cd0,&clip0);
#ifdef DEBUG
  //print_debug = TRUE;
  flush_debug();
#endif
}
*/
