/** \file eval.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "step.h"
#include "scan_print.h"
#include "cross_mutate.h"
#include "eval.h"

#define  LOG_LIBRARY   1

//#define  FAIL_BANK     1

// allowing trim mutations is always better
#define ALLOW_TRIM_MUTATIONS 1

//#define  ZERO_SEED       1

//#define CHOOSE_LEVEL_BY_PROPORTION 1

//#define ADAPTIVE_BREEDING  1


double sigmoid(double x);
Boolean anneal(double x);
double  leditd(int output_length, Floating output[],
               int target_length, Floating target[]);

void   add_to_library(Library *lib,Code *cd);
void      log_library(Library *lib,Code *cd);
Code   *mate_from_lib(Library *lib);
Code   *mate_from_lad(Ladder *lad,Level level);
Level    choose_level(Ladder *lad,Library *lib);

//int  cycle=0; // evolutionary cycle (when searching)

Long ncomp=0; // total number of fitness comparisons
Long neval=0; // total number of items evaluated

#ifdef MAX_TRACE
#define FAIL_TRACE  1
int  print_trace = FALSE;
#endif

char *level_code = (char *)"=SIPTKB:CJ12345";
char *status_code[10] = {
      (char *)"FLD",(char *)"SUP",(char *)"FST",(char *)"SRT",(char *)"EQL",
      (char *)"LNG",(char *)"SLW",(char *)"INF",(char *)"PNL",(char *)"RJT"};
char *status_name[10] = {
      (char *)" FAILED  ",(char *)"SUPERIOR ",(char *)" FASTER  ",
      (char *)" SHORTER ",(char *)"  EQUAL  ",(char *)" LONGER  ",
      (char *)" SLOWER  ",(char *)"INFERIOR ",(char *)" PENALTY ",(char *)" REJECT  "};
char *level_name[BLOCK+1]={
      (char *)" COPY ",(char *)" GRAD ",(char *)"INTERP",(char *)" TUNE ",
      (char *)" TRIM ",(char *)"POINT ",(char *)" BAR  ",(char *)"BRANCH",
      (char *)" CELL ",(char *)"BLOCK "};

long eval_seed;
#ifdef USE_MPI
int my_rank;
#endif

/********************************************************//**
   Return pointer to a new, empty candidate with a new, empty agent.
*/
Candidate *new_candidate()
{
  Candidate *can;
  can = (Candidate *)malloc(sizeof(Candidate));
  check_null(can);
  can->agt = new_agent();
  reset_score( can );
  return(can);
}

/********************************************************//**
   Compare b_can to a_can and set b_can->status to SUPERIOR,
   FASTER,SHORTER,EQUAL,LONGER,SLOWER,INFERIOR,PENALTY or REJECT.
   Return TRUE if b is better than a, FALSE if it is worse.
   If not final, return TRUE if b still has a chance of
   being better than a, FALSE otherwise.
*/
Boolean better_candidate(
                         Candidate *b_can,
                         Candidate *a_can,
                         double  target_cost,
                         Boolean final
                        )
{
  Score *a,*b;
  double a_steps=0.0,b_steps=0.0,a_score=0.0,b_score=0.0;
  int a_codons,b_codons;

  // first deal with the case where a or b is NULL
  if( a_can == NULL || b_can == NULL ) {
    if( b_can != NULL ) { // check for REJECT or PENALTY
      b = &b_can->score;
      if( b->reject ) {
        b_can->status = REJECT;
      }
      else if( b->penalty > 0.0 ) {
        b_can->status = PENALTY;
      }
      else {
        b_can->status = SUPERIOR;
      }
    }
    return( !final );  // can't replace NULL
  }

  a = &a_can->score;
  b = &b_can->score;

  // if result is final, the candidate completing a
  // larger number of trials is necessarily better
  if( final &&(   b->num_trials != a->num_trials
               || b->successful != a->successful )) {
    if(      b->num_trials  >  a->num_trials
       ||(   b->num_trials ==  a->num_trials
          && b->successful && !a->successful )) {
      b_can->status = SUPERIOR; return TRUE;
    }
    else { //( b->num_trials < a->num_trials ) {
      if( b->reject ) {
        b_can->status = REJECT;
      }
      else if( b->penalty > a->penalty ) {
        b_can->status = PENALTY;
      }
      else {
        b_can->status = INFERIOR;
      }
      return FALSE;
    }
  }

  // cases where a,b differ in reject or penalty
  if( a->reject || b->reject ) {
    if( a->reject && !b->reject ) {
      b_can->status = SUPERIOR; return TRUE;
    }
    if( b->reject && !a->reject ) {
      b_can->status = REJECT;   return FALSE;
    }
    if( a->reject && b->reject ) {
      if( b->penalty < a->penalty ) {
        b_can->status = SUPERIOR; return TRUE;
      }
      //( b->penalty > a->penalty ) handled below
    }
  }
  if( a->penalty != b->penalty ) {
    if( b->penalty < a->penalty ) {
      b_can->status = SUPERIOR;
      return( TRUE );
    }
    else { // ( b->penalty > a->penalty )
      b_can->status = PENALTY;  return FALSE;
    }
  }
  // now we can assume a,b are same in reject,penalty

  // check if a or b has achieved a perfect score
  if( final ) {
    Boolean a_perfect = perfect_score( a_can,target_cost );
    Boolean b_perfect = perfect_score( b_can,target_cost );
    if( a_perfect && !b_perfect ) {
      b_can->status = INFERIOR; return FALSE;
    }
    if( b_perfect && !a_perfect ) {
      b_can->status = SUPERIOR; return TRUE;
    }
  }

  a_codons = a_can->agt->cd->last_codon;
  b_codons = b_can->agt->cd->last_codon;
  //if(   b_can->agt->cd->mlevel != TRIM
  //   && a_codons > 16 && b_codons > 16 ) {
  a_score  = SPACE_FACTOR * a_codons;
  b_score  = SPACE_FACTOR * b_codons;
  //}

  if( !final ) {
    if( a->reject ) { // therefore also b->reject
      return( TRUE );
    }
    if( a->num_trials > 0 ) {
      a_steps = a->num_steps/(double)a->num_trials;
      b_steps = b->num_steps/(double)a->num_trials;
      a_score += TIME_FACTOR * a_steps;
      b_score += TIME_FACTOR * b_steps;
    }
    if(       b_score <= a_score
       ||(    b_can->agt->cd->level.m <= POINT
          &&( a->penalty > 0.0 || b->cost < a->cost ))) {
      return( TRUE );
    }
    else {
      return( FALSE );
    }
  }

  if( b->reject ) {
    b_can->status = (b_score == a_score) ?  EQUAL  :
                   ((b_score <  a_score) ? SHORTER : LONGER );
    return( anneal( b_score - a_score ));
    //return( b_score < a_score ? BETTER : WORSE );
  }
  // now we can assume ( !b->reject && !a->reject )
  if( a->num_trials > 0 ) { // num_trials MUST BE EQUAL????
    a_steps = a->num_steps/(double)a->num_trials;
    a_score += TIME_FACTOR * a_steps;
  }
  if( b->num_trials > 0 ) {
    b_steps = b->num_steps/(double)b->num_trials;
    b_score += TIME_FACTOR * b_steps;
  }
  //if( b->penalty == 0.0 && b->cost > target_cost ) {
  if( b->cost > target_cost ) {
    a_score += COST_FACTOR * a->cost;
    b_score += COST_FACTOR * b->cost;
  }
  if(    b->cost == a->cost
     ||( b->penalty == 0.0 && b->cost <= target_cost )) {
    if( b_steps == a_steps ) {// assume identical execution
      b_can->status = (b_score == a_score) ?  EQUAL  :
                     ((b_score <  a_score) ? SHORTER : LONGER );
    }
    else {                    // assume identical output
      b_can->status = (b_score == a_score) ?  EQUAL  :
                     ((b_score <  a_score) ?  FASTER : SLOWER );
    }
  }
  else {                      // different outputs
    b_can->status = (b_score == a_score) ?  EQUAL   :
                   ((b_score <  a_score) ? SUPERIOR : INFERIOR );
  }
  //if( final ) {
  if( b_score < a_score ||( b_score == a_score && random_bit())) {
    return( TRUE );
  }
    //}
  else if(    b_can->agt->cd->level.m <= POINT
          &&( b->penalty > 0.0 || b->cost < a->cost )) {
    // POINT and TUNE mutations can make the code slightly
    // longer, provided they strictly reduce the cost
    return( anneal( b_score - a_score ));
  }
  else {
    return( FALSE );
  }
}

/********************************************************//**
   Compute sigmoid function of x.
*/
double sigmoid( double x )
{
  if( x > 0.0 ) {
    return(       1.0/( 1.0 + exp(-x )));
  }
  else {
    return( 1.0 - 1.0/( 1.0 + exp( x )));
  }
}

/********************************************************//**
   Return TRUE or FALSE, according to an annealing function.
*/
Boolean anneal( double x )
{
  return( x < 0.0 || random_uniform() < sigmoid( -x ));
}

/********************************************************//**
   Free all the space occupied by candidate.
*/
void free_candidate( Candidate *can )
{
  if( can != NULL ) {
    if( can->agt != NULL ) {
      if( can->agt->cd != NULL ) {
        free_code( can->agt->cd );
      }
      compress_agent( can->agt );
      free( can->agt );
    }
    free( can );
  }
}

/********************************************************//**
   Reset num_trials, num_steps, reject, penalty and cost.
*/
void reset_score( Candidate *can )
{
  memset(&can->score,0,sizeof(Score));
  can->score.successful       = TRUE;
  can->score.all_outputs_same = TRUE;
  can->status = PENALTY;
  /*
  can->score.num_trials = 0;
  can->score.num_steps  = 0;
  can->score.mismatch   = 0;
  can->score.reject     = FALSE;
  can->score.penalty    = 0.0;
  can->score.cost       = 0.0;
  */
}

/********************************************************//**
   Print score, penalty or reject.
*/
void print_score( Candidate *can, FILE *fp )
{
  if( can->score.reject ) {
    fprintf( fp,"   R   ");
  }
  else if( can->score.penalty > 0.0 ) {
    fprintf(fp,"P=%1.f ",can->score.penalty);
  }
  else {
    fprintf(fp,"    ");
  }
  fprintf(fp,"%7.2f",can->score.cost);
}

/********************************************************//**
   Return TRUE if no penalty, non-reject and target cost achieved.
*/
int perfect_score( Candidate *can, double target_cost )
{
  return( !can->score.reject && can->score.penalty == 0.0
                             && can->score.cost <= target_cost );
}

/********************************************************//**
   Compare the actual and target output,
   and update the penalty, cost and mismatch appropriately. 
*/
void update_score(
                  Candidate *can,
                  Channel   *out,
                  Channel   *tgt,
                  int        cost_type
                 )
{
  Floating *output = &out->val[out->si];
  Floating *target = &tgt->val[tgt->si];
  int output_length = out->length[out->im+1];
  int target_length = tgt->length[tgt->im+1];
  int mismatch = 0;
  double cost  = 0.0;
  int k;
  if( cost_type == LED || cost_type == NONE ) {
    if( cost_type == LED ) { // Levenstein Edit Distance
      cost=leditd(output_length,output,target_length,target);
      mismatch = ( cost > 0.0 );
    }
    if( output_length < 0 && target_length >= 0 ) {
      can->score.penalty += 1.0;
      mismatch = 1;
    }
  }
  else if( output_length != target_length ) {
    can->score.penalty += (double)abs(output_length - target_length);
    mismatch = 1;
  }
  else {
    mismatch = 0; // ?????????????
    for( k = 0; k < output_length; k++ ) {
      if( !isfinite( output[k] )) {
        can->score.penalty += 1.0;
        mismatch = 1;
      }
      else {
        switch( cost_type ) {
          // GENERALIZE TO MULTIPLE OUTPUTS ??????????
        case KLD: // Kulback-Leibler Divergence (cross entropy)
          if( output[k] <= 0.0 ) {
            if( target[k] > 0.0 ) {
              can->score.penalty += 1.0;
              cost -= output[k];
              mismatch = 1;
            }
          }
          else if( output[k] >= 1.0 ) {
            if( target[k] < 1.0 ) {
              can->score.penalty += 1.0;
              cost += output[k] - 1.0;
              mismatch = 1;
            }
          }
          else {
            if( target[k] > 0.0 && target[k] < 1.0 ) {
              cost +=   target[k] *log(    target[k])
                  +(1.0-target[k])*log(1.0-target[k]);
            }
            cost  -=    target[k] *log(    output[k])
                  +(1.0-target[k])*log(1.0-output[k]);
            if(  ( target[k] > 0.5 && output[k] <= 0.5 )
               ||( target[k] < 0.5 && output[k] >= 0.5 )) {
              mismatch = 1;
            }
          }
          break;

        case KLD2: // KLD re-scaled between -1 and +1
          if( output[k] <= -1.0 ) {
            if( target[k] > -1.0 ) {
              can->score.penalty += 1.0;
              mismatch = 1;
            }
          }
          else if( output[k] >= 1.0 ) {
            if( target[k] < 1.0 ) {
              can->score.penalty += 1.0;
              mismatch = 1;
            }
          }
          else {
            if( target[k] > -1.0 && target[k] < 1.0 ) {
              cost += (0.5*(1.0+target[k]))*log(0.5*(1.0+target[k]))
                     +(0.5*(1.0-target[k]))*log(0.5*(1.0-target[k]));
            }
            cost -= (0.5*(1.0+target[k]))*log(0.5*(1.0+output[k]))
                   +(0.5*(1.0-target[k]))*log(0.5*(1.0-output[k]));
            if(  ( target[k] > 0.0 && output[k] <= 0.0 )
               ||( target[k] < 0.0 && output[k] >= 0.0 )) {
              mismatch = 1;
            }
          }
          break;

        case SGN: /* L2 but capped at 1,-1 */
          if( target[k] <= -1.0 ) {
            if( output[k] > -1.0 ) {
              cost += ( output[k] + 1.0 )
                     *( output[k] + 1.0 );
            }
          }
          else if( target[k] >= 1.0 ) {
            if( output[k] < 1.0 ) {
              cost += ( 1.0 - output[k] )
                     *( 1.0 - output[k] );
            }
          }
          else {
            cost += ( output[k] - target[k] )
                   *( output[k] - target[k] );
          }
          if(  ( target[k] > 0.0 && output[k] <= 0.0 )
             ||( target[k] < 0.0 && output[k] >= 0.0 )) {
            mismatch = 1;
          }
          break;

        case SQR:
          cost += ( output[k] - target[k] )
                 *( output[k] - target[k] );
          mismatch |= ( output[k] != target[k] );
          break;

        case LIN:
          cost += fabs( output[k] - target[k] );
          mismatch |= ( output[k] != target[k] );
          break;
        }
      }
    }
  }
#ifdef ASSERT
  assert( cost >= 0.0 );
#endif
  can->score.num_trials++;
  can->score.cost      += cost;
  can->score.mismatch  += mismatch;
}

/********************************************************//**
   Compute modified Levenshtein Edit Distance.
*/
double leditd( int output_length, Floating output[],
               int target_length, Floating target[] )
{
  double led0[4096];
  double led1[4096];
  int i,j;
    
  if( output_length <= 0 ) {
    return((double)abs(output_length-target_length));
  }
  //TODO ensure target_length <= 4096
  for( j=0; j <= target_length; j++ ) {
    led0[j] = (double) j;
  }
  for( i=1; i <= output_length; i++ ) {
    led1[0] = (double) i;
    for( j=1; j <= target_length; j++ ) {
      if( fabs(output[i-1]-target[j-1]) < 1000000.0 ) {
        led1[j] = led0[j-1] + 1.0
          - 1.0/(1.0+log(1.0+fabs(output[i-1]-target[j-1]))/log(4.0));
      }
      else {
        led1[j] = led0[j-1] + 1.0;
      }
      if( led1[j] > led0[j] + 1.0 ) {
        led1[j] = led0[j] + 1.0;
      }
      if( led1[j] > led1[j-1] + 1.0 ) {
        led1[j] = led1[j-1] + 1.0;
      }
    }
    for( j=0; j <= target_length; j++ ) {
      led0[j] = led1[j];
    }
  }
  return( led1[target_length] );
}

/********************************************************//**
   Return pointer to a new, empty library.
*/
Library *new_library( int max_code )
{
  Library *lib;
  lib = (Library *)malloc(sizeof(Library));
  check_null( lib );
  lib->code = NULL;
  lib->code_base = 0;
  lib->max_code = max_code;
  lib->num_code = 0;
  lib->index = -1;

  return( lib );
}

/********************************************************//**
   Open the named file and scan its code to library
*/
void scan_code_to_library( Library *lib, char *libname )
{
  FILE *fp;
  Code *cd;
  fp = fopen_check( libname,(char *)"r" );
  cd = scan_code( fp );
  while( cd != NULL ) {
    add_to_library( lib,cd );
    cd = scan_code( fp );
  }
  fclose( fp );
}

/********************************************************//**
  Insert new code into library at the specified index
 (replacing any code that was there previously)
  set lib->index to index if insertion was successful,
                    -1 otherwise.
*/
void insert_library( Library *lib, int index, Code *cd )
{
  if( index >= lib->max_code ) {
    free_code( cd );
    lib->index = -1;
  }
  else {
    if( lib->code == NULL ) {
      lib->code=(Code **)malloc(lib->max_code*sizeof(Code *));
      check_null( lib->code );
      memset(lib->code,0,lib->max_code*sizeof(Code *));
    }
    free_code( lib->code[index] );
    lib->code[index] = cd;
    if( index >= lib->num_code ) {
      lib->num_code = index+1;
    }
    lib->index = index;
  }
}

/********************************************************//**
   Print all the code in the library
*/
void print_library( Library *lib, FILE *fout )
{
  Code *cd;
  int index;
  if( lib != NULL && lib->code != NULL ) {
    for( index=0; index < lib->num_code; index++ ) {
      cd = lib->code[index];
      if( cd == NULL ) {
        fprintf(fout, "<null>\n");
      }
      else if( cd->num_cells == 1 ) {
        print_cell( cd,0,fout );
      }
      else {
        print_code( cd,fout );
      }
    }
  }
}

/********************************************************//**
   Free all the code is the library
*/
void clear_library( Library *lib )
{
  int l;
  if( lib != NULL ) {
    for( l=0; l < lib->num_code; l++ ) {
      free_code( lib->code[l] );
      lib->code[l] = NULL;
    }
    free( lib->code );
    lib->code = NULL;
    lib->num_code = 0;
    lib->index = -1;
  }
}

/********************************************************//**
  If max_index is zero, free code. Otherwise, free the
  next item in the library (using a circular ordering)
  and replace it with the new code.
  Return the index of the newly added code.
*/
void add_to_library( Library *lib, Code *cd )
{
  if( lib != NULL ) {
    int index = (lib->index+1)% lib->max_code;
    insert_library( lib,index,cd );
  }
}

/********************************************************//**
   Insert new code into library in a logarithmic fashion,
   and return the index of the newly inserted code.
*/
void log_library( Library *lib, Code *cd )
{
  if( lib == NULL ) {
    free_code( cd );
  }
  else {
    int index = 0;
    if( lib->code != NULL ) {
      while( index < lib->max_code && lib->code[index] != NULL ) {
        if(          lib->code[index]->last_codon == cd->last_codon
           && strcmp(lib->code[index]->codon,cd->codon) == 0 ) {
          free_code( cd ); // identical code already in library
          return;
        }
        index++;
      }
      index = 0;
      while(   index < lib->max_code-1
            && lib->code[index] != NULL
            &&(random_bit() || random_bit() || random_bit())) {
        index++;
      }
    }
    insert_library( lib,index,cd );
  }
}

/********************************************************//**
   Return pointer to a new, empty ladder, including a
   codebank with space for max_bank codes per level.
*/
Ladder * new_ladder( int max_bank )
{
  Ladder *lad;
  int g,m,s;
  lad = (Ladder *)malloc(sizeof(Ladder));
  check_null( lad );
  lad->max_bank   = max_bank; // max number of codes in bank, per level
  /*
    total_candidates[s][m][g] is the total number of candidates at step s,
    with mutation level m and transgenic level g.
    For m >= BAR, total_sup_inf[s][m][g] is the number of candidates that are
    either superior or initially inferior; total_superior[s][m][g] is the
    number that are (eventually) superior (i.e. SUPERIOR, FASTER or SHORTER)
  */
  memset(lad->bank,0,MAX_LEVEL*sizeof(Library *));

  lad->adapting = TRUE;
  for( m=0; m < MAX_LEVEL; m++ ) {
    //lad->bank[m] = new_library( max_bank );
    for( s=0; s < MAX_LEVEL; s++ ) {
      for( g=0; g < 3; g++ ) {
        lad->total_candidates[s][m][g]  = 1;
        lad->total_descendants[s][m][g] = 1;
        lad->total_sup_inf[s][m][g]     = 1;
        lad->total_superior[s][m][g]    = 1;
      }
    }
  }
#ifdef CHOOSE_LEVEL_BY_PROPORTION
  /*
  for( s = 0; s < MAX_LEVEL; s++ ) {
    for( g=0; g < 3; g++ ) {
      for( m=0; m < BAR; m++ ) {
        lad->total_descendants[s][m][g]    =   1;
      }
      lad->total_descendants[s][BAR][g]    =  16;
      lad->total_descendants[s][BRANCH][g] =  16;
      lad->total_descendants[s][CELL][g]   =  64;
      lad->total_descendants[s][BLOCK][g]  = 128;
      for( m = BLOCK+1; m < MAX_LEVEL; m++ ) {
        lad->total_descendants[s][m][g]=2*lad->total_descendants[s][m-1][g];
      }
    }
  }
  */
#endif
#ifdef PRINT_HISTOGRAM
  memset(lad->hist,0,MAX_LEVEL*MAX_LEVEL*3*10*17*sizeof(int));
#endif
  for( s=0; s < MAX_LEVEL; s++ ) {
    lad->can[s] = new_candidate();
  }
  lad->num_children[CHAMP] = 0; //  champ has no children yet
  lad->s = CHAMP-1;
  return( lad );
}

/********************************************************//**
   Print all candidates in the ladder, together with their
   level, size, time, penalty and cost.
*/
void print_ladder( Ladder *lad, FILE *fp )
{
  Candidate *can;
  Code  *cd;
  int s;

  for( s=0; s <= lad->s; s++ ) {
    can = lad->can[s];
    if( s <= lad->s ) {
      fprintf(fp,"->[");
    }
    else {
      fprintf(fp,"..[");
    }
    if( can != NULL ) {
      cd = can->agt->cd;
      if( cd != NULL ) {
        fprintf(fp,"%c,%d,",level_code[cd->level.m],cd->last_codon);
      }
      fprintf(fp,"%d/%d,",can->score.num_steps,can->score.num_trials);

      if( can->score.reject || can->score.penalty > 0.0 ) {
        fprintf(fp,"-,");
      }
      else {
        fprintf(fp,"%1.2f,",can->score.cost);
      }
      fprintf(fp,"%s",status_code[can->status]);
    }
    fprintf(fp,"]");
  }
  fprintf(fp,"\n");
}

/********************************************************//**
   Free the ladder, its codebank and all its candidates.
*/
void free_ladder( Ladder *lad )
{
  int m,s;
  reset_codebank( lad );
  for( s=0; s < MAX_LEVEL; s++ ) {
    free_candidate( lad->can[s] );
  }
  for( m=0; m < MAX_LEVEL; m++ ) {
    if( lad->bank[m] != NULL ) {
      free( lad->bank[m] );
      lad->bank[m] = NULL;
    }
  }
  free( lad );
}

/********************************************************//**
   Reset ladder to contain only one candidate,
   with the specified code.
*/
void reset_ladder( Ladder *lad, Code *cd )
{
  int s;
  for( s=0; s < MAX_LEVEL; s++ ) {
    if( lad->can[s] == NULL ) {
      lad->can[s] = new_candidate();
    }
    else {
      reset_score( lad->can[s] );
      //lad->can[s]->status = PENALTY;
      if( lad->can[s]->agt->cd != NULL ) {
        free_code( lad->can[s]->agt->cd );
        lad->can[s]->agt->cd = NULL;
      }
    }
  }
  lad->can[CHAMP]->agt->cd = cd;
#ifdef RESERVE_CHAMP
  lad->can[0]->agt->cd = cross_mutate( cd,NULL,COPY,FALSE );
#endif
  lad->s = CHAMP;
#ifdef MAX_TRACE
  if( globals.verbose ) {
    lad->ti[CHAMP] = sprint_code( cd,lad->ts );
    //lad->ts + sprintf(lad->ts," %s\n",cd->codon);
  }
#endif
}

/********************************************************//**
   Free the contents and all space occupied by the codebank.
*/
void reset_codebank( Ladder *lad )
{
  int m;
  for( m=0; m < MAX_LEVEL; m++ ) {
    clear_library( lad->bank[m] );
  }
}

/********************************************************//**
   Free the contents and all space occupied by the codebank.
*/
void random_codebank( Ladder *lad, int num_code )
{
  Code *cd0=NULL;
  Code *cd;
  int index;
  int m;
  if( num_code > lad->max_bank ) {
    num_code = lad->max_bank;
  }
  if( lad->can[CHAMP]->agt != NULL ) {
    cd0 = lad->can[CHAMP]->agt->cd;
  }
  for( m=POINT; m <= BLOCK; m++ ) {
    if( lad->bank[m] == NULL ) {
      lad->bank[m] = new_library( lad->max_bank );
    }
    for( index = 0; index < num_code; index++ ) {
      cd = random_code(cd0,1);
      //print_code(cd,stdout);
      insert_library( lad->bank[m],index,cd );
    }
  }
}

/********************************************************//**
*/
void procreate( Ladder *lad, Library *lib )
{
  Code *cd1;
  Level level = choose_level( lad,lib );
  if( level.g == TRANSGENIC ) {
    cd1 = mate_from_lib( lib );
  }
  else {
    //cd1 = mate_from_lad( lad,POINT );
    cd1 = mate_from_lad( lad,level );
  }
  breed( lad,cd1,level );
}

/********************************************************//**
   Choose code (uniformly) randomly from library.
*/
Code *mate_from_lib( Library *lib )
{
  // hacky solution to ensure that empty slots
  // are not chosen (this problem comes up with multiple processes)
  int l;
  
  int a,b; // bounds
  if (   lib->code_base == 0
      ||( my_random()%2 == 0 && (lib->num_code != lib->code_base))) {
    // choose from shared code
    a = lib->code_base;
    b = lib->num_code;
  } else {
    // choose from code base
    a = 0;
    b = lib->code_base;
  }
  
  l = my_random()%(b-a);
  while (lib->code[a+l] == NULL) {
    l = (l + 1) % (b-a);
  }
  
  return( lib->code[a+l] );
}

/********************************************************//**
  Choose breeding partner from a (uniformly) randomly
  chosen level within the codebank of this ladder.
  If there is no codebank, or the codebank contains no
  code at the chosen level, return the top candidate
  from the ladder.
*/
Code *mate_from_lad( Ladder *lad,Level level )
{
  int option,num_options=0;
  int m;

  for( m = 0; m < MAX_LEVEL; m++ ) {
    if( lad->bank[m] != NULL && lad->bank[m]->num_code > 0 ) {
      num_options++;
    }
  }

  if( num_options > 0 ) {
    option = random()% num_options;
    for( m = 0; m < MAX_LEVEL; m++ ) {
      if( lad->bank[m] != NULL && lad->bank[m]->num_code > 0 ) {
        option--;
        if( option < 0 ) {
          return( mate_from_lib( lad->bank[m] ));
        }
      }
    }
  }
  // if codebank is empty, return current candidate
  return( top(lad)->agt->cd );
}

/********************************************************//**
   Choose mutation level randomly, weighted by the reciprocal
   of the expected number of offspring (descendants).
*/
Level choose_level( Ladder *lad, Library *lib )
{
#ifdef CHOOSE_LEVEL_BY_PROPORTION
  static double proportion[MAX_LEVEL]={0.0};
#endif
  Candidate *top;
  Code  *cd;
  Level  level;
  double x[MAX_LEVEL][3]={{0.0}};
  double portion;
  double z;
  int mlevel_max;
  int s,m,g; // step, mutation level

#ifdef CHOOSE_LEVEL_BY_PROPORTION
  if( proportion[POINT] == 0.0 ) {
    proportion[TUNE]   = 0.5;
    proportion[TRIM]   = 0.5;
    proportion[POINT]  = 2.0;
    proportion[BAR]    = 0.5;
    proportion[BRANCH] = 0.2;
    proportion[CELL]   = 0.2;
    proportion[BLOCK]  = 0.032;
    for( m = BLOCK+1; m < MAX_LEVEL; m++ ) {
      proportion[m] = 0.4*proportion[m-1];
    }
  }
#endif
  top = lad->can[lad->s];
  cd  = top->agt->cd;
  s = lad->s+1; // new candidate will be at step s
  /*
  switch( cd->level.m ) {
   case BLOCK:  mlevel_max = CELL;          break;
   case JUMP: case BRANCH: // case CELL:
                mlevel_max = BAR;           break;
   case POINT:  mlevel_max = TUNE;          break;
   case BAR:    mlevel_max = POINT;         break;
   default:     mlevel_max = cd->level.m-1; break;
  }
  */
  mlevel_max = cd->level.m-1;
  if( cd->num_param == 0 ) {
    x[TUNE][NON_ALIGNED] = 0.0; // no params to tune
  }
  else {
#ifdef CHOOSE_LEVEL_BY_PROPORTION
    x[TUNE][NON_ALIGNED] = proportion[TUNE];
#else
    x[TUNE][NON_ALIGNED] = 1.0;
#endif
  }
#ifdef NO_TUNING
  if( lad->s > 0 ) {
    x[TUNE][NON_ALIGNED] = 0.0;
  }
#endif
  x[TRIM][NON_ALIGNED] = x[TUNE][NON_ALIGNED];
#ifdef ALLOW_TRIM_MUTATIONS
  if( !top->score.reject && top->score.penalty == 0.0 ) {
#ifdef CHOOSE_LEVEL_BY_PROPORTION
    x[TRIM][NON_ALIGNED]  = x[TUNE][NON_ALIGNED] + proportion[TRIM];
#else
    x[TRIM][NON_ALIGNED]  = x[TUNE][NON_ALIGNED] + 1.0;
#endif
  }
#endif
#ifdef CHOOSE_LEVEL_BY_PROPORTION
  x[POINT][NON_ALIGNED] = x[TRIM][NON_ALIGNED] + proportion[POINT];
#else
  x[POINT][NON_ALIGNED] = x[TRIM][NON_ALIGNED] + 1.0;
#endif
  x[POINT][TRANSGENIC] = x[POINT][NON_ALIGNED];
  // TUNE, TRIM, POINT can only be NON_ALIGNED
  for( m = POINT+1; m <= mlevel_max; m++ ) {
    x[m][ALIGNED] = x[m-1][TRANSGENIC];
    for( g = ALIGNED; g <= TRANSGENIC; g++ ) {
      if(  ( g == ALIGNED     && m == BRANCH ) //|| m == JUMP ))
         ||( g == NON_ALIGNED && m >= BLOCK && power2[m-BLOCK] >= cd->num_cells)
         ||( g == TRANSGENIC  &&( lib == NULL || lib->num_code == 0 ))
           // ||( m == JUMP       && cd->num_cells == 1 )
         ||( top->status == REJECT && m == BRANCH )) {
              // || top->status == PENALTY )
              // &&( m == BRANCH || m == JUMP ))) {
        portion = 0.0;
      }
      else {
        //        portion = proportion[m]  *  lad->total_candidates[s][m][g]
        ///(double)(lad->total_candidates[s][m][g]+lad->total_descendants[s][m][g]);
        portion=sqrt(lad->total_candidates[s][m][g] /(double)(
                     lad->total_candidates[s][m][g]+lad->total_descendants[s][m][g]));
#ifdef CHOOSE_LEVEL_BY_PROPORTION
        portion *= proportion[m];
#endif
      }
      if( g == 0 ) {
        x[m][g] = x[m-1][2] + portion;
      }
      else {
        x[m][g] = x[m][g-1] + portion;
      }
    }
  }
  z = random_uniform() * x[mlevel_max][2];
  m = TUNE - 1;
  g = 2;
  while( x[m][g] <= z ) {
    m++;
    g = 0;
    while( g < 2 && x[m][g] <= z ) {
      g++;
    }
  }
  level.m = m;
  level.g = g;
  return( level );
}

/********************************************************//**
   move code at step s to codebank[m], where m is its
   mutation level, and set code to NULL.
*/
void move_to_codebank( Ladder *lad, int s )
{
  Code *cd;
  int m;

  // if identical code remains in ladder, free this code
  if( s >= CHAMP && strcmp(top(lad)->agt->cd->codon,
                           pop(lad)->agt->cd->codon) == 0 ) {
    free_code( lad->can[s]->agt->cd );
    lad->can[s]->agt->cd = NULL;
    return;
  }
  if( s < CHAMP ) { // copy champ to special level of codebank
    Level level_copy={COPY,NON_ALIGNED};
    cd = cross_mutate(lad->can[CHAMP]->agt->cd,NULL,level_copy);
     m = lad->can[CHAMP]->agt->cd->level.m+1;
  }
  else {
    cd = lad->can[s]->agt->cd; // move top or pop code to codebank
    m = top(lad)->agt->cd->level.m;
    lad->can[s]->agt->cd = NULL;
  }
  if(    m <  POINT || lad->max_bank == 0
     ||( s == CHAMP && lad->can[s+1]->agt->cd->level.m < POINT )) {
    free_code( cd );
  }
  else {
    if( lad->bank[m] == NULL ) {
      lad->bank[m] = new_library( lad->max_bank );
    }
#ifdef LOG_LIBRARY
    log_library( lad->bank[m],cd );
#else
    add_to_library( lad->bank[m],cd );
#endif
  }
}

/********************************************************//**
<pre>
   +-----+  +-----+  +-----+
 ->| pop |->| top |->| lop |
   +-----+  +-----+  +-----+
      |        |        x
      |        |
   +-----+  +-----+  +-----+
 ->|     |->| pop |->| top |
   +-----+  +-----+  +-----+
</pre>
   Let the current top become the pop, and generate a new
   top by mutating the current top, using cd1 and mlevel.
*/
void breed( Ladder *lad, Code *cd1, Level level )
{
  Candidate *top = lad->can[lad->s];
  Code *cd,*cd0;
  int g,m,s;              // mutation level, step

  cd0 = top->agt->cd;      // primary parent
  if( level.m >= BLOCK ) {
    // block must not exceed size of secondary parent
    while( cd1->num_cells < power2[level.m-BLOCK] ) {
      level.m--;
    }
  }

  cd = cross_mutate( cd0,cd1,level );
  m = cd->level.m; // cross_mutate() can change mutation level
  g = cd->level.g; // cross_mutate() can change kinship level
  //printf("m=%d,g=%d\n",m,g);
  lad->num_children[lad->s]++; // one more child for candidate at step s
  lad->s++;
  s = lad->s;                     // new candidate will have step s
  lad->num_children[s]    = 0;    // new candidate has no children yet
  lad->num_descendants[s] = 0;    // self not counted as a descendant
  free_code(lad->can[s]->agt->cd);// free previous code at step s
  lad->can[s]->agt->cd = cd;      //   install new code at step s
  reset_score( lad->can[s] );
#ifdef MAX_TRACE
  if( globals.verbose ) {
    if( lad->ti[s-1] + 1 + cd->last_codon < lad->ts + MAX_TRACE ) {
      lad->ti[s] = lad->ti[s-1] + sprintf(lad->ti[s-1],"%c",
        (level.g == TRANSGENIC) ? level_code[m] : level_code[m]+32);
      lad->ti[s] = lad->ti[s] + sprintf(lad->ti[s],"%d",level.g);
      lad->ti[s] = sprint_code( cd,lad->ti[s] );
      //lad->ti[s] = lad->ti[s] + sprintf(lad->ti[s],"%s",cd->codon);
    }
  }
#endif
  if( m <= POINT ) {
    lad->max_children[s] = 0;
  }
  else {
    lad->max_children[s] = lad->total_sup_inf[s][m][g]
                          /lad->total_superior[s][m][g];
    if( lad->max_children[s] > MAX_CHILD ) {
      lad->max_children[s] = MAX_CHILD;
    }
    else {
      int min_max_children = 1;
      if( s - CHAMP <= 3 ) {
        //min_max_children = power2[11 - 3*(s-CHAMP)];
        //min_max_children = power2[8 - 2*(s-CHAMP)];
        min_max_children = power2[7 - 2*(s-CHAMP)];
      }
      //if( g != ALIGNED ) {
      //  min_max_children *= 2;
      //}
      if( lad->max_children[s] < min_max_children ) {
        lad->max_children[s] = min_max_children;
      }
    }
  }
}

/********************************************************//**
   Return the top item from the ladder.
*/
Candidate *top( Ladder *lad )
{
  if( lad->s >= CHAMP ) {
    return( lad->can[lad->s] );
  }
  else {
    return( NULL );
  }
}

/********************************************************//**
   Return the second top item from the ladder.
*/
Candidate *pop( Ladder *lad )
{
  if( lad->s > CHAMP ) {
    return( lad->can[lad->s-1] );
  }
  else {
    return( NULL );
  }
}

/********************************************************//**
<pre>
   +-----+  +-----+  +-----+  +-----+
 ->|     |->| pop |->| top |->| lop |
   +-----+  +-----+  +-----+  +-----+
      |        |        |        |
      |        |        |        x
   +-----+  +-----+  +-----+
 ->| pop |->| top |->| lop |
   +-----+  +-----+  +-----+
</pre>
   Top item can only be removed if it is rejected,
   or if it has a penalty (when pop has no penalty)
   or if it has exceeded max_children.
   Remove top from ladder, so that pop becomes new top,
   top becomes new lop, and old lop is deleted.
   Return TRUE if top item was culled, FALSE otherwise.
*/
int cull_top( Ladder *lad, FILE *fo )
{
  Candidate *top,*pop;
  int top_m,top_g;
  int s=lad->s;
#ifdef PRINT_HISTOGRAM
  int top_status = lad->can[s]->status;
#endif

  if( s < CHAMP ) { // ladder is empty
    return FALSE;
  }
  else if( lad->can[s]->agt->cd == NULL ) {
    free_code( lad->can[s+1]->agt->cd );
    lad->can[s+1]->agt->cd = NULL;
    s--;
    return TRUE;
  }
  else if( s == CHAMP ) {
    return FALSE;
  }
  top   = lad->can[s];
  pop   = lad->can[s-1];
  top_m = top->agt->cd->level.m;
  top_g = top->agt->cd->level.g;

  switch( top->status ) {
  case FAILED: case SUPERIOR: case FASTER: case SHORTER:
    // debugging purposes; these cases should not occur
    fprintf(stderr,"top->status == %s\n",status_code[top->status]);
    exit(1);
  }

  // if top has exceeded max_children or is otherwise unfit,
  // change its status to FAILED so it can actually be culled
  if(      top->score.reject
     ||(   top->score.penalty  > 0.0
        && pop->score.penalty == 0.0 )
     ||  lad->num_children[s] >= lad->max_children[s] ) {
    switch( top->status ) {
     case EQUAL: case LONGER: case SLOWER: case INFERIOR:
       top->status = FAILED;
       break;
    }
  }

  switch( top->status ) {
  case REJECT: case PENALTY: case FAILED:
  case EQUAL:  case LONGER:  case SLOWER:
    // in these cases, top is actually culled and TRUE is returned
    if( lad->adapting && !pop->score.reject
                      &&  pop->score.penalty == 0 ) {
      // update statistics only if non-reject and no penalty
      int m = top_m, g = top_g;
      lad->total_candidates[s][m][g]++;
      lad->total_descendants[s][m][g] += lad->num_descendants[s];
      if( lad->num_children[s] > 0 ) {
        lad->total_sup_inf[s][m][g]++;
      }
      lad->num_descendants[s-1] += 1 + lad->num_descendants[s];
#ifdef PRINT_HISTOGRAM
      {
        // maintain histogram of log of the number of children
        int r=0;
        while( r < 16 && lad->num_children[s] >= power2[r] ) {
          r++;
        }
        lad->hist[s][top_m][top_g][top_status][r]++;
        //printf("hist[%d][%d][%d][%d]++\n",s-1,top_m,top->status,r);
      }
#endif
    }

#ifdef MAX_TRACE
#ifdef FAIL_TRACE
    if( globals.verbose ) {
      if( lad->ti[s] + 4 < lad->ts + MAX_TRACE ) {
        lad->ti[s-1] = lad->ti[s]
          + sprintf(lad->ti[s],"%s\n",status_code[top->status]);
      }
    }
#endif
#endif
    free_code( lad->can[s+1]->agt->cd ); // previous lop
    lad->can[s+1]->agt->cd = NULL;
    if( !top->score.reject && top->score.penalty <= pop->score.penalty ){
#ifdef RESERVE_CHAMP
      if( lad->s > CHAMP+1 ) {
        move_to_codebank( lad,lad->s );
      }
#else
      move_to_codebank( lad,lad->s );
#endif
    }
    lad->s--;
#ifdef MAX_TRACE
#ifdef FAIL_TRACE
    if( globals.verbose ) {
      if( lad->s == CHAMP ) {
        if( fo != NULL ) {
          sprintf(lad->ti[CHAMP],"%c",'\0');
          fprintf(fo,"%s\n",lad->ts);
        }
        lad->ti[CHAMP]=lad->ts+sprintf(lad->ts," %s\n",pop->agt->cd->codon);
      }
    }
#endif
#endif
    return TRUE;
    break;

  default: // top is not (yet) eligible to be culled
#ifdef MAX_TRACE
    if( globals.verbose ) {
      if( lad->num_children[s] == 0 ) {
        if( lad->ti[s] + 4 < lad->ts + MAX_TRACE ) {
            lad->ti[s] +=
             sprintf(lad->ti[s],"%s\n",status_code[top->status]);
        }
      }
    }
#endif
    return FALSE;
    break;
  }
}

/********************************************************//**
<pre>
   +-----+  +-----+  +-----+  +-----+
 ->|     |->| pop |->| top |->| lop |
   +-----+  +-----+  +-----+  +-----+
      |            \/            |
      |            /\            x
   +-----+  +-----+  +-----+
 ->| pop |->| top |->| lop |
   +-----+  +-----+  +-----+
</pre>
   First, mutation and kinship level of top and pop are swapped.
   Then, pop is culled from the ladder,
   by swapping it with top, so that top moves up one rung,
   and pop becomes the new lop (the old lop is deleted).
*/
void top_replace_pop( Ladder *lad, FILE *fo )
{
  if( lad->s > CHAMP ) {
    Candidate *top,*pop;
    int top_m,pop_m; // mutation level
    int top_g,pop_g; // kinship  level
    int n,s=lad->s;
    free_code( lad->can[s+1]->agt->cd ); // previous lop
    lad->can[s+1]->agt->cd = NULL;
    top = lad->can[s];
    pop = lad->can[s-1];

    top_m = top->agt->cd->level.m;
    pop_m = pop->agt->cd->level.m;
    top_g = top->agt->cd->level.g;
    pop_g = pop->agt->cd->level.g;

    //if( s > CHAMP ) {
    // top inherits focus cell(s) from pop
    top->agt->cd->Fc = pop->agt->cd->Fc;
    top->agt->cd->Nc = pop->agt->cd->Nc;

    // cull focus bars from levels less than pop_m
    n=0;
    while( n < top->agt->cd->num_focus_bars
            && top->agt->cd->Fb[n].m >= pop_m ) {
      n++;
    }
    top->agt->cd->num_focus_bars = n;
      //}
    if( lad->adapting && !pop->score.reject
                      &&  pop->score.penalty == 0 ) {
      int m = top_m;
      int g = top_g;

      lad->total_candidates[s][m][g]++;
      lad->total_descendants[s][m][g] += lad->num_descendants[s];
      if( top->status == SUPERIOR ) {
        lad->total_superior[s][m][g]++;
        lad->total_sup_inf[s][m][g]++;
      }
      else if( lad->num_children[s] > 0 ) {
        lad->total_sup_inf[s][m][g]++;
      }
      lad->num_descendants[s-1] += 1 + lad->num_descendants[s];
#ifdef PRINT_HISTOGRAM
      {
        // maintain histogram of log of the number of children
        int r=0;
        while( r < 16 && lad->num_children[s] >= power2[r] ) {
          r++;
        }
        lad->hist[s][top_m][top_g][top->status][r]++;
        //printf("hist[%d][%d][%d][%d]++\n",s-1,top_m,top->status,r);
      }
#endif
    }
    //switch( top->status ) {
    //case SUPERIOR: case SLOWER: case INFERIOR: case PENALTY: case REJECT:
#ifdef RESERVE_CHAMP
    if( s-1 > CHAMP ) {
      move_to_codebank( lad, s-1 );
    }
#else
    move_to_codebank( lad, s-1 );
#endif

    // top inherits mutation and kinship level from pop
    top->agt->cd->level.m = pop_m;
    top->agt->cd->level.g = pop_g;

    lad->can[s]   = pop;
    lad->can[s-1] = top;

    //if( s > CHAMP ) {
    lad->s--;

#ifdef MAX_TRACE
    if( globals.verbose ) {
      if( lad->ti[s] + 4 < lad->ts + MAX_TRACE ) {
        lad->ti[s-1] = lad->ti[s]
          + sprintf(lad->ti[s],"%s\n",status_code[top->status]);
      }
      if( lad->s == CHAMP ) {
        if( fo != NULL ) {
          sprintf(lad->ti[CHAMP],"%c",'\0');
#ifdef FAIL_TRACE
          fprintf(fo,"%s\n",lad->ts);
#else
          if(  (  pop->score.reject && !top->score.reject )
             ||( !top->score.reject &&  top->score.penalty < pop->score.penalty)
             ||( !top->score.reject &&  top->score.penalty == 0
               && top->score.cost <= pop->score.cost - 0.5 )) {
            fprintf(fo,"%s\n",lad->ts);
          }
#endif
        }
        lad->ti[CHAMP] = sprint_code( top->agt->cd,lad->ts );
      }
    }
#endif
  }
}

/********************************************************//**
   Compress all agents in the ladder.
*/
void compress_ladder( Ladder *lad )
{
  int s;
  for( s=0; s < MAX_LEVEL; s++ ) {
    if( lad->can[s] != NULL ) {
      compress_agent( lad->can[s]->agt );
    }
  }
}

#ifdef PRINT_HISTOGRAM
/********************************************************//**
   Print a histogram of how many candidates with what number
   of children were in which category.
*/
void print_hist( Ladder *lad,FILE *fo )
{
  int max_r,sum;
  int g,m,r,s,t;
  for( s = CHAMP+1; s < MAX_LEVEL; s++ ) {
    for( m = TUNE; m < MAX_LEVEL; m++ ) {
      for( g=0; g < 3; g++ ) {
        max_r = -1;
        for( r=0; r < 17; r++ ) {
          sum = 0;
          for( t=0; t < 10; t++ ) {
            sum += lad->hist[s][m][g][t][r];
          }
          if( sum > 0 ) {
            max_r = r;
          }
        }
        if( max_r >= 0 ) {
          fprintf(fo,"---------+-");
          for( r=0; r <= max_r; r++ ) {
            fprintf(fo,"-----");
          }
          fprintf(fo,"\n");
          if( m < BLOCK ) {
            fprintf(fo,"  %s", level_name[m] );
          }
          else {
            fprintf(fo," BLOCK%2d", power2[m-BLOCK]);
          }
          fprintf(fo," | ");
          for( r=0; r <= max_r; r++ ) {
            fprintf(fo,"%5d",power2[r]-1);
          }
          fprintf(fo,"\n---------+-");
          for( r=0; r <= max_r; r++ ) {
            fprintf(fo,"-----");
          }
          fprintf(fo,"\n");
          for( t=1; t < 10; t++ ) {
            fprintf(fo,"%s|",status_name[t]);
            fprintf(fo,"%6d",lad->hist[s][m][g][t][0]);
            for( r=1; r <= max_r; r++ ) {
              fprintf(fo,"%5d",lad->hist[s][m][g][t][r]);
            }
            fprintf(fo,"\n");
          }
        }
      }
    }
  }
  fprintf(fo,"---------+----------------------------");
  fprintf(fo,"---------------------------------------\n");
}
#endif

/********************************************************//**
   Print relevant statistics for each mutation level.
*/
void print_stats( Ladder *lad,FILE *fo )
{
  int g,m,s;

  for( s = CHAMP+1; s < MAX_LEVEL; s++ ) {
    for( m = POINT+1; m < MAX_LEVEL; m++ ) {
      for( g=0; g < 3; g++ ) {
        if( lad->total_candidates[s][m][g] > 1 ) {
          if( m < BLOCK ) {
            fprintf(fo,"%2d %s ", s-CHAMP, level_name[m] );
          }
          else {
            fprintf(fo,"%2d BLOCK%2d", s-CHAMP, power2[m-BLOCK]);
          }
          fprintf(fo,"(%6lld/%8lld)=%1.3f",
              lad->total_candidates[s][m][g],lad->total_descendants[s][m][g],
              lad->total_candidates[s][m][g]
    /(double)(lad->total_candidates[s][m][g]+lad->total_descendants[s][m][g]));

          fprintf(fo," |(%6d/%5d)=%4d\n",
              lad->total_sup_inf[s][m][g],lad->total_superior[s][m][g],
              lad->total_sup_inf[s][m][g]/lad->total_superior[s][m][g]);
        }
      }
    }
  }
}

/********************************************************//**
   Generate num_seeds random longs and store them into seed[]
*/
void generate_seeds(
                    int  max_seeds,
                    long seed[]
                   )
{
  int r;
  for( r=1; r <= max_seeds; r++ ) {
    seed[r] = my_random();
  }
#ifdef ZERO_SEED
  seed[1] = 0;
#endif
}

/********************************************************//**
   Shuffle seeds into a random order
*/
void shuffle_seeds(
                   int  max_seeds,
                   long seed[]
                  )
{
  int  *rank;
  long *seed0;
  int r,n;

  seed0 = (long *)malloc((max_seeds+1)*sizeof(long));
  check_null( seed0 );
  memcpy( seed0,seed,(max_seeds+1)*sizeof(long));

  rank = (int *)malloc((max_seeds+1)*sizeof(int));
  check_null( rank );
  memset( rank,0,(max_seeds+1)*sizeof(int));

  for( n=1; n <= max_seeds; n++ ) {
    r = 1 + random()% max_seeds;
    while( rank[r] != 0 ) {
      r = 1 + random()% max_seeds;
    }
    rank[r] = n;
  }
  for( r=1; r <= max_seeds; r++ ) {
    seed[r] = seed0[rank[r]];
  }
  free( rank );
  free( seed0 );
}

#ifdef RESERVE_CHAMP
/********************************************************//**
*/
void lop_replace_reserve( Ladder *lad )
{
  Candidate *can = lad->can[2];
  Code *cd = can->agt->cd;

  cd->level.m = lad->can[CHAMP]->agt->cd->level.m;
  cd->level.g = ALIGNED;
  cd->Fc = 0;
  cd->Nc = cd->num_cells;
  cd->num_focus_bars = 0;

  move_to_codebank( lad, 0 );

  lad->can[2] = lad->can[0];
  lad->can[0] = can;
}

/********************************************************//**
*/
void swap_reserve_champ( Ladder *lad )
{
  Candidate *can = lad->can[1];
  lad->can[1] = lad->can[0];
  lad->can[0] = can;
}
#endif

/********************************************************//**
   Return the item most recently removed from the ladder.
*/
Candidate *lop( Ladder *lad )
{
  return( lad->can[lad->s+1] );
}

/********************************************************//**
*/
void print_termination(Ladder *lad,
                       Long    ncomp,
                       Long    neval,
                       int     epoch,
                       int     task,
                       FILE   *fo
                      )
{
  fprintf(fo, "Task %d terminating early\n", task);
  fprintf(fo, "Epoch: %d\n", epoch); 

  fprintf(fo, "HI ");
  fprintf(fo, "/%d,",top(lad)->score.num_trials);
  print_score( top(lad),fo );
  fprintf(fo, " %lld %lld\n",ncomp,neval);
  if( top(lad)->agt->cd->num_cells == 1 ) {
    print_cell( top(lad)->agt->cd,0,fo );
  }
  else {
    print_code( top(lad)->agt->cd,fo );
  }
  fprintf(fo, "\n");
  fflush(fo);
}
