/** \file super.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "step.h"
#include "scan_print.h"
#include "point.h"
#include "cross_mutate.h"
#include "eval.h"
#include "super.h"

#include "message.h"

#include <time.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#define   MAX_FEATURES   1024
//#define   MAX_STEP     100000
#define   MAX_STEP     10000

#define   SWAP_FAILS     1

#define   MAX_TASK        64

void       super_options(char line[]);
void     super_prep_data(Boolean incremental,Search_Param *par,
                         Super *sup,int *rank,FILE *fp);
void   super_prep_ladder(Search_Param *par,Super *sup,Ladder *lad,FILE *fo);
void    super_trim_phase(Search_Param *par,Ladder *lad,Super *sup,
                         Channel *agt_cfg,Channel *out,int *rank,FILE *fo);
Boolean     evolve_super(Boolean incremental,int trim_int,Ladder *lad,
                         Super *sup,Library *crit,Channel *agt_cfg,
                         Channel *out,int rank[],Library *lib,FILE *fo);
Boolean       trim_super(int trim_comparisons,Ladder *lad,Super *sup,
                         Channel *agt_cfg,Channel *out,int rank[],FILE *fo);
Boolean    compare_super(Ladder *lad,Super *sup,Library *crit,Channel *agt_cfg,
                         Channel *out,int rank[],int set,int num_trials,FILE *fo);
void          crit_super(Candidate *can,Super *sup,Library *crit,Agent *agt,
                         Channel *agt_cfg,int *tset,int r,Channel *out,int o);
//double         crit_eval(Candidate *can,Channel *agt_cfg,
//                         Super *sup,int m,Channel *out,int o);
void          run_critic(Agent *agt,Channel *agt_cfg,
                         Super *sup,int m,Channel *out,int o);
int     choose_cost_type(Channel *target,int total_items);
void       shuffle_super(Super *sup);
void print_train_test_fitness(char *str,Super *sup,Ladder *lad,
                         Channel *agt_cfg,Channel *out,FILE *fp);
void print_train_fitness(Super *sup,Ladder *lad,FILE *fo );
void  print_test_fitness(Super *sup,Ladder *lad,Channel *agt_cfg,
                         Channel *out,FILE *fo);
int    get_test_mismatch(Super *sup,Ladder *lad,Channel *agt_cfg,Channel *out);

#ifdef USE_MPI
extern int my_rank;
#endif

/********************************************************//**
   Apply tuning algorithm to evolve a candidate
   using the specified training data.
*/
void tune_super(
                Search_Param *par,
                Super   *sup,
                Channel *agt_cfg,
                Channel *out,
                Ladder  *lad
               )
{
  int     *rank;
  Code    *cd1;
  FILE    *fo = par->silent ? NULL : stdout;
  Level    level_tune  ={ TUNE, NON_ALIGNED};
  Level    level_interp={INTERP,NON_ALIGNED};
  Boolean  top_is_better;
  int tune_comp   = 10000;
  int refine_comp =  2000;
  int prev_trials,num_trials;
  int epoch=0;
  int cmp;

  rank = (int *)malloc((sup->train_items+1)*sizeof(int));
  check_null( rank );

  if( !par->silent ) {
    printf("tune super...\n");
    fflush( stdout );
  }
  super_prep_data( par->incremental,par,sup,rank,stdout );

  //is_allowed[(int)BRB] = (    sup->cost_type == LED
  //                        || !sup->inputs_same_length );

  // NO register_features option
  sup->register_features = FALSE;

  if( par->incremental ) {
    top(lad)->score.num_trials = par->min_trials;
  }
  top(lad)->score.successful = FALSE;

  if( par->incremental ) {
    num_trials = top(lad)->score.num_trials;
  }
  else {
    num_trials = sup->train_items;
  }
  // this happens twice the first time through the loop, but that's ok
  compare_super(lad,sup,NULL,agt_cfg,out,rank,TRAIN,num_trials,fo);

  while(       epoch < par->max_epoch
        &&(    top(lad)->score.num_trials < sup->train_items
           || !top(lad)->score.successful )) {
    epoch++;

    num_trials = top(lad)->score.num_trials;
    compare_super(lad,sup,NULL,agt_cfg,out,rank,TRAIN,num_trials,fo);

    prev_trials = num_trials;

    cmp=0;

    while(       cmp < tune_comp
          &&(    top(lad)->score.num_trials < sup->train_items
             || !top(lad)->score.successful )) {
      cmp++;

      breed( lad,NULL,level_tune );

      top_is_better = compare_super(lad,sup,NULL,agt_cfg,out,
                                    rank,TRAIN,num_trials,fo);

      if( top(lad)->score.num_trials > num_trials ) {// incremental
        num_trials = top(lad)->score.num_trials; // increment trials
        bring_to_front( rank,num_trials );
        rank[0] = rank[num_trials];
        cmp = tune_comp; // break out of loop
      }
      if( top_is_better ) {
        if( globals.interp_rate == 1.0 ) {
          top_replace_pop( lad,fo );
        }
        else {
          cd1 = cross_mutate(pop(lad)->agt->cd,
                             top(lad)->agt->cd,level_interp);
          free_code( top(lad)->agt->cd );
          top(lad)->agt->cd = cd1;
          top_replace_pop( lad,fo );
          num_trials = top(lad)->score.num_trials;
          //reset_score(top(lad));
          top_is_better = compare_super(lad,sup,NULL,agt_cfg,out,
                                        rank,TRAIN,num_trials,fo);
          num_trials = top(lad)->score.num_trials;
        }
      }
      else {
        cull_top( lad,fo );
      }
    }
#ifdef SWAP_FAILS
    if( par->incremental && num_trials == prev_trials
                         && num_trials <  sup->train_items ) {
      // if no new items added, swap newest item with next item
      int next_item      = rank[num_trials+1];
      rank[num_trials+1] = rank[1];
      rank[1]            = next_item;
    }
#endif
#ifdef RESHUFFLE_FAILS
    if( check_reshuffle(lad,par->min_trials,prev_trials,
                        num_trials,sup->train_items)) {
      printf("XX %d",sup->train_set[rank[1]]);
      printf(" %d\n",sup->train_set[rank[num_trials+1]]);
      shuffle_super( sup );
    }
#endif
    if( !par->silent ) {
      print_train_test_fitness("KI",sup,lad,agt_cfg,out,stdout);
      if( top(lad)->agt->cd->num_cells == 1 ) {
        print_cell( top(lad)->agt->cd,0,stdout );
      }
      else {
        print_code( top(lad)->agt->cd,stdout );
      }
      printf("\n");
      fflush( stdout );
    }
  }


  if( !par->silent ) {
    if(perfect_score(top(lad),sup->item_cost*sup->train_items)){
      print_train_test_fitness("FI",sup,lad,agt_cfg,out,stdout);
    }
    else {
      print_train_test_fitness("GI",sup,lad,agt_cfg,out,stdout);
    }
    print_code( top(lad)->agt->cd,stdout );
    // tune for a fixed number of comparisons
    printf("REFINING...\n");
    fflush( stdout );
  }
  num_trials = sup->train_items;
  sup->item_cost = 0.0;
  top_is_better=compare_super(lad,sup,NULL,agt_cfg,out,rank,TRAIN,num_trials,fo);
  cmp = 0;

  while(       cmp < refine_comp
        &&( !perfect_score( top(lad),0.0 ))) {

    cmp++;

    breed( lad,NULL,level_tune );

    top_is_better = compare_super(lad,sup,NULL,agt_cfg,out,
                                  rank,TRAIN,num_trials,fo);
    if( top_is_better ) {
      if( globals.interp_rate == 1.0 ) {
        top_replace_pop( lad,fo );
      }
      else {
        cd1 = cross_mutate( pop(lad)->agt->cd,
                            top(lad)->agt->cd,level_interp );
        free_code( top(lad)->agt->cd );
        top(lad)->agt->cd = cd1;
        top_replace_pop( lad,fo );
        top_is_better = compare_super(lad,sup,NULL,agt_cfg,out,
                                      rank,TRAIN,num_trials,fo);
      }
    }
    else {
      cull_top( lad,fo );
    }
  }

  if( !par->silent ) {
    print_train_test_fitness("TI",sup,lad,agt_cfg,out,stdout);
  }

  //if( sup->test_items > 0 ) { // return cost on TEST set
  //  reset_score(top(lad));
  //  test_super( top(lad),sup,agt_cfg,out,TEST );
  //}

  //free_channel( out );
  free( rank );
}

/********************************************************//**
   Apply Ladder algorithm to evolve a candidate
   using the specified training data.
*/
void search_super(
                  Search_Param *par,
                  Super   *sup,
                  Channel *agt_cfg,
                  Channel *out,
                  Ladder  *lad,
                  Library *lib
                 )
{
  Ladder  *lad_crit=NULL;
  Library *crit=NULL;
  int     *rank;
  FILE    *fo = par->silent ? NULL : stdout;
  Boolean still_searching,reshuffled;
  int prev_trials;
  int epoch;

  rank = (int *)malloc((sup->train_items+1)*sizeof(int));
  check_null( rank );

  if( !par->silent ) {
    printf("search super...\n");
  }
  super_prep_data( par->incremental,par,sup,rank,stdout );

  is_allowed[(int)BRB] = (    sup->cost_type == LED
                          || !sup->inputs_same_length );

  super_prep_ladder( par,sup,lad,stdout );

  //if( par->adversarial ) {
  //if( 1 == 0 ) {
  {
    Code *crit_code = empty_code(10,par->num_cells,64,64);
    int m;
    lad_crit = new_ladder(256);
    reset_ladder( lad_crit,crit_code );
    super_prep_ladder( par,sup,lad_crit,stdout );
    m = crit_code->level.m+1;
    lad_crit->bank[m] = new_library(256);
    crit = lad_crit->bank[m];
  }

  epoch = 0;

  still_searching = TRUE;

  while( still_searching ) {

    reshuffled = FALSE;

    while( still_searching & !reshuffled ) {

      epoch++;

      prev_trials = top(lad)->score.num_trials;
      evolve_super(par->incremental,par->trim_interim,
                   lad,sup,crit,agt_cfg,out,rank,lib,fo);

      still_searching =                 epoch < par->max_epoch
            &&                     lad->neval < globals.max_eval
            &&(    top(lad)->score.num_trials < sup->train_items
               || !top(lad)->score.successful );
    /*{
      int m;
      for( m=0; m < MAX_LEVEL; m++ ) {
        if( lad->bank[m] != NULL && lad->bank[m]->code != NULL ) {
          printf("m=%d,\n",m);
          print_library( lad->bank[m],stdout );
        }
      }
    }*/
      adjust_max_child( lad,prev_trials,sup->train_items );

#ifdef RESHUFFLE_FAILS
      reshuffled = check_reshuffle(lad,par->min_trials,prev_trials,
                       top(lad)->score.num_trials,sup->train_items);
      if( reshuffled ) {
        printf("XX %d",sup->train_set[rank[1]]);
        printf(" %d\n",sup->train_set[rank[prev_trials+1]]);
        shuffle_super( sup );
        move_to_codebank( lad,CHAMP-1 );
        lad->max_child_per_epoch += 2000;
      }
#endif
    }
  }

  super_trim_phase( par,lad,sup,agt_cfg,out,rank,fo );

  printf("*** TASK COMPLETED AT EPOCH %d ***\n",epoch );

  free( rank );
}

/********************************************************//**
   Apply Ladder algorithm to evolve a candidate
   using the specified training data.
*/
void super_trim_phase(
                      Search_Param *par,
                      Ladder  *lad,
                      Super   *sup,
                      Channel *agt_cfg,
                      Channel *out,
                      int     *rank,
                      FILE    *fo
                     )
{
  if ( lad->neval < globals.max_eval ) {
    if( !par->silent ) {
      print_train_test_fitness("FI",sup,lad,agt_cfg,out,stdout);
      print_code( top(lad)->agt->cd,stdout );
#ifdef PRINT_HISTOGRAM
      if( globals.verbose ) {
        print_hist( lad,stdout );
      }
#endif
      print_stats( lad,stdout );
      printf("TRIMMING...\n");
      fflush( stdout );
    }
    sup->item_cost = 0.0;
    trim_super( par->trim_final,lad,sup,agt_cfg,out,rank,fo );
    //trim_super( 100000,lad,sup,agt_cfg,out,rank,stdout,FALSE );
    //trim_super( 20000,lad,sup,agt_cfg,out,rank,FALSE );

    if( !par->silent ) {
      print_train_test_fitness("TI",sup,lad,agt_cfg,out,stdout);
    }
  }
  else {
    printf("max_eval (%lld) exceeded\n",globals.max_eval);
    if( !par->silent ) {
      print_train_test_fitness("MI",sup,lad,agt_cfg,out,stdout);
    }
  }
  //if( sup->test_items > 0 ) { // return cost on TEST items
  //  reset_score(top(lad));
  //  test_super( top(lad),sup,agt_cfg,out,TEST );
  //}

  //free_channel( out );
}

/********************************************************//**
   Apply Ladder algorithm to evolve a candidate
   using the specified training data.
*/
void multi_super(
                 Search_Param *super_par,
                 Code    *agt_code0,
                 Library *lib
                )
{
  // these are all for local tasks only:
  FILE            *fo[MAX_TASK]={NULL};
  Channel    *agt_cfg[MAX_TASK]={NULL};
  Code      *agt_code[MAX_TASK]={NULL};
  Super          *sup[MAX_TASK];
  Ladder         *lad[MAX_TASK];
  int           *rank[MAX_TASK];
  Boolean incremental[MAX_TASK]={FALSE};
  Boolean      solved[MAX_TASK]={FALSE};

  Code    *champ;
  Channel *out;
  Level    level_copy={COPY,NON_ALIGNED};
  Boolean  fin;
  Boolean  keep_going=TRUE;
  Boolean  halt = FALSE;
  char     line[1024];
  char     word[1024];
  int max_train_items=0;
  int num_tasks;
  int prev_trials;
  int epoch=0;
  int d,j,k,i;

  Search_Param param = {
     NULL,  // file_in
     NULL,  // filename
     NULL,  // output_dir
    FALSE,  // is_tune
    FALSE,  // reactive
    FALSE,  // is_multi
    FALSE,  // silent
    FALSE,  // incremental
        2,  // min_trials
     1000,  // max_trials
        0,  // num_instances
     NONE,  // task
   100000,  // max_epoch
    10000,  // trim_interim
   100000,  // trim_final
        1,  // num_cells
     NONE,  // cost_type
     -1.0   // item_cost
  };

  if( lib != NULL ) {
    lib->code_base = lib->num_code;
  }
  for( d=0; d < MAX_TASK; d++ ) {
    fo[d] = stdout;
  }

  // count lines to determine total number of tasks.
  num_tasks = 0;
  while(     fgets( line,1024,super_par->file_in ) > 0
        && ( j = next_word( 0,line,word )) > 0 ) {
    num_tasks++;
  }
  rewind( super_par->file_in );

  out = new_channel();

  // default sub-interval is the entire range of tasks
  // for 1-process case this is correct
  int interval_start = 0;
  int interval_end = num_tasks - 1;
  
#ifdef USE_MPI
  // Initialize the MPI environment
  // MPI_Init(NULL, NULL);
  // Find out rank, size
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int num_processes;
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

  if (num_processes > 1) {
    globals.multi_process = TRUE;
  }

// 0 is master, so don't count it as a worker process
#if NETWORK_TYPE == STAR_NETWORK
  if (globals.multi_process) {
    num_processes--;
    my_rank--;
  }
#endif
  
  // determine sub-interval for this process
  // (which tasks this process will handle)
  
  // [this could be done centrally, but having
  // each processor compute this independently
  // saves the need for extra message passing]
  
  int quo = num_tasks / num_processes;
  int rem = num_tasks % num_processes;
  
  interval_end = -1;
  for( i = 0; i <= my_rank; i++ ) {
    interval_start = interval_end+1;
    interval_end = interval_start + quo + ((i >= rem) ? -1 : 0);
  }
#endif
  
  // number of tasks to be performed on this process:
  int num_local_tasks = interval_end - interval_start + 1;
  
  for( d=0; d < num_local_tasks; d++ ) {
    if( super_par->output_dir == NULL || super_par->output_dir[0] == '\0') {
      // no output folder, use stdout
      fo[d] = stdout;
    }
    else {
      i = d + interval_start; // overall index
      char filename[400];
      sprintf(filename, "%s/TASK_%d.txt",super_par->output_dir,i);
      fo[d] = fopen(filename, "w");
    
      if (fo[d] == NULL) {
        printf("Couldn't open %s\n", filename);
      }
    
      fprintf(fo[d], "multi_instance_super / task %d of %d.\n", i, num_tasks);
#ifdef USE_MPI
      fprintf(fo[d], "(on process %d of %d)\n", my_rank, num_processes);
#endif
      fflush(fo[d]);
    }
  }

#ifdef USE_MPI
  if( my_rank == 0 )
#endif
  printf("multi_super: running %d tasks in total.\n", num_tasks);
  
#ifdef USE_MPI
  printf("Process %d assigned %d tasks.\n", my_rank, num_local_tasks);
#endif

  fflush( stdout );

  d = -1;
  i = -1;
  // process each line
  while(     fgets( line,1024,super_par->file_in ) > 0
        && ( j = next_word( 0,line,word )) > 0 ) {

    i++;

    if( interval_start <= i && i <= interval_end ) {
      // this is always the case when using only 1 process

      d = i - interval_start;
      // d is the index of the task relative to this process,
      // i is the overall index of the task

      param.num_cells = -1.0; // sentinel
      param.item_cost = super_par->item_cost;
      param.cost_type = super_par->cost_type;
      incremental[d] = FALSE;

      fprintf(fo[d], "This task: %s\n",line);
      param.file_in = fopen_check( word,(char *)"r" );

      sup[d] = scan_super( param.file_in );
      fclose( param.file_in );

      //fo[d] = fopen_check( word,"w" ); // ?????!!!!!!!!

      // scan parameters
      while(( j = next_word( j,line,word )) > 0 ) {

        if( word[0] != '-' ) {
          super_options( line );
        }
        if( word[1] == 'n' ) {   // incremental learning
          incremental[d] = TRUE;
          if( ( k = next_word( j,line,word )) > 0
             && isdigit( word[0] )) {
            param.min_trials = atoi(word);
            j = k;
          }
          else {
            param.min_trials = 2;
          }
        }
        else {
          switch( word[1] ) {

          case 'i':
            if(( j = next_word( j,line,word )) < 0 ) {
              super_options( line );
            }
            param.file_in = fopen_check( word,(char *)"r" );
            agt_code[d] = scan_code( param.file_in );
            fclose( param.file_in );
            print_code( agt_code[d],fo[d] );
            break;

          case 'o':
            if(( j = next_word( j,line,word )) < 0 ) {
              super_options( line );
            }
            fo[d] = fopen_check( word,(char *)"w" );
            break;

          case 'a':
            agt_cfg[d] = new_channel();
            while(( j = next_word( j,line,word )) > 0 && word[0] != '-' ) {
              param.file_in = fopen_check( word,(char *)"r" );
              scan_inputs( agt_cfg[d], param.file_in );
              fclose( param.file_in );
            }
            if( word[0] == '-' ) { // back up to previous word
              while( j > 0 && !isspace( line[j-1] )) {
                j--;
              }
            }
            break;

          case 'c':
            if(( j = next_word( j,line,word )) < 0 ) {
              super_options( line );
            }
            param.num_cells = atoi( word );
            break;

          case 'f':
            if(( j = next_word( j,line,word )) < 0 ) {
              super_options( line );
            }
            if( strcmp(word,"gled") == 0 ) {
              param.cost_type = LED;
            }
            else if( strcmp(word,"lin") == 0 ) {
              param.cost_type = LIN;
            }
            else if( strcmp(word,"sqr") == 0 ) {
              param.cost_type = SQR;
            }
            else {
              super_options( line );
            }
            break;

          case 't':
            if(( j = next_word( j,line,word )) < 0 ) {
              super_options( line );
            }
            sscanf( word,"%lf",&(param.item_cost));
            break;

          default:
            super_options( line );
            break;
          }
        }
      }

      if( agt_code[d] == NULL ) {
        if( param.num_cells > 0 ) {
          agt_code[d] = empty_code(10,param.num_cells,64,64);
        }
        else {
          agt_code[d] = cross_mutate(agt_code0,NULL,level_copy);
        }
      }

      lad[d] = new_ladder(256);
      reset_ladder( lad[d],agt_code[d] );

      rank[d] = (int *)malloc((sup[d]->train_items+1)*sizeof(int));
      check_null( rank[d] );

      super_prep_data(incremental[d],&param,sup[d],rank[d],fo[d]);

      if( param.item_cost < 0.0 && super_par->item_cost >= 0.0 ) {
        sup[d]->item_cost = super_par->item_cost; // default
      }
      if( sup[d]->train_items > max_train_items ) {
        max_train_items = sup[d]->train_items;
      }

      super_prep_ladder( &param,sup[d],lad[d],fo[d] );

      //fgets( line,1024,fi );
    }
  }
  fclose( super_par->file_in );
  free_code( agt_code0 );

  // num_tasks = d+1;

  while( epoch < super_par->max_epoch && keep_going ) {

    epoch++;

    keep_going = FALSE;
    //printf("epoch = %d\n",epoch);
    //fflush( stdout );

    for( d=0; d < num_local_tasks; d++ ) {
      // d is index relative to this task.
      i = d + interval_start; // overall index

      if( !solved[d] ) {

        if( halt && globals.terminate_early ) {
          print_termination(lad[d],lad[d]->ncomp,
                            lad[d]->neval,epoch,i,fo[d]);
          solved[d] = TRUE;
          continue;
        }

        keep_going = TRUE;

        is_allowed[(int)BRB] = (    sup[d]->cost_type == LED
                                || !sup[d]->inputs_same_length );
        clear_message(out,1);

        prev_trials = top(lad[d])->score.num_trials;
        halt |= evolve_super(incremental[d],super_par->trim_interim,
                      lad[d],sup[d],NULL,NULL,out,rank[d],lib,fo[d]);

        adjust_max_child( lad[d],prev_trials,sup[d]->train_items );

#ifdef RESHUFFLE_FAILS
        if( check_reshuffle(lad[d],param.min_trials,prev_trials,
                      top(lad[d])->score.num_trials,sup[d]->train_items)){
          printf("XX %d",sup[d]->train_set[rank[d][1]]);
          printf(" %d\n",sup[d]->train_set[rank[d][prev_trials+1]]);
          shuffle_super( sup[d] );
          move_to_codebank( lad[d],CHAMP-1 );
          lad[d]->max_child_per_epoch += 2000;
        }
#endif
        champ = cross_mutate( top(lad[d])->agt->cd,NULL,level_copy );
#ifdef SHARING
        insert_library( lib,lib->code_base+i,champ );
#endif        
#ifdef USE_MPI
        if (globals.multi_process) {
          broadcast_code(champ, i, my_rank, num_processes);
        }
#endif
        //if( top(lad[d])->score.num_trials > prev_trials_[d] ) {
        //  prev_trials_[d] = top(lad[d])->score.num_trials;
        //}

        // TODO: alter output when max_eval exceeded
        if(        lad[d]->neval >= globals.max_eval
           ||(   top(lad[d])->score.num_trials >= sup[d]->train_items
              && top(lad[d])->score.successful )) {

          solved[d] = TRUE;
          print_train_test_fitness("FI",sup[d],lad[d],agt_cfg[d],out,fo[d]);
          print_code( top(lad[d])->agt->cd,fo[d] );
          fflush( fo[d] );
          
#ifdef PRINT_HISTOGRAM
          if( globals.verbose ) {
            print_hist( lad[d],fo[d] );
          }
#endif
          print_stats( lad[d],fo[d] );

          Boolean first_to_halt = TRUE;

          if( globals.terminate_early ) {
            // broadcast halt message.
#ifdef USE_MPI
            fprintf(fo[d],"Signalling other processes to terminate.\n");
            first_to_halt = broadcast_halt(my_rank, num_processes);
#endif

            // stop local processes:
            halt = TRUE;
          }
          if (first_to_halt) {

            fprintf(fo[d],"TRIMMING...\n");

            fflush( fo[d] );

            sup[d]->item_cost = 0.0;

            fin = trim_super(super_par->trim_final,lad[d],sup[d],
                             agt_cfg[d],out,rank[d],fo[d]);
          }
          else {
            fin = FALSE;
          }
          
          if( !fin ){ // Trimming terminated early
            print_termination(lad[d],lad[d]->ncomp,lad[d]->neval, 
                              epoch, d+interval_start, fo[d]);
          }
          else {
            print_train_test_fitness("TI",sup[d],lad[d],agt_cfg[d],out,fo[d]);
            print_code( top(lad[d])->agt->cd,fo[d] );
            
            fprintf(fo[d],"*** TASK %d SOLVED AT EPOCH %d ***\n",
                          d+interval_start,epoch );
            print_code( top(lad[d])->agt->cd,fo[d] );

            fflush( fo[d] );
          }

          reset_codebank(lad[d]);
        }
      }
    }
  }

#ifdef USE_MPI
#if NETWORK_TYPE == STAR_NETWORK
  signal_finished();
#endif
#endif

  for( d=0; d < num_local_tasks; d++ ) {

    if( fo[d] != stdout ) {
      fclose( fo[d] );
    }

    free_super( sup[d] );
    free_ladder( lad[d] );
    free( rank[d] );
  }
  free_channel( out );
#ifdef USE_MPI
  //MPI_Finalize();
#endif
}

/********************************************************//**
   Prepare parameters for supervised training.
*/
void super_prep_data(
                     Boolean incremental,
                     Search_Param *par,
                     Super *sup,
                     int   *rank,
                     FILE  *fp
                    )
{
  int r;

  if( !par->silent ) {
    fprintf(fp,"ITEMS %d %d\n",sup->train_items,sup->test_items);
    fflush( fp );
  }
  if( incremental ) {
    shuffle_super( sup );
  }
  if( par->item_cost >= 0.0 ) {
    sup->item_cost = par->item_cost;
  }
  if( par->cost_type != NONE ) {
    sup->cost_type = par->cost_type;
  }
  if( globals.verbose ) {
    print_super( sup,fp );
  }
  for( r=0; r <= sup->train_items; r++ ) {
    rank[r] = r;
  }
}

/********************************************************//**
   Prepare ladder for supervised training
*/
void super_prep_ladder(
                       Search_Param *par,
                       Super  *sup,
                       Ladder *lad,
                       FILE   *fp
                      )
{
  lad->max_child_per_epoch = 10000;
  lad->num_fails = 0;
  lad->ncomp = 0;
  lad->neval = 0;

  if( sup->register_features ) {
    if( sup->input->length[1] > 10 ) {
      top(lad)->agt->cd->num_reg = sup->input->length[1];
      if( !par->silent ) {
        fprintf(fp,"num_reg = %d\n",top(lad)->agt->cd->num_reg);
        fflush(fp);
      }
    }
  }
  random_codebank( lad,128 );
  if( par->incremental ) {
    top(lad)->score.num_trials = par->min_trials;
  }
  top(lad)->score.successful = FALSE;

  if( !par->silent ) {
    print_code( top(lad)->agt->cd,fp );
  }
}

/********************************************************//**
   Print error message and fail
*/
void super_options( char line[] )
{
  printf("Unrecognized option:\n");
  printf("%s",line);
  printf("Available options:\n");
  fprintf(stderr,"     -o output.out     (output file)\n");
  fprintf(stderr,"     -n [num_trials]   (incremental)\n");
  fprintf(stderr,"     -p param1.in ..   (env or seq params)\n");
  fprintf(stderr,"     -a agent.cfg ..        (agent config)\n");
  fprintf(stderr,"     -i <file.hrc>     (initial code)\n");
  fprintf(stderr,"     -c <num_cells>\n");
  fprintf(stderr,"     -f <cost_type>\n");
  fprintf(stderr,"     -t <target_cost>\n");
  exit( 1 );
}

/********************************************************//**
   Apply Ladder algorithm to evolve a candidate
   using the specified training data.
*/
Boolean evolve_super(
                     Boolean  incremental,
                     int      trim_int,
                     Ladder  *lad,
                     Super   *sup,
                     Library *crit,
                     Channel *agt_cfg,
                     Channel *out,
                     int      rank[],
                     Library *lib,
                     FILE    *fo
                    )
{
  double  target_cost;
  Boolean top_is_better;
  int prev_trials,num_trials;
  int child=0;

  if( incremental ) {
    num_trials = top(lad)->score.num_trials;
    compare_super( lad,sup,crit,agt_cfg,out,rank,TRAIN,num_trials,fo );
    if( top(lad)->score.num_trials > num_trials && fo != NULL ) {
      num_trials = top(lad)->score.num_trials;
      fprintf( fo,"SI ");
      print_train_fitness( sup,lad,fo );
      fprintf( fo,"\n");
    }
  }
  else {
    num_trials = sup->train_items;
    compare_super(lad,sup,crit,agt_cfg,out,rank,TRAIN,num_trials,fo);
  }
  target_cost = num_trials * sup->item_cost;
  prev_trials = num_trials;

  Boolean halt = FALSE;
  
  while(  ( child < lad->max_child_per_epoch || pop(lad) != NULL )
        &&(    top(lad)->score.num_trials < sup->train_items
           || !top(lad)->score.successful )) {
    child++;
    
#ifdef USE_MPI
    if( globals.multi_process ) {
      // check for waiting messages,
      // update library with champs from
      // other processes
      // Also check for halt message
      //   if (cmp%256 == 0) {
      halt |= update_library(lib);
      if (halt) {
        break;
      }
      //  }
    }
#endif

#ifdef DEBUG
    //di += sprintf(&ds[di],"cmp=%d\n",cmp);
#endif
    
    procreate( lad,lib );

    //reset_score(top(lad)); // already in compare_super()
    //prev_trials = num_trials; ??????!!!!!!!!!!
    top_is_better = compare_super(lad,sup,crit,agt_cfg,out,
                                  rank,TRAIN,num_trials,fo);

    if( top(lad)->score.num_trials > num_trials ) {// incremental
      num_trials = top(lad)->score.num_trials; // increment trials
      bring_to_front( rank,num_trials );
      target_cost = num_trials * sup->item_cost;
      child = lad->max_child_per_epoch; // break out of loop
    }

    while( top_is_better ) {
      top_replace_pop( lad,fo );
      top_is_better = better_candidate(top(lad),pop(lad),target_cost,TRUE);
    }

    while( cull_top( lad,fo ))
      ;
  }

  if(    num_trials > prev_trials
     &&  num_trials < sup->train_items
     && !top(lad)->score.reject
     &&  top(lad)->agt->cd->last_codon > 1 ) {
    trim_super( trim_int,lad,sup,agt_cfg,out,rank,fo );
  }

  if( fo != NULL ) {
    print_train_test_fitness("KI",sup,lad,agt_cfg,out,fo);
    if( top(lad)->agt->cd->num_cells == 1 ) {
      print_cell( top(lad)->agt->cd,0,fo );
    }
    else {
      print_code( top(lad)->agt->cd,fo );
    }
    fprintf( fo,"\n");
    fflush( fo );
  }
#ifdef SWAP_FAILS
  if( incremental && num_trials == prev_trials
                  && num_trials <  sup->train_items ) {
    // if no new items added, swap newest item with next item
    int next_item      = rank[num_trials+1];
    rank[num_trials+1] = rank[1];
    rank[1]            = next_item;
  }
#endif
  return halt;
}

/********************************************************//**
   Try to simplify (shorten) the code if we can.
   return FALSE if terminated early by another process.
*/
Boolean trim_super(
                   int      trim_comparisons,
                   Ladder  *lad,
                   Super   *sup,
                   Channel *agt_cfg,
                   Channel *out,
                   int      rank[],
                   FILE    *fo
                  )
{
  Code   *champ;
  Level   level_copy={COPY,NON_ALIGNED};
  Level   level_trim={TRIM,NON_ALIGNED};
  Boolean top_is_better;
  Boolean halt = FALSE;
  double  champ_cost;
  int num_trials = top(lad)->score.num_trials;
  int mismatch,new_mismatch;
  int cmp=0;
  int s=0, n=0;

  // if trim_comparisons < 0, adaptive trim is used
  // i.e. t_c = -20 -> keep trimming until 20 consecutive unsuccessful trims

  //print_code(top(lad)->agt->cd,stdout);

  champ = cross_mutate( top(lad)->agt->cd,NULL,level_copy );

  compare_super(lad,sup,NULL,agt_cfg,out,rank,TRAIN,num_trials,fo);
  champ_cost = top(lad)->score.cost;

  mismatch = get_test_mismatch( sup,lad,agt_cfg,out );

  int t = abs(trim_comparisons);

  if( trim_comparisons < 0 && fo != NULL ) {
    fprintf(fo, "(adaptive)\n");
  }

  while( !halt && cmp   < t
               && lad->neval < globals.max_eval + 100000000 - num_trials){

    breed( lad,NULL,level_trim );
    top_is_better = compare_super(lad,sup,NULL,agt_cfg,out,
                                  rank,TRAIN,num_trials,fo);
    if( top_is_better ) {
      top_replace_pop( lad,fo );
      // replace champ only if cost strictly lower,
      // or cost equal and code is shorter
      if(      top(lad)->score.cost <  champ_cost
         ||(   top(lad)->score.cost == champ_cost
            && top(lad)->agt->cd->last_codon < champ->last_codon )) {
        free_code( champ );
        champ = cross_mutate( top(lad)->agt->cd,NULL,level_copy );
        champ_cost = top(lad)->score.cost;
        new_mismatch = get_test_mismatch( sup,lad,agt_cfg,out );
        if(( globals.verbose || new_mismatch != mismatch )&& fo != NULL ) {
          print_train_test_fitness("LI",sup,lad,agt_cfg,out,fo);
          if( top(lad)->agt->cd->num_cells == 1 ) {
            print_cell( top(lad)->agt->cd,0,fo );
          }
          else {
            print_code( top(lad)->agt->cd,fo );
          }
          fprintf(fo,"\n");
          fflush( fo );
        }
        mismatch = new_mismatch;
      }
      s++;
      if (trim_comparisons < 0) { // adaptive
        cmp = -1; // reset counter after successful trim
      }
    }
    else {
      cull_top( lad,fo );
    }

#ifdef USE_MPI
#if NETWORK_TYPE == ALL_TO_ALL
    // TODO test only every few iterations to save time
    if( globals.multi_process && globals.terminate_early ) {
      // check for halt message
      halt |= test_trim_halt(my_rank);
    }
#endif
#endif
    cmp++;
    n++;
  }
  if( fo != NULL ) {
    fprintf(fo, "%d successful out of %d trims\n",s,n);
  }
#if NETWORK_TYPE == ALL_TO_ALL
  if( halt ) {
    return( FALSE ); // signals early termination
  }
#endif

  // install the champ at the top of the ladder
  free_code( top(lad)->agt->cd );
  top(lad)->agt->cd = champ;

  return( TRUE );
}

/********************************************************//**
   Print current train and test fitness, ncomp and neval
*/
void print_train_test_fitness(
                              char    *str,
                              Super   *sup,
                              Ladder  *lad,
                              Channel *agt_cfg,
                              Channel *out,
                              FILE    *fp
                             )
{
  fprintf( fp,"%s ",str );
  print_train_fitness( sup,lad,fp );
  print_test_fitness( sup,lad,agt_cfg,out,fp );
  printf(" %lld %lld\n",lad->ncomp,lad->neval);
}

/********************************************************//**
   Test top(lad) against pop(lad) on items from the specified set
  (TRAIN or TEST) in super data, store the outputs to *out,
   update can->score and return result of comparison.
*/
Boolean compare_super(
                      Ladder  *lad,
                      Super   *sup,
                      Library *crit,
                      Channel *agt_cfg,
                      Channel *out,
                      int      rank[],
                      int      set,
                      int      num_trials,
                      FILE    *fo
                     )
{
  Candidate *can=top(lad);
  Agent  *agt_crit=NULL;
  Score   prev_score,this_score;
  int    *tset; // train or test set
  Boolean top_is_better;
  Boolean final=FALSE;
  double  target_cost;
  int     level_m = can->agt->cd->level.m;
  int     prev_trials = num_trials;
  int     max_trials;
  int r=0;

  lad->ncomp++;

  if( set == TRAIN ) {
    tset       = sup->train_set;
    max_trials = sup->train_items;
  }
  else {
    tset       = sup->test_set;
    max_trials = sup->test_items;
  }
  target_cost = num_trials * sup->item_cost;

  reset_score( can );

  while( !final ) {
    r++;
    prev_score = can->score;
    eval_super( lad,can,sup,agt_cfg,out,tset,rank[r],r );
    crit_super( can,sup,crit,agt_crit,agt_cfg,tset,rank[r],out,r );
    if( can->score.num_trials > num_trials ) {
      if( num_trials == prev_trials && fo != NULL ) {
        this_score = can->score;
        can->score = prev_score;
        fprintf( fo,"II ");
        print_train_fitness( sup,lad,fo );
        fprintf( fo,"\n");
        can->score = this_score;
      }
      num_trials  = can->score.num_trials; // increment trials
      target_cost = num_trials * sup->item_cost;
    }
    final=(     r >= max_trials || can->score.reject
           ||(  r >= num_trials
              &&(    level_m == TRIM
                 || !perfect_score(can,target_cost)))
           ||(       level_m  < BAR
              && !better_candidate(can,pop(lad),target_cost,FALSE)));
  }
  if( num_trials > prev_trials ) {
    if( num_trials < max_trials ) {
      int rank_r;
      // keep swapping failed item with the following item,
      // so long as agent successfully processes the following item
      can->score = prev_score;
      r = num_trials-1;
      do {
        r++;
        if( r < max_trials ) {
          rank_r    = rank[r];
          rank[r]   = rank[r+1];
          rank[r+1] = rank_r;
        }
        prev_score = can->score;
        eval_super( lad,can,sup,agt_cfg,out,tset,rank[r],r );
        crit_super( can,sup,crit,agt_crit,agt_cfg,tset,rank[r],out,r );
        target_cost = r * sup->item_cost;
      } while(   r < max_trials
              && perfect_score(can,target_cost));
      num_trials = r;
    }
    if( num_trials > prev_trials+1 && fo != NULL ) {
      this_score = can->score;
      can->score = prev_score;
      fprintf( fo,"JI ");
      print_train_fitness( sup,lad,fo );
      fprintf( fo,"\n");
      can->score = this_score;
    }
  }
  can->score.successful = perfect_score(can,target_cost);
  can->score.all_outputs_same = all_items_same( out );
  if( r > 1 && can->score.all_outputs_same ) {
    can->score.penalty += 1.0;
  }
  top_is_better = better_candidate(can,pop(lad),target_cost,TRUE);

  fflush(fo);
  return( top_is_better );
}

/********************************************************//**
   Test candidate agent on all items in specified set
  (TRAIN or TEST) in super data, store the outputs to *out
   and update can->score.
*/
void test_super(
                Candidate *can,
                Super   *sup,
                Channel *agt_cfg,
                Channel *out,
                int set
               )
{
  int *tset;
  int  num_items;
  int r;

  if( set == TRAIN ) {
    num_items = sup->train_items;
    tset      = sup->train_set;
  }
  else {
    num_items = sup->test_items;
    tset      = sup->test_set;
  }
  for( r=1; r <= num_items; r++ ) {
    eval_super( NULL,can,sup,agt_cfg,out,tset,r,r );
  }
  can->score.all_outputs_same = all_items_same(out);
  if( num_items > 1 && can->score.all_outputs_same ) {
    can->score.penalty += 1.0;
  }
}

/********************************************************//**
   Test candidate agent on item r in the specified set
  (TRAIN or TEST) in super data and store the output to *out(o).
   Update can->score.num_trials,num_steps, cost and penalty.
*/
void eval_super(
                Ladder  *lad,
                Candidate *can,
                Super   *sup,
                Channel *agt_cfg,
                Channel *out,
                int     *tset,
                int r,
                int o
               )
{
  Agent   *agt    = can->agt;
  Channel *agt_in = agt_cfg;
  Score   *scr    = &can->score;
  int op;
  int m,n;
  if( scr->reject) {
    return;
  }
  if( lad != NULL ) {
    lad->neval++;
  }
  clear_message(out,o);
  m = tset[r];
  reset_agent( agt );

  if( sup->register_features ) {
    int i=0; // load features into registers
    while( i < agt->cd->num_reg && i < sup->input->length[m] ) {
      agt->reg[i] = sup->input->val[sup->input->index[m]+i];
      i++;
    }
    sup->input->max_im = 0;  // no inputs allowed
  }
  else {
    n = 1 + (m-1)*sup->inputs_per_item;
    cue_message( sup->input, n );
    sup->input->max_im = n + sup->inputs_per_item - 1;
  }
  cue_message( agt_cfg,1 );

  while( agt->running && agt->step <  MAX_STEP
         //           && agt->call_stack_size < MAX_CALL_STACK
                      && out->om   <= o ) {
    op = step( agt, agt_in, out );
    //print_state( agt,agt_in,out,stdout );
    if(   op == INP && agt_in != sup->input
       && input_exhausted( agt_in )) {
      if( sup->register_features ) {
        agt_in = NULL;
      }
      else {
        agt_in = sup->input;
        agt->bstack[agt->bp] = fetch_input( agt_in );
      }
      //print_state( agt,agt_in,out,stdout );
    }
    //print_state(agt,sup->input,out,stdout);
    //keych = getchar();
  }
  //out->wi = out->index[out->om];
  //print_message( out,n,stdout );

  scr->num_steps += agt->step;
  if( agt->running ) { // time limit or call stack exceeded
    scr->reject  = TRUE;
  }
  cue_message( sup->target, m );
  cue_message( out, o );
  update_score( can,out,sup->target,sup->cost_type );
}

/********************************************************//**
   Test each critic on item r in the specified set
  (TRAIN or TEST) in super data and evaluate the output *out(o).
   Update can->score.cost
*/
void crit_super(
                Candidate *can,
                Super   *sup,
                Library *crit,
                Agent   *agt,
                Channel *agt_cfg,
                int     *tset,
                int r,
                Channel *out,
                int o
               )
{
  Score *scr = &can->score;
  double val;
  int k,l,m;

  m = tset[r];
  if( scr->reject || crit == NULL || crit->num_code == 0
                  || same_message( sup->target,m,out,o )) {
    return;
  }
  for( l=0; l < crit->num_code; l++  ) {
    agt->cd = crit->code[l];
    run_critic( agt,agt_cfg,sup,m,out,o );
    // scr->num_steps += agt->step;
    if( agt->running ) { // time limit or call stack exceeded
      //scr->reject  = TRUE;
      // remove this code from the library of critics
      free_code( crit->code[l] );
      crit->num_code--;
      for( k = l; k < crit->num_code; k++ ) {
        crit->code[k] = crit->code[k+1];
      }
      l--;
    }
    else if( agt->stack_items > 0 ) {
      //task->score.cost = top_of_stack( env );
      val = top_of_stack( agt );
      if( isfinite( val )) {
        if( val < 0.0 ) {
          scr->cost += 1.0;
        }
        else if ( val < 1.0 ) {
          scr->cost += ( 1.0 - val );
        }
      }
    }
    // add val to cost
  }
}

/********************************************************//**
   Run the critic (specified in can->agt->cd) on input m of super
   and evaluate out(o) as a putative answer.
   Return the evaluation, which is the top item on the stack.
*/
void run_critic(
                Agent   *agt,
                Channel *agt_cfg,
                Super   *sup,
                int      m,
                Channel *out,
                int      o
               )
{
  Channel *agt_in = agt_cfg;
  int op;
  int n;

  if( agt->cd == NULL ) {
    return;
  }
  cue_message( out, o );
  reset_agent( agt );

  if( sup->register_features ) {
    int i=0; // load features into registers
    while( i < agt->cd->num_reg && i < sup->input->length[m] ) {
      agt->reg[i] = sup->input->val[sup->input->index[m]+i];
      i++;
    }
    sup->input->max_im = 0;  // no inputs allowed
  }
  else {
    n = 1 + (m-1)*sup->inputs_per_item;
    cue_message( sup->input, n );
    sup->input->max_im = n + sup->inputs_per_item - 1;
  }
  cue_message( agt_cfg,1 );

  while( agt->running && agt->step < MAX_STEP ) {
         //           && agt->call_stack_size < MAX_CALL_STACK
         //           && op != OUT ) {
    op = step( agt, agt_in, NULL );
    //print_state( agt,agt_in,out,stdout );
    if( op == INP && agt_in != out && input_exhausted( agt_in )) {
      if( agt_in == agt_cfg && !sup->register_features ) {
        agt_in = sup->input;
      }
      else {
        agt_in = out;
      }
      agt->bstack[agt->bp] = fetch_input( agt_in );
    }
    // penalize for output ??
    //print_state(agt,sup->input,out,stdout);
    //keych = getchar();
  }
  //out->wi = out->index[out->om];
  //print_message( out,n,stdout );



  //cue_message( sup->target, m );
  //cue_message( out, o );
}

/********************************************************//**
   Choose the cost type (LED, KLD, KLD2, LIN or SQR)
   based on whether the targets differ in length,
   are between 0 and 1, etc.
*/
int choose_cost_type( Channel *target, int total_items )
{
  int targets_differ_in_length=FALSE;
  int all_targets_between_zero_and_one=TRUE;
  int all_targets_between_one_and_minus_one=TRUE;
  int target_length=0;
  int cost_type;
  int i,m;
  if( total_items > 0 ) {
    target_length = target->length[1];
    for( m=1; m <= total_items; m++ ) {
      targets_differ_in_length |=
        ( target->length[m] != target_length );
      for( i = target->index[m]; i < target->index[m+1]; i++ ){
        all_targets_between_zero_and_one &=
          ( target->val[i] >=  0.0 && target->val[i] <= 1.0 );
        all_targets_between_one_and_minus_one &=
          ( target->val[i] >= -1.0 && target->val[i] <= 1.0 );
      }
    }
  }
  if( targets_differ_in_length || total_items == 1 ) {
    cost_type = LED;
  }
  else if( all_targets_between_zero_and_one ) {
    cost_type = KLD;
  }
  else if( all_targets_between_one_and_minus_one ) {
    cost_type = SGN;
  }
  else if( target_length == 1 ) {
    cost_type = LIN;
  }
  else {
    cost_type = SQR;
  }
  return( cost_type );
}

/********************************************************//**
   Return pointer to a new, empty super.
*/
Super * new_super()
{
  Super *sup;
  sup = (Super *)malloc(sizeof(Super));
  check_null(sup);
  sup->input  = NULL;
  sup->target = NULL;
  sup->item_cost = 0.0;
  sup->cost_type = LIN;
  sup->inputs_per_item = 1;
  sup->inputs_same_length = TRUE;
  sup->register_features = FALSE;
  sup->train_set = NULL;
  sup->test_set  = NULL;
  sup->train_items = 0;
  sup->test_items  = 0;
  sup->total_items = 0;
  return( sup );
}

/********************************************************//**
   Scan several inputs from fi, until EOF, and store them
   successively into buf[], along with their input_lengths.
   Return the number of inputs scanned.
*/
Super * scan_super( FILE *fi )
{
  char line[MAX_LINE];
  Super *sup;
  int is_target = FALSE;
  int is_test = FALSE;
  int ch;
  int i,r,s;

  sup = new_super();
  sup->input  = new_channel();
  sup->target = new_channel();

  ch = getc( fi );
  ungetc( ch, fi );
  if( isdigit( ch )) {
    fscanf( fi,"%d",&sup->inputs_per_item );
    while( ch != '\n' ) {// scan to end of line
      ch = getc( fi );
    }    
  }
  while( fgets(line,MAX_LINE,fi) > 0 ) {
    ch = getc( fi );
    if( ch != '\n' ) {
      ungetc( ch, fi );
      line[strlen(line)-1] = '\0';
    }
    if( is_target ) {
      scan_next_input( sup->target,line );
      is_target = FALSE;
    }
    else {
      scan_next_input( sup->input,line );
      if( sup->input->length[sup->input->om] < 0 ) {
        clear_message( sup->input,sup->input->om );
        is_test = TRUE;
      }
      else { 
        sup->inputs_same_length &=
          (   sup->input->length[sup->input->om]
           == sup->input->length[1] );
        for( i=1; i < sup->inputs_per_item; i++ ) {
          fgets(line,MAX_LINE,fi);
          ch = getc( fi );
          if( ch != '\n' ) {
            ungetc( ch, fi );
            line[strlen(line)-1] = '\0';
          }
          scan_next_input( sup->input,line );
          sup->inputs_same_length &=
            (   sup->input->length[sup->input->om]
             == sup->input->length[i+1] );
        }
        is_target = TRUE;
        if( is_test ) {
          sup->test_items++;
        }
        else {
          sup->train_items++;
        }
      }
    }
  }
  sup->total_items = sup->train_items + sup->test_items;
  sup->train_set=(int *)malloc((sup->total_items+1)*sizeof(int));
  check_null(sup->train_set);
  sup->test_set =(int *)malloc((sup->total_items+1)*sizeof(int));
  check_null(sup->test_set);
  for( r=1; r <= sup->train_items; r++ ) {
    sup->train_set[r] = r;
  }
  for( s=1; s <= sup->test_items; s++ ) {
    sup->test_set[s] = s + sup->train_items;
  }
  sup->cost_type = choose_cost_type(sup->target,sup->total_items);
  switch( sup->cost_type ) {
   case KLD: case KLD2:
     sup->item_cost = 0.1;
     break;
   case LED:
     sup->item_cost = 0.0;
     break;
   default:
     sup->item_cost = 0.2;
  }
  
  if(   sup->inputs_per_item == 1 && sup->inputs_same_length
     && sup->cost_type != LED     && sup->input->length[1] > 6 ){
    sup->register_features = TRUE;
  }
  return( sup );
}

/********************************************************//**
   Free the space used by super.
*/
void free_super( Super *sup )
{
  if( sup != NULL ) {
    free_channel( sup->input );
    free_channel( sup->target );
    if( sup->train_set != NULL ) {
      free( sup->train_set );
    }
    if( sup->test_set != NULL ) {
      free( sup->test_set );
    }
    free( sup );
  }
}

/********************************************************//**
   Return TRUE if all items are identical; FALSE otherwise.
*/
int all_items_same( Channel *chan )
{
  int all_same;
  int m,i;
  if( chan->om < 2 ) {
    return FALSE;
  }
  all_same = TRUE;
  for( m=2;(m <= chan->om && all_same); m++ ){
    if( chan->length[m] != chan->length[1] ) {
      all_same = FALSE;
    }
    else {
      for( i=0; i < chan->length[1]; i++ ) {
        if(   chan->val[chan->index[m]+i]
           != chan->val[chan->index[1]+i] ) {
          all_same = FALSE;
        }
      }
    }
  }
  return( all_same );
}

/********************************************************//**
   Shuffle input items into a random order
*/
void shuffle_super( Super *sup )
{
  int *rank;
  int r,n;

  rank = (int *)malloc((sup->train_items+1)*sizeof(int));
  check_null( rank );
  memset( rank,0,(sup->train_items+1)*sizeof(int));

  for( n=1; n <= sup->train_items; n++ ) {
    r = 1 + random()% sup->train_items;
    while( rank[r] != 0 ) {
      r = 1 + random()% sup->train_items;
    }
    rank[r] = sup->train_set[n];
  }
  for( r=1; r <= sup->train_items; r++ ) {
    sup->train_set[r] = rank[r];
  }
  free( rank );
}

/********************************************************//**
   Print all training inputs and target outputs.
*/
void print_super( Super *sup, FILE *fo )
{
  int i,m,n,r,s;
  fprintf( fo,"Supervised Data:\n"); // ????????????
  fprintf( fo,"%d\n",sup->inputs_per_item );
  for( r=1; r <= sup->train_items; r++ ) {
    m = sup->train_set[r];
    for( i=0; i < sup->inputs_per_item; i++ ) {
      n = i + 1 + (m-1)*sup->inputs_per_item;
      print_message( sup->input, n, fo );
    }
    print_message( sup->target, m, fo );
  }
  fprintf(fo,"\\\n");
  for( s=1; s <= sup->test_items; s++ ) {
    m = sup->test_set[s];
    for( i=0; i < sup->inputs_per_item; i++ ) {
      n = i + 1 + (m-1)*sup->inputs_per_item;
      print_message( sup->input, n, fo );
    }
    print_message( sup->target, m, fo );
  }
}

/********************************************************//**
*/
void print_train_fitness( Super *sup, Ladder *lad, FILE *fo )
{
  fprintf(fo,"/%d,",top(lad)->score.num_trials);
  print_score( top(lad),fo );
  if( sup->cost_type == KLD || sup->cost_type == SGN ) {
    fprintf(fo," %3d  ",top(lad)->score.mismatch);
  }
}

/********************************************************//**
*/
void print_test_fitness(
                        Super   *sup,
                        Ladder  *lad,
                        Channel *agt_cfg,
                        Channel *out,
                        FILE    *fo
                       )
{
  Score top_score;
  if( sup->test_items > 0 ) {
    top_score = top(lad)->score;
    reset_score(top(lad));
    test_super( top(lad),sup,agt_cfg,out,TEST );
    print_score(top(lad),fo );
    if( sup->cost_type == KLD || sup->cost_type == SGN ) {
      fprintf(fo," %3d",top(lad)->score.mismatch);
    }
    top(lad)->score = top_score;
  }
}

/********************************************************//**
*/
int get_test_mismatch(
                      Super   *sup,
                      Ladder  *lad,
                      Channel *agt_cfg,
                      Channel *out
                     )
{
  Score top_score;
  int mismatch=0;
  if( sup->test_items > 0 ) {
    top_score = top(lad)->score;
    reset_score(top(lad));
    test_super( top(lad),sup,agt_cfg,out,TEST );
    mismatch = top(lad)->score.mismatch;
    top(lad)->score = top_score;
  }
  return( mismatch );
}

/********************************************************//**
   Return TRUE if items m,n are identical; FALSE otherwise.

int identical_items( Channel *chan, int m, int n )
{
  int i;
  if( chan->om < m || chan->om < n ) {
    return FALSE;
  }
  if( chan->length[m] != chan->length[n] ) {
    return FALSE;
  }
  for( i=0; i < chan->length[m]; i++ ) {
    if(   chan->val[chan->index[m]+i]
       != chan->val[chan->index[n]+i] ) {
      return FALSE;
    }
  }
  return TRUE;
}
*/

/********************************************************//**
  Shuffle the input features into a random order

void shuffle_features( Super *sup )
{
  Channel *in = sup->input;
  Floating val[MAX_FEATURES+1];
  int rank[MAX_FEATURES+1];
  int input_length;
  int m,n,r;
  if( sup->inputs_per_item > 1 ) {
    fprintf(stderr,"Cannot shuffle: %d inputs per item.\n",
            sup->inputs_per_item );
  }
  else if( sup->total_items > 0 ) {
    input_length = in->length[1];
    for( m=2; m <= sup->total_items; m++ ) {
      if( in->length[m] != input_length ) {
        fprintf(stderr,"Cannot shuffle: inputs differ in length.\n");
        return;
      }
    }
    for( r=0; r < input_length; r++ ) {
      rank[r] = -1;
    }
    for ( n=0; n < input_length; n++ ) {
      r = random() % input_length;
      while ( rank[r] >= 0) {
        r = random() % input_length;
      }
      rank[r] = n;
    }
    for( m=1; m <= sup->total_items; m++ ) {
      for( r=0; r < input_length; r++ ) {
        val[r] = in->val[in->index[m]+rank[r]];
      }
      for( r=0; r < input_length; r++ ) {
        in->val[in->index[m]+r] = val[r];
      }
    }
  }
}
*/

