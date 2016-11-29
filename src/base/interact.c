/** \file interact.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "step.h"
#include "scan_print.h"
#include "cross_mutate.h"
#include "eval.h"
#include "interact.h"

#include "message.h"

#include <unistd.h>

#ifdef USE_MPI
#include <time.h>
#include <mpi.h>
#endif

#define   SWAP_FAILS       1

#define   MAX_TASK        64

Boolean  evolve_interact(Boolean incremental,int trim_int,Ladder *lad,
                         Channel *env_cfg,Channel *agt_cfg,Channel *env_out,
                         Channel *agt_out,Channel *agt_out_prev,
                         int reactive,int rank[],long seed[],int max_trials,
                         double item_cost,Library *lib,FILE *fo);

Boolean    trim_interact(int trim_comparisons,Ladder *lad,Channel *env_cfg,
                         Channel *agt_cfg,Channel *env_out,Channel *agt_out,
                         Channel *agt_out_prev,int reactive,int rank[],
                         long seed[],int max_trials,double item_cost,FILE *fo);

Boolean compare_interact(Ladder *lad,void *env_state,
                         Channel *env_cfg,Channel *agt_cfg,Channel *env_out,
                         Channel *agt_out,Channel *agt_out_prev,int reactive,
                         int rank[],long seed[],int num_trials,int max_trials,
                         double item_cost,FILE *fo);

void    eval_inter_react(Candidate *can,void *env_state,
                         Channel *env_cfg,Channel *agt_cfg,Channel *env_out,
                         Channel *agt_out,Channel *agt_out_prev,int r);

void    eval_inter_recur(Candidate *can,void *env_state,
                         Channel *env_cfg,Channel *agt_cfg,Channel *env_out,
                         Channel *agt_out,Channel *agt_out_prev,int r);

extern long eval_seed;
#ifdef USE_MPI
extern int my_rank;
#endif

/********************************************************//**
   Apply tuning algorithm to evolve a candidate
   for the specified task.
*/
void tune_interact(
                   Search_Param *par,
                   long     seed[],
                   Channel *env_cfg,
                   Channel *agt_cfg,
                   Ladder  *lad
                  )
{
  void    *env_state;
  Channel *env_out;
  Channel *agt_out;
  Channel *agt_out_prev;
  int     *rank;
  Code    *cd1;
  FILE    *fo = par->silent ? NULL : stdout;
  Level    level_tune  ={ TUNE, NON_ALIGNED};
  Level    level_interp={INTERP,NON_ALIGNED};
  Boolean  top_is_better;
  int tune_comparisons=10000;
  int prev_trials,num_trials;
  int epoch=0;
  int cmp;
  int r;

  if( par->incremental ) {
    top(lad)->score.num_trials = par->min_trials;
  }
  top(lad)->score.successful = FALSE; // force initial loop cycle

  env_state = new_env_state();

  env_out      = new_channel();
  agt_out      = new_channel();
  agt_out_prev = new_channel();

  rank = (int *)malloc((par->max_trials+1)*sizeof(int));
  check_null( rank );
  for( r=0; r <= par->max_trials; r++ ) {
    rank[r] = r;
  }

  generate_seeds( par->max_trials,seed );

  if( par->incremental ) {
    num_trials = top(lad)->score.num_trials;
  }
  else {
    num_trials = par->max_trials;
  }
  // this happens twice the first time through the loop, but that's ok
  compare_interact(lad,env_state,env_cfg,agt_cfg,env_out,agt_out,
                   agt_out_prev,par->reactive,rank,seed,
                   num_trials,par->max_trials,par->item_cost,fo);

  while(       epoch < par->max_epoch
        &&(    top(lad)->score.num_trials < par->max_trials
           || !top(lad)->score.successful )) {
    epoch++;

    num_trials = top(lad)->score.num_trials;
    compare_interact(lad,env_state,env_cfg,agt_cfg,env_out,agt_out,
                     agt_out_prev,par->reactive,rank,seed,
                     num_trials,par->max_trials,par->item_cost,fo);

    prev_trials = num_trials;

    cmp=0;

    while(       cmp < tune_comparisons
          &&(    top(lad)->score.num_trials < par->max_trials
             || !top(lad)->score.successful )) {
      cmp++;

      breed( lad,NULL,level_tune );

      top_is_better = compare_interact(lad,env_state,env_cfg,agt_cfg,env_out,
                                  agt_out,agt_out_prev,par->reactive,rank,seed,
                                  num_trials,par->max_trials,par->item_cost,fo);

      if( top(lad)->score.num_trials > num_trials ) {// incremental
        num_trials = top(lad)->score.num_trials; // increment trials
        bring_to_front( rank,num_trials );
        cmp = tune_comparisons; // break out of loop
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
          top_is_better=compare_interact(lad,env_state,env_cfg,agt_cfg,
                                  env_out,agt_out,agt_out_prev,par->reactive,rank,
                                  seed,num_trials,par->max_trials,par->item_cost,fo);
          num_trials = top(lad)->score.num_trials;
        }
      }
      else {
        cull_top( lad,fo );
      }
    }
#ifdef SWAP_FAILS
    if( par->incremental && num_trials == prev_trials
                         && num_trials  <  par->max_trials ) {
      // if no new items added, swap newest item with next item
      int next_item      = rank[num_trials+1];
      rank[num_trials+1] = rank[1];
      rank[1]            = next_item;
    }
#endif
#ifdef RESHUFFLE_FAILS
    if(check_reshuffle(lad,par->min_trials,prev_trials,
                       num_trials,par->max_trials)) {
      shuffle_seeds( par->max_trials,seed );
    }
#endif
    if( !par->silent ) {
      printf("KI ");
      if( par->incremental ) {
        printf("/%d,",top(lad)->score.num_trials);
      }
      print_score( top(lad),stdout );
      printf("\n");
      if( top(lad)->agt->cd->num_cells == 1 ) {
        print_cell( top(lad)->agt->cd,0,stdout );
      }
      else {
        print_code( top(lad)->agt->cd,stdout );
      }
      fflush( stdout );
    }
  }
  if( !par->silent ) {
    printf("FI ");
    if( par->incremental ) {
      printf("/%d,",top(lad)->score.num_trials);
    }
    print_score( top(lad),stdout );
    printf("\n");
    fflush( stdout );
  }

  free_channel( env_out );
  free_channel( agt_out );
  free_channel( agt_out_prev );
  free_env_state( env_state );
  free( rank );
}

/********************************************************//**
   Apply Ladder algorithm to evolve a candidate
   for the specified task.
*/
void search_interact(
                     Search_Param *par,
                     long     seed[],
                     Channel *env_cfg,
                     Channel *agt_cfg,
                     Ladder  *lad,
                     Library *lib
                    )
{
  Channel *env_out;
  Channel *agt_out;
  Channel *agt_out_prev;
  int     *rank;
  FILE    *fo = par->silent ? NULL : stdout;
  int prev_trials;
  int epoch=0;
  int r;

  lad->max_child_per_epoch = 10000;

  if( par->reactive ) {
    is_allowed[(int)BRB] = FALSE;
  }
  if( par->incremental ) {
    top(lad)->score.num_trials = par->min_trials;
  }
  top(lad)->score.successful = FALSE; // force initial loop cycle

  random_codebank( lad,128 );

  env_out      = new_channel();
  agt_out      = new_channel();
  agt_out_prev = new_channel();

  rank = (int *)malloc((par->max_trials+1)*sizeof(int));
  check_null( rank );
  for( r=0; r <= par->max_trials; r++ ) {
    rank[r] = r;
  }

  while(                            epoch < par->max_epoch
        &&                     lad->neval < globals.max_eval
        &&(    top(lad)->score.num_trials < par->max_trials
           || !top(lad)->score.successful )) {

    epoch++;

    prev_trials = top(lad)->score.num_trials;
    evolve_interact(par->incremental,par->trim_interim,lad,env_cfg,
               agt_cfg,env_out,agt_out,agt_out_prev,par->reactive,
               rank,seed,par->max_trials,par->item_cost,lib,fo);

    adjust_max_child( lad,prev_trials,par->max_trials );

#ifdef RESHUFFLE_FAILS
    if(check_reshuffle(lad,par->min_trials,prev_trials,
                top(lad)->score.num_trials,par->max_trials)) {
      shuffle_seeds( par->max_trials,seed );
      move_to_codebank( lad,CHAMP-1 );
      lad->max_child_per_epoch += 2000;
    }
#endif
  }
  
  if ( lad->neval < globals.max_eval ) {
    if( !par->silent ) {
      printf("FI ");
      if( par->incremental ) {
        printf("/%d,",top(lad)->score.num_trials);
      }
      print_score( top(lad),stdout );
      printf(" %lld %lld\n",lad->ncomp,lad->neval);
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

    trim_interact(par->trim_final,lad,env_cfg,NULL,env_out,agt_out,
     agt_out_prev,par->reactive,rank,seed,par->max_trials,par->item_cost,fo);

    if( !par->silent ) {
      printf("TI ");
      print_score( top(lad),stdout );
      printf(" %lld %lld\n",lad->ncomp,lad->neval);
      printf("*** TASK SOLVED AT EPOCH %d ***\n",epoch );
    }
  } else {
    printf("max_eval (%lld) exceeded\n", globals.max_eval);
    
    if( !par->silent ) {
      printf("MI ");
      print_score( top(lad),stdout );
      printf(" %lld %lld\n",lad->ncomp,lad->neval);
    }
  }

  free_channel( env_out );
  free_channel( agt_out );
  free_channel( agt_out_prev );
  free( rank );
}

/********************************************************//**
   Apply Ladder algorithm in parallel to evolve a candidate
   for the specified task.
*/
void multi_instance_interact(
                             Search_Param *inter_par,
                             Channel *env_cfg,
                             Channel *agt_cfg,
                             Code    *agt_code0,
                             Library *lib
                            )
{
  // these are all for local tasks only:
  FILE            *fo[MAX_TASK]={NULL};
 //  Channel *agt_cfg[MAX_TASK]={NULL};
  Ladder         *lad[MAX_TASK];
  int           *rank[MAX_TASK];
  Boolean incremental[MAX_TASK]={FALSE};
  Boolean      solved[MAX_TASK]={FALSE};
  Channel *env_out;
  Channel *agt_out;
  Channel *agt_out_prev;
  Code    *cd0 = NULL;
  Code    *champ;
  Level    level_copy={COPY,NON_ALIGNED};
  Boolean  keep_going=TRUE;
  // int lib_code_base=0;
  
  Boolean halt;
  
  long    seed[MAX_SEEDS];
  //int max_trials=1000;
  int prev_trials;
  int epoch=0;
  int d,r,i;
  
  if( inter_par->reactive ) {
    is_allowed[(int)BRB] = FALSE;
  }
  if( lib != NULL ) {
    lib->code_base = lib->num_code;
  }
  for( d=0; d < MAX_TASK; d++ ) {
    lad[d]->max_child_per_epoch = 10000;
    random_codebank( lad[d],128 );
  }

  env_out      = new_channel();
  agt_out      = new_channel();
  agt_out_prev = new_channel();

  // default sub-interval is the entire range of tasks
  // for 1-process case this is correct
  int interval_start = 0;
  int interval_end = inter_par->num_instances - 1;
  
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
  
  int quo = inter_par->num_instances / num_processes;
  int rem = inter_par->num_instances % num_processes;
  
  interval_end = -1;
  for (i = 0; i <= my_rank; i++) {
    interval_start = interval_end+1;
    interval_end = interval_start + quo + ((i >= rem) ? -1 : 0);
  }
#endif
  
  // number of tasks to be performed on this process:
  int num_local_instances = interval_end - interval_start + 1;
  
  for( d=0; d < num_local_instances; d++ ) {
    if (inter_par->output_dir == NULL || inter_par->output_dir[0] == '\0') {
      // no output folder, use stdout
      fo[d] = stdout;
    } else {
      i = d + interval_start; // overall index
      char filename[400];
      sprintf(filename, "%s/TASK_%d.txt",inter_par->output_dir,i);
      fo[d] = fopen(filename, "w");
    
      if (fo[d] == NULL) {
        printf("Couldn't open %s\n", filename);
      }
    
      fprintf(fo[d], "multi_instance_interact / task %d of %d.\n",
              i,inter_par->num_instances);
#ifdef USE_MPI
      fprintf(fo[d], "(on process %d of %d)\n", my_rank, num_processes);
#endif
      fflush(fo[d]);
    }
  }

#ifdef USE_MPI
  if (my_rank == 0)
#endif
  printf("multi_instance_interact: running %d tasks in total.\n",
         inter_par->num_instances);
  
#ifdef USE_MPI
  printf("Process %d assigned %d tasks.\n", my_rank, num_local_instances);
#endif
  
  fflush( stdout );
  
  // initialise each instance
  for (d = 0; d < num_local_instances; d++) {
    cd0 = cross_mutate( agt_code0,NULL,level_copy );
    // same starting code for all ladders.
    // could be more flexible but would require input file
    
    lad[d] = new_ladder(256);

    reset_ladder( lad[d],cd0 );
    if( !inter_par->incremental ) {
      incremental[d] = FALSE;
    }
    else {
      incremental[d] = TRUE;
      top(lad[d])->score.num_trials = inter_par->min_trials;
    }

    top(lad[d])->score.penalty = 1.0;

    rank[d] = (int *)malloc((inter_par->max_trials+1)*sizeof(int));
    check_null( rank[d] );
    for( r=0; r <= inter_par->max_trials; r++ ) {
      rank[d][r] = r;
    }
  }
  free_code( agt_code0 );

  // seeds should be different between processes
  generate_seeds( inter_par->max_trials,seed );
  
  halt = FALSE;

  while( epoch < inter_par->max_epoch && keep_going ) {

    epoch++;

    keep_going = FALSE;

    for( d=0; d < num_local_instances; d++ ) {
      // d is index relative to this task.
      i = d + interval_start; // overall index
      
      if( !solved[d] ) {

        if( halt && globals.terminate_early ) {
          print_termination(lad[d],lad[d]->ncomp,lad[d]->neval,
                            epoch, i, fo[d]);
          solved[d] = TRUE;
          continue;
        }
        
        keep_going = TRUE;

        clear_message(agt_out,1);
        clear_message(agt_out_prev,1);

        prev_trials = top(lad[d])->score.num_trials;
        halt |= evolve_interact(incremental[d],inter_par->trim_interim,
                                lad[d],env_cfg,NULL,env_out,agt_out,agt_out_prev,
                                inter_par->reactive,rank[d],seed,inter_par->max_trials,
                                inter_par->item_cost,lib,fo[d]);
           // agt_cfg or NULL?

        champ = cross_mutate( top(lad[d])->agt->cd,NULL,level_copy );
        insert_library( lib,lib->code_base+i,champ );
        
#ifdef USE_MPI
        if (globals.multi_process) {
          broadcast_code(champ, i, my_rank, num_processes);
        }
#endif
        //fprintf(fo[d], "Library:\n");
        //print_library(lib, fo[d]);
        /*
        if ( incremental[d] ) {
          if( top(lad[d])->score.num_trials > num_trials_[d] ) {
            num_trials_[d] = top(lad[d])->score.num_trials;
          }
        }
        */
        adjust_max_child( lad[d],prev_trials,inter_par->max_trials );

#ifdef RESHUFFLE_FAILS
        if(check_reshuffle(lad[d],inter_par->min_trials,prev_trials,
                top(lad[d])->score.num_trials,inter_par->max_trials)) {
          shuffle_seeds( inter_par->max_trials,seed );
          move_to_codebank( lad[d],CHAMP-1 );
          lad[d]->max_child_per_epoch += 2000;
        }
#endif
        // TODO: alter output when max_eval is exceeded
        if(       lad[d]->neval >= globals.max_eval
           ||(   top(lad[d])->score.num_trials >= inter_par->max_trials
              && top(lad[d])->score.successful )) {

          solved[d] = TRUE;
          fprintf(fo[d],"FI ");
          if( incremental[d] ) {
            fprintf(fo[d], "/%d,",top(lad[d])->score.num_trials);
          }
          print_score( top(lad[d]),fo[d] );
          fprintf(fo[d]," %lld %lld\n",lad[d]->ncomp,lad[d]->neval);
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
            first_to_halt = broadcast_halt(my_rank, num_processes);
            fprintf(fo[d],"Signalling other processes to terminate.\n");
#endif
            // stop local processes:
            halt = TRUE;
          }
          Boolean fin;

          if (first_to_halt) {

            fprintf(fo[d],"TRIMMING...\n");

            fflush( fo[d] );

            fin = trim_interact(inter_par->trim_final,lad[d],env_cfg, // agt_cfg or NULL
                   agt_cfg,env_out,agt_out,agt_out_prev,inter_par->reactive,
                   rank[d],seed,inter_par->max_trials,inter_par->item_cost,fo[d]);
          }
          else {
            fin = FALSE;
          }
          
          if( !fin ){ // Trimming terminated early (i.e. didn't finish)
            print_termination(lad[d],lad[d]->ncomp,lad[d]->neval, 
                              epoch,d+interval_start,fo[d]);
          } else {
            // Trimming finished, thus this is the actual solution
            fprintf(fo[d],"TI ");
            print_score( top(lad[d]),fo[d] );
            fprintf(fo[d]," %lld %lld\n",lad[d]->ncomp,lad[d]->neval);
            print_code( top(lad[d])->agt->cd,fo[d] );
            
            fprintf(fo[d],"*** TASK %d SOLVED AT EPOCH %d ***\n",
                          d+interval_start,epoch );
            print_code( top(lad[d])->agt->cd,fo[d] );

            reset_codebank(lad[d]);

            fflush( fo[d] );
          }          
        }
      }
    }
  }

#ifdef USE_MPI
#if NETWORK_TYPE == STAR_NETWORK
  signal_finished();
#endif
#endif

  free_channel( env_out );
  free_channel( agt_out );
  free_channel( agt_out_prev );

  for( d=0; d < num_local_instances; d++ ) {
    if( fo[d] != stdout ) {
      fclose( fo[d] );
    }
    free_ladder( lad[d] );
    free( rank[d] );
  }
}

/********************************************************//**
   Evolve ladder for the specified number of comparisons,
   or until the number of trials is incremented.
   Returns false if terminated early by another process
*/
Boolean evolve_interact(
                        Boolean  incremental,
                        int      trim_int,
                        Ladder  *lad,
                        Channel *env_cfg,
                        Channel *agt_cfg,
                        Channel *env_out,
                        Channel *agt_out,
                        Channel *agt_out_prev,
                        Boolean  reactive,
                        int      rank[],
                        long     seed[],
                        int      max_trials,
                        double   item_cost,
                        Library *lib,
                        FILE    *fo
                       )
{
  void   *env_state;
  Boolean top_is_better;
  int prev_trials,num_trials;
  int child=0;

  env_state = new_env_state();

  if( incremental ) {
    num_trials = top(lad)->score.num_trials;
    compare_interact(lad,env_state,env_cfg,agt_cfg,env_out,agt_out,
        agt_out_prev,reactive,rank,seed,num_trials,max_trials,item_cost,fo);
    if( top(lad)->score.num_trials > num_trials && fo != NULL ) {
      num_trials = top(lad)->score.num_trials;
      fprintf(fo, "SI ");
      fprintf(fo, "/%d,",top(lad)->score.num_trials);
      print_score( top(lad),fo );
      fprintf(fo, "\n");
    }
    num_trials = top(lad)->score.num_trials;
  }
  else {
    num_trials = max_trials;
    compare_interact(lad,env_state,env_cfg,agt_cfg,env_out,agt_out,
        agt_out_prev,reactive,rank,seed,num_trials,max_trials,item_cost,fo);
  }
  prev_trials = num_trials;
  
  Boolean halt = FALSE;

  while(  ( child < lad->max_child_per_epoch || pop(lad) != NULL )
        &&(    top(lad)->score.num_trials < max_trials
           || !top(lad)->score.successful )) {
#ifdef USE_MPI
    if (globals.multi_process) {
      // check for waiting messages,
      // update library with champs from
      // other processes
      // Also check for halt message
      //if (cmp%256 == 0) { // (faster without this)
      halt |= update_library(lib);
      if (halt) {
        break;
      }
      //}
    }
#endif

    child++;
    
    //printf("cmp = %d\n", cmp);

    procreate( lad,lib );

    prev_trials = num_trials;

    top_is_better = compare_interact(lad,env_state,env_cfg,agt_cfg,
                              env_out,agt_out,agt_out_prev,reactive,
                              rank,seed,num_trials,max_trials,item_cost,fo);

    if( incremental ) {
      if( top(lad)->score.num_trials > num_trials ) {
        num_trials = top(lad)->score.num_trials;
        bring_to_front( rank,num_trials );
        child = lad->max_child_per_epoch; // break out of loop
      }
    }
    while( top_is_better ) {
      top_replace_pop( lad,fo );
      top_is_better = better_candidate(top(lad),pop(lad),0.0,TRUE);
    }
    while( cull_top( lad,fo ))
      ;
  }
  free_env_state( env_state );

  if( trim_int > 0 &&  num_trials > prev_trials
                   &&  num_trials <  max_trials
                   && !top(lad)->score.reject
                   &&  top(lad)->agt->cd->last_codon > 1 ) {
    trim_interact(trim_int,lad,env_cfg,agt_cfg,env_out,agt_out,
            agt_out_prev,reactive,rank,seed,max_trials,item_cost,fo );
  }

  if( fo != NULL ) {
    fprintf(fo, "KI ");
    fprintf(fo, "/%d,",top(lad)->score.num_trials);
    print_score( top(lad),fo );
    fprintf(fo, " %lld %lld\n",lad->ncomp,lad->neval);
    if( top(lad)->agt->cd->num_cells == 1 ) {
      print_cell( top(lad)->agt->cd,0,fo );
    }
    else {
      print_code( top(lad)->agt->cd,fo );
    }
    fprintf(fo, "\n");
    fflush( fo );
  }
#ifdef SWAP_FAILS
  if( incremental && num_trials == prev_trials
                  && num_trials <   max_trials ) {
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
   Returns False if terminated early by another process.
*/
Boolean trim_interact(
                      int      trim_comparisons,
                      Ladder  *lad,
                      Channel *env_cfg,
                      Channel *agt_cfg,
                      Channel *env_out,
                      Channel *agt_out,
                      Channel *agt_out_prev,
                      int      reactive,
                      int      rank[],
                      long     seed[],
                      int      max_trials,
                      double   item_cost,
                      FILE    *fo
                     )
{

  // if trim_comparisons < 0, adaptive trim is used
  // i.e. t_c = -20 -> keep trimming until 20 consecutive unsuccessful trims

  void *env_state;
  Boolean top_is_better;
  Level level_trim={TRIM,NON_ALIGNED};
  int num_trials = top(lad)->score.num_trials;
  int cmp;
  Boolean halt = FALSE;

  env_state = new_env_state();
 
  if( trim_comparisons != 0 ) {
    compare_interact(lad,env_state,env_cfg,agt_cfg,env_out,agt_out,
        agt_out_prev,reactive,rank,seed,num_trials,num_trials,item_cost,fo);
  }

  int t = abs( trim_comparisons );

  if( trim_comparisons < 0 && fo != NULL ) {
    fprintf(fo, "(adaptive)\n");
  }
  cmp = 0;

  int n = 0;
  int s = 0;
  while (cmp < t && !halt) {
    
    breed( lad,NULL,level_trim );
    while( top(lad)->agt->cd->last_codon >= pop(lad)->agt->cd->last_codon ){
      cull_top( lad,fo );
      breed( lad,NULL,level_trim );
    }
    top_is_better=compare_interact(lad,env_state,env_cfg,agt_cfg,env_out,agt_out,
             agt_out_prev,reactive,rank,seed,num_trials,num_trials,item_cost,fo);
    if( top_is_better ) {
      top_replace_pop( lad,fo );
      if (trim_comparisons < 0) { // adaptive
        cmp = -1; // reset counter after successful trim
      } 
      s++;
    }
    else {
      cull_top( lad,fo );
    }
#ifdef USE_MPI
#if NETWORK_TYPE == ALL_TO_ALL
    if( globals.multi_process && globals.terminate_early ) {
      // check for halt message
      halt |= test_trim_halt(my_rank);
    }
#endif
#endif
    cmp++;
    n++;
  }

  free_env_state( env_state );

  if( fo != NULL ) {
    fprintf(fo, "%d successful out of %d trims\n",s,n);
  }
#ifdef USE_MPI
#if NETWORK_TYPE == ALL_TO_ALL
  if( halt ) { 
    return( FALSE );
  }
#endif
#endif

  return( TRUE );
}

/********************************************************//**
   Compare top(lad) against pop(lad) on the same task,
   for at least min_seeds and so long as the agent achieves
   a perfect score (no reject, no penalty, and zero cost).
   Return TRUE if the top agent is better; FALSE otherwise.
*/
Boolean compare_interact(
                         Ladder  *lad,
                         void    *env_state,
                         Channel *env_cfg,
                         Channel *agt_cfg,
                         Channel *env_out,
                         Channel *agt_out,
                         Channel *agt_out_prev,
                         int      reactive,
                         int      rank[],
                         long     seed[],
                         int      num_trials,
                         int      max_trials,
                         double   item_cost,
                         FILE     *fo
                        )
{
  Candidate *can=top(lad);
  Score   prev_score,this_score;
  Level   level = can->agt->cd->level;
  Boolean top_is_better;
  Boolean final = FALSE;
  double  target_cost = 0.0;
  int prev_trials = num_trials;
  int r=0;

  lad->ncomp++;

  reset_score( can );
  clear_message(agt_out,1);
  clear_message(agt_out_prev,1);

  while( !final ) {

    r++;
    prev_score = can->score;
    eval_seed = seed[rank[r]];
    init_evaluation_state(seed[rank[r]]);

    lad->neval++;
    if( reactive ) {
      eval_inter_react( can,env_state,env_cfg,agt_cfg,
                        env_out,agt_out,agt_out_prev,r );
    }
    else {
      eval_inter_recur( can,env_state,env_cfg,agt_cfg,
                        env_out,agt_out,agt_out_prev,r );
    }

    if( can->score.num_trials > num_trials ) {
      if( num_trials == prev_trials && fo != NULL ) {
        this_score = can->score;
        can->score = prev_score;
        fprintf(fo, "II ");
        fprintf(fo, "/%d,",can->score.num_trials);
        print_score( top(lad),fo );
        fprintf(fo, "\n");
        fflush(fo);
        can->score = this_score;
      }
      num_trials   = can->score.num_trials; // increment trials
    }

    if( item_cost >= 0.0 ) {
      target_cost = num_trials * item_cost;
      can->score.successful &= (can->score.cost < target_cost);
    }

    final=(     r >= max_trials || can->score.reject
           ||(  r >= num_trials
              &&(    level.m == TRIM
                 || !can->score.successful ))
           ||(       level.m  < BAR
              && !better_candidate(can,pop(lad),target_cost,FALSE)));
  }
  if( num_trials > prev_trials ) {
    if( num_trials < max_trials ) {
      int rank_r;
      // keep swapping failed seed with the following seed,
      // so long as agent successfully processes the following seed
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
        eval_seed = seed[rank[r]];
        init_evaluation_state(seed[rank[r]]);

        lad->neval++;
        if( reactive ) {
          eval_inter_react( can,env_state,env_cfg,agt_cfg,
                            env_out,agt_out,agt_out_prev,r );
        }
        else {
          eval_inter_recur( can,env_state,env_cfg,agt_cfg,
                            env_out,agt_out,agt_out_prev,r );
        }
        if( item_cost >= 0.0 ) {
          target_cost = r * item_cost;
          can->score.successful &= (can->score.cost < target_cost);
        }
      } while( r < max_trials && can->score.successful );

      num_trials = r;
    }
    if( num_trials > prev_trials+1 && fo != NULL ) {
      this_score = can->score;
      can->score = prev_score;
      fprintf(fo, "JI ");
      fprintf(fo, "/%d,",can->score.num_trials);
      print_score( top(lad),fo );
      fprintf(fo, "\n");
      fflush(fo);
      can->score = this_score;
    }
  }
  if( can->score.all_outputs_same &&( reactive || num_trials > 1 )) {
    //printf("all outputs same, adding penalty\n");
    can->score.penalty += 1.0;
  }

  restore_mutation_state();// otherwise same agent bred repeatedly

  top_is_better = better_candidate(can,pop(lad),target_cost,TRUE);

  return( top_is_better );
}

/********************************************************//**
*/
void test_interact(
                   Candidate *can,
                   void      *env_state,
                   Channel   *env_cfg,
                   Channel   *agt_cfg,
                   Channel   *env_out,
                   Channel   *agt_out,
                   Channel   *agt_out_prev,
                   int        reactive,
                   long       seed[],
                   int        max_trials
                  )
{
  int r=0;

  clear_message(agt_out,1);
  clear_message(agt_out_prev,1);

  while( r < max_trials && can->score.successful ) {
    r++;
    eval_seed = seed[r];
    init_evaluation_state(seed[r]);

    //lad->neval++;
    if( reactive ) {
      eval_inter_react( can,env_state,env_cfg,agt_cfg,
                        env_out,agt_out,agt_out_prev,r );
    }
    else {
      eval_inter_recur( can,env_state,env_cfg,agt_cfg,
                        env_out,agt_out,agt_out_prev,r );
    }
  }
  if( can->score.all_outputs_same &&( reactive || max_trials > 1 )) {
    //printf("all outputs same, adding penalty\n");
    can->score.penalty += 1.0;
  }

  restore_mutation_state();// otherwise same agent bred repeatedly
}

/********************************************************//**
   Allow can->agt to interact with environment (env_state) and update
   cat->score.cost,num_trials,num_steps,reject,all_outputs_same
*/
void eval_inter_react(
                      Candidate *can,
                      void    *env_state,
                      Channel *env_cfg,
                      Channel *agt_cfg,
                      Channel *env_out,
                      Channel *agt_out,
                      Channel *agt_out_prev,
                      int r
                     )
{
  Channel *agt_in=agt_cfg;
  Agent   *agt =  can->agt;
  Score   *scr = &can->score;
  Score   *env_scr;
  Boolean  env_running = TRUE;
  int max_agt_step;
  int op = NONE;

  //neval++;

  max_agt_step = AGT_STEPS_PER_MOVE;

  // Let environment run, taking inputs from env_cfg
  // until it writes its first output.
  cue_message(env_cfg,1);
  clear_message(env_out,1);
  if( !env_reset( env_state,env_cfg,env_out,r )) {
    env_scr = env_get_score( env_state );
    env_scr->reject = TRUE; // failed to produce output
    return;
  }

  // start agent running
  scr->num_trials++;
  cue_message( agt_cfg,1 );
  reset_agent( agt );

  // now agent and environment will take turns,
  // until one of them halts or exceeds max steps
  while( env_running && agt->running
                     && agt->step <= max_agt_step) {

    // reset agent at each time step
    agt_in = agt_cfg;
    cue_message( agt_cfg,1 );
    scr->num_steps += agt->step;
    reset_agent( agt );

    cue_message(env_out,env_out->om);    // latest message from env
    clear_message(agt_out,agt_out->om+1);// clear any partially written message
    op = NONE;
    // Run the agent until it outputs
    while( op != OUT && agt->running
                     && agt->step <= max_agt_step ) {
      op = step( agt,agt_in,agt_out );
        // when agt_cfg is exhausted, switch to chan
      if( op == INP && agt_in != env_out && input_exhausted(agt_in)) {
        agt_in = env_out;
        agt->bstack[agt->bp] = fetch_input( agt_in );
      }
    }
#ifdef DEBUG
    //print_state(agt,agt_in,agt_out,stdout);
#endif

    if(  op == OUT ) { // agent successfully produced output
      cue_message(agt_out,agt_out->om); // latest message from agent
      if( agt_out->om >= 2 ) {
        if(    scr->all_outputs_same
           && !same_message(agt_out,1,agt_out,agt_out->om )) {
          scr->all_outputs_same = FALSE;
        }
        if( agt_out->om > 1 ) {
          shift_message(agt_out,1,agt_out->om);
        }
      }
      env_running = env_continue( env_state,agt_out,env_out );

      if( env_out->om > 1 ) {
        shift_message(env_out,1,env_out->om);
      }      
    }
    else {
      scr->penalty += 1.0;
      agt->running = FALSE;
    }
  }

  scr->num_steps += agt->step;

  if( agt->step > max_agt_step ) {
    scr->reject = TRUE;
  }
  else {
    if( env_running ) {
      // if agent halts, environment will keep running
      // but inputs will fail and outputs will be ignored
      while( env_continue( env_state,NULL,NULL ))
        ;
    }
  }

  env_scr = env_get_score( env_state );
  if( !env_scr->reject ) { // CHECK ????!!??
    if( isfinite( env_scr->cost )) {
      if( env_scr->cost < 0.0 ) {
        env_scr->cost = 0.0;
      }
      scr->cost += env_scr->cost;
    }
    else {
      scr->penalty += 1.0;
      //scr->reject = TRUE;
    }
  }
  scr->successful &= env_scr->successful;
  scr->penalty    += env_scr->penalty;
}

/********************************************************//**
   Allow can->agt to interact with environment (env_state)
   Record cost in can->cost, resource usage in can->num_trials,
   can->num_steps. Set can->reject in the case where
   env->error_flag is TRUE, or agt exceeds max_agt_step.
   Set task->reject if env exceeds max_env_step.
   Return TRUE if all agent outputs are identical, FALSE otherwise.
*/
void eval_inter_recur(
                      Candidate *can,
                      void    *env_state,
                      Channel *env_cfg,
                      Channel *agt_cfg,
                      Channel *env_out,
                      Channel *agt_out,
                      Channel *agt_out_prev,
                      int r
                     )
{
  Channel *agt_in=agt_cfg;
  Agent   *agt =  can->agt;
  Score   *scr = &can->score;
  Score   *env_scr;
  Boolean  env_running = TRUE;
  int max_agt_step;
  int op=NONE;

  //neval++;

  max_agt_step = 10 * AGT_STEPS_PER_MOVE;

  // Let environment run, taking inputs from env_cfg
  // until it writes its first output.
  cue_message(env_cfg,1);
  clear_message(env_out,1);
  if( !env_reset( env_state,env_cfg,env_out,r )) {
    return;    // environment failed to produce output
  }

  // Let agent run, taking inputs first from agt_cfg and then
  //  from env_out, until it writes its first output.
  scr->num_trials++;
  agt_in = agt_cfg;
  cue_message( agt_cfg,1 );
  cue_message( env_out,1 );
  clear_message( agt_out,1 );
  reset_agent( agt );
  op = NONE;
  while( op != OUT && agt->running
         //        && agt->call_stack_size < MAX_CALL_STACK
                   && agt->step <= max_agt_step ) {
    op = step( agt,agt_in,agt_out );
    // when agt_cfg is exhausted, switch to env_out
    if( op == INP && agt_in != env_out && input_exhausted(agt_in)) {
      agt_in = env_out;
      agt->bstack[agt->bp] = fetch_input( agt_in );
    }
  }

  // now environment and agent will take turns,
  // until one of them halts or exceeds max steps
  while( env_running && agt->running
         //          && agt->call_stack_size < MAX_CALL_STACK
                     && agt->step <= max_agt_step ) {

    cue_message(agt_out,agt_out->om); // latest message from agent

    env_running = env_continue( env_state,agt_out,env_out );

    if( env_out->om > 1 ) {
      shift_message( env_out,1,env_out->om );
    }

    if( env_running ) {

      // Run agent until it produces output
      max_agt_step += AGT_STEPS_PER_MOVE;
      cue_message(env_out,env_out->om);    // latest message from env
      clear_message(agt_out,agt_out->om+1);// clear any partially written message
      op = NONE;
      while( op != OUT && agt->running
             //        && agt->call_stack_size < MAX_CALL_STACK
                       && agt->step <= max_agt_step ) {
        op = step( agt,env_out,agt_out );
      }
#ifdef DEBUG
      //print_state(agt,agt_in,agt_out,stdout);
#endif
    }
  }

  scr->num_steps += agt->step;

  if(   agt->step > max_agt_step ) {
        //  || agt->call_stack_size >= MAX_CALL_STACK ) {
    scr->reject = TRUE;
  }
  else {
    if( env_running ) {
      // if agent halts, environment will keep running
      // but inputs will fail and outputs will be ignored
      while( env_continue( env_state,NULL,NULL ))
        ;
    }
  }

  scr->all_outputs_same &= matching_channel(agt_out,agt_out_prev);
  if( scr->all_outputs_same && agt_out->om > agt_out_prev->om ) {
     Channel tmp_chan;
     tmp_chan     = *agt_out;
    *agt_out      = *agt_out_prev;
    *agt_out_prev =  tmp_chan;
  }

  env_scr = env_get_score( env_state );
  if( !env_scr->reject ) { // CHECK ????!!??
    if( isfinite( env_scr->cost )) {
      if( env_scr->cost < 0.0 ) {
        env_scr->cost = 0.0;
      }
      scr->cost += env_scr->cost;
    }
    else {
      scr->reject = TRUE;
    }
  }
  scr->successful &= env_scr->successful;
  scr->penalty    += env_scr->penalty;
}

/********************************************************//**
   old version of eval_interact(),
   using shared channel between agent and environment

Boolean eval_interact_shared(
                             Candidate *can,
                             void      *env_state,
                             Channel   *env_cfg,
                             Channel   *agt_cfg,
                             Channel   *chan,
                             int        reactive,
                             FILE      *fo
                            )
{
  Channel *agt_in=agt_cfg;
  Agent   *agt =  can->agt;
  Score   *scr = &can->score;
  Score   *env_scr;
  Boolean  env_running = TRUE;
  Boolean  all_outputs_same = TRUE;
  int max_agt_step;
  int op=NONE;

  neval++;

  if( reactive ) {
    max_agt_step =    AGT_STEPS_PER_MOVE;
  }
  else {
    max_agt_step = 10*AGT_STEPS_PER_MOVE;
  }

  // Let environment run, taking inputs from env_cfg
  // until it writes its first output.
  cue_message(env_cfg,1);
  if( !env_reset( env_state,env_cfg,chan )) {
    return( TRUE ); // environment failed to produce output
  }
  
  // start agent running
  scr->num_trials++;
  cue_message( agt_cfg,1 );
  reset_agent( agt );

  // now agent and environment will take turns,
  // until one of them halts or exceeds max steps
  while( env_running && agt->running
                     && agt->step <= max_agt_step ) {
  
    if( reactive ) { // reset agent each time
      agt_in = agt_cfg;
      cue_message( agt_cfg,1 );
      scr->num_steps += agt->step;
      reset_agent( agt );
    }
    else {           // agent continues running
      max_agt_step += AGT_STEPS_PER_MOVE;
    }

    cue_message(chan,chan->om); // latest message from env
    op = NONE;
    // Run the agent until it outputs
    while( op != OUT && agt->running
                     && agt->step <= max_agt_step ) {
      op = step( agt,agt_in,chan );
        // when agt_cfg is exhausted, switch to chan
      if( op == INP && agt_in != chan && input_exhausted(agt_in)) {
        agt_in = chan;
        agt->bstack[agt->bp] = fetch_input( agt_in );
      }
    }
#ifdef DEBUG
    //print_state(agt,agt_in,chan,stdout);
#endif

    if(  op == OUT ) { // agent successfully produced output
      cue_message(chan,chan->om); // latest message from agent
      // shift messages back, to save space
      if( chan->om >= 4 ) {       // [0|env|agt|env|agt|..
        if( chan->om >= 6 ) {     // [0|env|agt|env|agt|env|agt|
          shift_message(chan,3,5);//           ,-------'
        }                         // [0|env|agt|env|agt|
        if( all_outputs_same & !same_message(chan,2,chan,4)) {
          all_outputs_same = FALSE;
        }
      }
      env_running = env_continue( env_state,chan,chan );
    }
  }
  if( chan->om % 2 == 1 ) {       // [0|env|..|agt|env
    // if agent fails to output, clear message from env
    clear_message(chan,chan->om); // [0|env|..|agt|
  }
  scr->num_steps += agt->step;

  if( agt->step > max_agt_step ) {
    scr->reject = TRUE;
  }
  else {
    if( env_running ) {
      // if agent halts, environment will keep running
      // but inputs will fail and outputs will be ignored
      while( env_continue( env_state,NULL,NULL ))
        ;
    }
  }

  env_scr = env_get_score( env_state );
  if( !env_scr->reject ) { // CHECK ????!!??
    if( isfinite( env_scr->cost )) {
      if( env_scr->cost < 0.0 ) {
        env_scr->cost = 0.0;
      }
      scr->cost += env_scr->cost;
    }
    else {
      scr->reject = TRUE;
    }
  }
  if( env_scr->penalty ) {
    scr->penalty += 1.0;
  }
  return( all_outputs_same );
}
*/

