/** \file debug_interact.c

// provides step-by-step simulation of an interaction between an environment and an agent
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "step.h"
#include "scan_print.h"
#include "cross_mutate.h"
#include "eval.h"
#include "interact.h"

#define AGT_STEPS_PER_MOVE 10000


void usage( char *argv_0 );

/********************************************************//**
   Interpret the specified hercl program,
   taking input from the specified file(s).
*/
int main( int argc, char *argv[] )
{
  Agent   *agt;
  FILE     *code_file=NULL;
  Channel *agt_in, *chan;
  //char line[MAX_LINE];
  int  freq[NUM_CODONS][NUM_CODONS]={{0}};
  int  op;
  //char ch;
  char key,command='s';
  long seed=0;
  int  c=NONE,k=NONE; // breakpoint cell and codon
  int  fp=NONE; // frame pointer
  int  a=1;
  int num_trials = 1;
  
  FILE *fi;

  Boolean reactive = FALSE;

  Boolean  env_running = TRUE;
  Boolean  all_outputs_same = TRUE;
  int max_agt_step;
  
  Channel *env_cfg = NULL; // task  configuration
  Channel *agt_cfg = NULL; // agent configuration
  
  init_codons();
  
  void *env_state = new_env_state();

  agt_in = agt_cfg;
  chan = new_channel();


  if( argc < 2 ) {
    usage(argv[0]);
  }
  while( a < argc) {

    if( argv[a][0] == '-' ) {
      if ( argv[a][1] == 'b' ) { // batch mode
        command = 'g';
        a++;
      }
      else if( argv[a][1] == 'd' ) {
        if( a < argc-1 ) {
          seed = atoi(argv[a+1]);
          a += 2;
        }
        else {
          usage(argv[0]);
        }
      }
      else if( argv[a][1] == 'r' ) {
        reactive = TRUE;
        a++;
      }
      else if( argv[a][1] == 'n' ) {
        a++;
        num_trials = atoi(argv[a]);
        a++;
      } else if (argv[a][1] == 'p') {
        env_cfg = new_channel();
        a++;
        printf("cfg env from %s\n", argv[a]);
        while( a < argc && argv[a][0] != '-' ) {
          fi = fopen_check( argv[a],"r" );
          scan_inputs( env_cfg,fi );
          fclose( fi );
          a++;
        }
      } else if (argv[a][1] == 'a') {
        agt_cfg = new_channel();
        a++;
        while( a < argc && argv[a][0] != '-' ) {
          fi = fopen_check( argv[a],"r" );
          scan_inputs( agt_cfg,fi );
          fclose( fi );
          a++;
        }
      } else {
        usage(argv[0]);
      }
    } else if( code_file == NULL ) {
      code_file = fopen(argv[a],"r");
      if( code_file == NULL ) {
        printf("File not found: %s\n",argv[a]);
        exit(1);
      }
      a++;
    }
  }

  if( code_file == NULL ) {
    usage(argv[0]);
  }

  agt = new_agent();

  if( command != 'g' ) {
    printf("Scanning...\n");
  }
  agt->cd = scan_code( code_file );
  fclose( code_file );

  //k = string_from_int( s,47 );
  //k = k + string_from_code( s+k,agt->cd );
  //printf("%s\n",s);
  //k = skip_int( s );
  //agt->cd = code_from_string( s+k );

  if( command != 'g' ) {
    printf("Printing...\n");
    print_code( agt->cd, stdout );
    aggregate_freq(freq,agt->cd);
    print_freq(freq,stdout);
  }
  
  if( reactive ) {
    max_agt_step =    AGT_STEPS_PER_MOVE;
  }
  else {
    max_agt_step = 10*AGT_STEPS_PER_MOVE;
  }

  srandom( seed );


  // Begin trials.
  int trial;
  for (trial = 0; trial < num_trials; trial++) {
    printf(" *** Trial %d ***\n", trial);
    
    env_running = TRUE;
    
  
    // Let environment run, taking inputs from env_cfg
    // until it writes its first output.
    cue_message(env_cfg,1);
    if( !env_reset( env_state,env_cfg,chan,trial )) {
      //clear_message(chan,1);
      // reject env if it halts without ever producing output
      // fail("Environment failed to produce output\n");
      //printf("Environment failed to produce output\n");
      return( TRUE ); // ?????????????????
    }
    
    // start agent running
    cue_message( agt_cfg,1 );

    reset_agent( agt );
    
    while( env_running && agt->running
                       && agt->step <= max_agt_step ) {

      if( reactive ) { // reset agent each time
        agt_in = agt_cfg;
        cue_message( agt_cfg,1 );
        reset_agent( agt );
      }
      else {           // agent continues running
        max_agt_step += AGT_STEPS_PER_MOVE;
      }

      cue_message(chan,chan->om); // latest message from env
      op = NONE;

      while( op != OUT && agt->running
                       && agt->step <= max_agt_step ) {

        if( command != 'g' ) {
          Snode *sf = agt->call_stack + agt->fp;
          if(    command == 's' || command == 'h' || command == 'b'
             ||( command == 'c' && sf->c == c && sf->k == k )
             ||( command == 'n' && agt->fp <= fp )) {
            if( command == 'h' || command == 'b' ) {
              command = 's';
            }
            else {
              print_state( agt, agt_in, chan, stdout );
            }
            printf("? ");
            key = getchar();
            while( !isalpha(key) && key != '\n' && key != EOF ) {
              key = getchar();
            }
            if( key == 's' || key == 'c' || key == 'b' || key == 'h'
                           || key == 'n' || key == 'q' || key == 'g' ) {
              command = key;
            }
            if( command == 'q' ) {
              printf("Bye!\n");
              return 0;
            }
            if( command == 'b' ) {
              scanf("%d",&c);
              scanf("%d",&k);
              printf("Breakpoint at cell %d, codon %d\n",c,k);
              k += agt->cd->bar[agt->cd->cell[c]];
              while( k > 0 && is_dot_or_digit(agt->cd->codon[k-1])) {
                k--;
              }
            }
            if( command == 'n' ) {
              fp = agt->fp;
            }
            else {
              fp = NONE;
            }
            if( command == 'h' ) {
              printf(" s - Step\n");
              printf(" n - Next\n");
              printf(" b <c> <k> - set Breakpoint at cell c, codon k\n");
              printf(" c - Continue (until breakpoint)\n");
              printf(" g - Go (until program halts)\n");
            }
            while( key != '\n' ) {
              key = getchar();
            }
          }
        }

        if( command != 'b' && command != 'h' ) {
          op = step( agt, agt_in, chan );
          
          // when agt_cfg is exhausted, switch to chan
          if( op == INP && agt_in != chan && input_exhausted(agt_in)) {
            agt_in = chan;
            agt->bstack[agt->bp] = fetch_input( agt_in );
          }
          
          if( op == OUT ) {
            print_message( chan,chan->om,stdout );
          //  clear_message( chan,1 );
          }
        }
      }

//      if( command != 'g' ) {
//        print_state( agt,agt_in,chan,stdout );
//      }
      
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

    
    if( env_running ) {
      // if agent halts, environment will keep running
      // but inputs will fail and outputs will be ignored
      while( env_continue( env_state,NULL,NULL ))
        ;
    }

    Score   *env_scr;
    env_scr = env_get_score( env_state );
    if( !env_scr->reject ) {
      if( isfinite( env_scr->cost )) {
        if( env_scr->cost < 0.0 ) {
          env_scr->cost = 0.0;
        }
      }
      else {
        env_scr->reject = TRUE;
      }
    }
    
    if (env_scr->reject) {
      printf("(reject)\n");
    }
    printf("Score: P = %0.2lf, C = %0.2lf, success = %s\n", env_scr->penalty, env_scr->cost, env_scr->successful ? "True" : "False");
    
    printf("Total number of steps = %lld\n",agt->step );
  }
  
  free_env_state(env_state);

  //free_channel( agt_in );
  free_channel( chan );
  free_code( agt->cd );
  compress_agent( agt );
  free( agt );

  return 0;
}

/********************************************************//**
   Print usage information and exit
*/
void usage( char *argv_0 )
{
  fprintf(stderr,"Usage: %s [-d <seed>] [-b] file.hrc [-a agt.cfg] [-p env.cfg] [-r] [-n <num_trials>]\n",argv_0);
  exit(1);
}


