/** \file hercsearch.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "step.h"
#include "scan_print.h"
#include "eval.h"
#include "super.h"
#include "interact.h"
#include "inter_hercl.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

void parse_args(int argc,char *argv[],Search_Param *param,long *p_seed,
     Channel **p_env_cfg,Channel **p_agt_cfg,Code **p_agt_code,Library *lib);
void        usage(char *argv_0);
void blank_spaces(FILE *fo,int num_spaces);


/********************************************************//**
*/
int main( int argc, char *argv[] )
{
  Library   *lib      = NULL;  // code library
  Channel   *agt_cfg  = NULL;  // agent configuration
  Channel   *env_cfg  = NULL;  // environment configuration
  Code      *env_code = NULL;  // environment code
  Code      *agt_code = NULL;  // initial code
  Super     *sup;
  Ladder    *lad;
  long      *seed;
  long  main_seed = 0;
  int p;

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
  1000000,  // max_epoch
    10000,  // trim_interim
   100000,  // trim_final
        1,  // num_cells
     NONE,  // cost_type
     -1.0   // item_cost
  };

#ifdef USE_MPI
  int rc,my_rank;
  rc = MPI_Init(&argc,&argv);
  if( rc != MPI_SUCCESS ) {
    printf("Error starting MPI program. Terminating.\n");
    MPI_Abort( MPI_COMM_WORLD, rc );
  }
#endif

  for( p=0; p < argc; p++ ) {
    printf(" %s",argv[p]);
  }
  printf("\n");

  lib  = new_library( MAX_LIB );
  init_codons();

  is_allowed[(int)PAY] = FALSE;
  is_allowed[(int)RAN] = FALSE;

#ifdef INTEGERS_ONLY
  is_allowed[(int)ASN] = FALSE;
  is_allowed[(int)EXP] = FALSE;
  is_allowed[(int)INT] = FALSE;
  is_allowed[(int)LOG] = FALSE;
  is_allowed[(int)RCP] = FALSE;
  is_allowed[(int)SQT] = FALSE;
  is_allowed[(int)TNH] = FALSE;
  is_allowed[(int)PLR] = FALSE;
  is_allowed[(int)TRG] = FALSE;
#endif

#ifdef NO_BRANCHES
  is_allowed[(int)BRB] = FALSE;
  is_allowed[(int)BRF] = FALSE;
  is_allowed[(int)EQL] = FALSE;
  is_allowed[(int)GRT] = FALSE;
  is_allowed[(int)FLS] = FALSE;
#endif
  parse_args(argc,argv,&param,&main_seed,&env_cfg,&agt_cfg,&agt_code,lib);

  if( param.filename != NULL ) {
    param.file_in = fopen_check( param.filename,(char *)"r" );
  }
  if( main_seed == 0 ) {
    main_seed = time( NULL );
  }
#ifdef USE_MPI  
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  main_seed += my_rank;
#endif
  printf( "seed=%ld\n",main_seed );
  init_mutation_state( main_seed );

  if( agt_code == NULL ) {
    if( param.is_tune ) {
      fail((char *)"initial code must be specified for tuning\n");
    }
    agt_code = empty_code(10,param.num_cells,64,64);
  }

  if( param.task == SUPER ) {
    if( param.is_multi ) {
      param.max_epoch = 100000;
      multi_super( &param,agt_code,lib );
    }
    else {
      Channel *out = new_channel();
      sup = scan_super( param.file_in );
      fclose( param.file_in );
      if( param.is_tune ) {
        param.max_epoch = 1000000;
        lad = new_ladder(0);           // don't need codebank
        reset_ladder( lad,agt_code );
        tune_super( &param,sup,agt_cfg,out,lad );
      }
      else {
        param.max_epoch = 100000;
        lad = new_ladder(256);
        reset_ladder( lad,agt_code );
        random_codebank( lad,128 );
        search_super( &param,sup,agt_cfg,out,lad,lib );
      }
      print_code( top(lad)->agt->cd,stdout );
      fflush( stdout );
      free_channel( out );
      free_super( sup );
      free_ladder( lad );
    }
  }
  /*
  else if( task == Predict ) { // NOT FULLY IMPLEMENTED
    if( is_multi ) {
      //
    }
    else {
      search_predict( 10000,fi,cd,inc,item_cost,lib );
    }
  }
  */
  else if( param.task == INTERACT ) {
    if( param.is_multi ) {
      // ??????????
    }
    else {
      FILE *fp = fopen_check( param.filename,(char *)"r" );
      env_code = scan_code( fp );
      fclose( fp );
      set_env_code( env_code );
      printf("Environment code:\n");
      print_code( env_code, stdout );
      fflush( stdout );
      seed = (long *)malloc((param.max_trials+1)*sizeof(long));
      check_null( seed );
      generate_seeds( param.max_trials,seed );
      if( param.is_tune ) {
        param.max_epoch = 1000000;
        lad = new_ladder(0);           // don't need codebank
        reset_ladder( lad,agt_code );
        tune_interact( &param,seed,env_cfg,agt_cfg,lad );
      }
      else {
        param.max_epoch = 100000;
        lad = new_ladder(256);
        reset_ladder( lad,agt_code );
        search_interact( &param,seed,env_cfg,agt_cfg,lad,lib );
      }
      print_code( top(lad)->agt->cd,stdout );
      fflush( stdout );
      free( seed );
      free_ladder( lad );
    }
  }
  else {
    usage( argv[0] );
  }

  // free_param(&param);

  clear_library( lib );
  free( lib );
  free_channel( agt_cfg );
  free_channel( env_cfg );
#ifdef USE_MPI
  MPI_Finalize();
#endif
  
  printf("Done.\n");

  return 0;
}

/********************************************************//**
   Parse command-line arguments
*/
void parse_args(
                int   argc,
                char *argv[],
                Search_Param *par,
                long     *p_seed,
                Channel **p_env_cfg,
                Channel **p_agt_cfg,
                Code    **p_agt_code,
                Library  *lib
               )
{
  FILE *fp = NULL;
  char libname[256];
  char    line[256];
  int p=1;
  while( p < argc ) {

    if( argv[p][0] != '-' || argv[p][1] == 'h' ) {
      usage( argv[0] );
    }
    if( argv[p][1] == 'n' ) {  // incremental learning
      par->incremental = TRUE;
      p++;
      if( p < argc && isdigit( argv[p][0] )) {
        par->min_trials = atoi(argv[p]);
        p++;
      }
      else {
        par->min_trials = 2;
      }
    }
    else if( argv[p][1] == 'k' ) { // max_trials
      p++;
      if( p < argc && isdigit(argv[p][0])) {
        par->max_trials = atoi( argv[p] );
        p++;
      }
      else {
        usage( argv[0] );
      }
    }
    else if( argv[p][1] == 'g' ) { // max_eval;
      p++;
      if( p < argc && isdigit(argv[p][0])) {
        globals.max_eval = atoll( argv[p] );
        p++;
      }
      else {
        usage( argv[0] );
      }
    }
    else if( argv[p][1] == 'v' ) {
      globals.verbose = TRUE;
      p++;
    }
    else if( argv[p][1] == 'u' ) {
      par->is_tune = TRUE;
      p++;
      if( p < argc && argv[p][0] != '-' ) {
        sscanf( argv[p],"%lf",&globals.tune_magnitude );
        p++;
      }
      if( p < argc && argv[p][0] != '-' ) {
        sscanf( argv[p],"%lf",&globals.interp_rate );
        p++;
      }
    }
    else if( argv[p][1] == 'r' ) {
      if( argv[p][2] != 'e' ) {
        usage( argv[0] );
      }
      if( argv[p][3] == 'a' ) { // react
        par->reactive = TRUE;
      }
      else {                    // recur
        par->reactive = FALSE;
      }
      p++;
    }
    else if( argv[p][1] == 'j' ) {
      globals.terminate_early = TRUE;
      p++;
    }
    else if( argv[p][1] == 'x' ) {
      globals.no_share = TRUE;
      p++;
    }
    else if( p+2 > argc ) {
      usage( argv[0] );
    }
    else {
      switch( argv[p][1] ) {

      case 's': case 'S':
        if( par->task != NONE ) {
          usage( argv[0] );
        }
        par->task = SUPER;
        par->filename =  argv[p+1]; 
        par->is_multi =( argv[p][1] == 'S' );
        p += 2;
        break;

      case 'q': case 'Q':
        if( par->task != NONE ) {
          usage( argv[0] );
        }
        par->task = PREDICT;
        par->filename =  argv[p+1]; 
        par->is_multi =( argv[p][1] == 'Q' );
        p += 2;
        break;

      case 'e': case 'E':
        if( par->task != NONE ) {
          usage( argv[0] );
        }
        par->task = INTERACT;
        par->filename =  argv[p+1]; 
        par->is_multi =( argv[p][1] == 'E' );
        p += 2;
        break;

      case 'l':
        scan_code_to_library( lib,argv[p+1] );
        p += 2;
        break;

      case 'L':
        fp = fopen_check( argv[p+1],(char *)"r" );
        while( fgets( line,256,fp ) > 0 ) {
          if( sscanf( line,"%256s",libname ) > 0 ) {
            scan_code_to_library( lib,libname );
          }
        }
        fclose(fp);
        p += 2;
        break;

      case 'd':
        *p_seed = atoi( argv[p+1] );
        p += 2;
        break;

      case 'z':
        par->trim_interim = atoi( argv[p+1] );
        p += 2;
        break;

      case 'p':
        *p_env_cfg = new_channel();
        p++;
        while( p < argc && argv[p][0] != '-' ) {
          fp = fopen_check( argv[p],(char *)"r" );
          scan_inputs( *p_env_cfg,fp );
          fclose( fp );
          p++;
        }
        break;

      case 'a':
        *p_agt_cfg = new_channel();
        p++;
        while( p < argc && argv[p][0] != '-' ) {
          fp = fopen_check( argv[p],(char *)"r" );
          scan_inputs( *p_agt_cfg,fp );
          fclose( fp );
          p++;
        }
        break;

      case 'i':
        fp = fopen_check( argv[p+1],(char *)"r" );
        *p_agt_code = scan_code( fp );
        fclose( fp );
        p += 2;
        break;

      case 'c':
        par->num_cells = atoi( argv[p+1] );
        if( par->num_cells <= 0 ) {
          usage( argv[0] );
        }
        p += 2;
        break;

      case 'f':
        if( strcmp(argv[p+1],"gled") == 0 ) {
          par->cost_type = LED;
        }
        else if( strcmp(argv[p+1],"lin") == 0 ) {
          par->cost_type = LIN;
        }
        else if( strcmp(argv[p+1],"sqr") == 0 ) {
          par->cost_type = SQR;
        }
        else {
          usage( argv[0] );
        }
        p += 2;
        break;
        
      case 'o':
        par->output_dir = argv[p+1];
       // check_null(par->output_dir);
       // strcpy(par->output_dir, argv[p+1]);
        p += 2;
        break;

      case 't':
        sscanf( argv[p+1],"%lf",&(par->item_cost));
        p += 2;
        break;

      default:
        usage( argv[0] );
        break;
      }
    }
  }
}

/********************************************************//**
   Print usage information and exit
*/
void usage( char *argv_0 )
{
  fprintf(stderr,"Usage: %s ( -s data.in\n",argv_0);
  blank_spaces(stderr,strlen(argv_0));
  fprintf(stderr,"        ( -S data.list\n");
  blank_spaces(stderr,strlen(argv_0));
  fprintf(stderr,"        ( -q seq.in\n");
  blank_spaces(stderr,strlen(argv_0));
  fprintf(stderr,"        ( -Q seq.list\n");
  blank_spaces(stderr,strlen(argv_0));
  fprintf(stderr,"        ( -e env.hrc\n");
  blank_spaces(stderr,strlen(argv_0));
  fprintf(stderr,"        ( -E env.list\n");
  fprintf(stderr,"       -l lib.hrc  ..\n");
  fprintf(stderr,"       -L lib.list ..\n");
  fprintf(stderr,"       -h                     (help)\n");
  fprintf(stderr,"       -u [alpha] [magnitude] (tune)\n");
  fprintf(stderr,"       -d <seed>\n");
  fprintf(stderr,"       -n [num_trials]        (incremental)\n");
  fprintf(stderr,"       -k <max_trials>        (max trials)\n");
  fprintf(stderr,"       -react                 (reactive)\n");
  fprintf(stderr,"       -recur                 (recurrent)\n");
  fprintf(stderr,"       -x                     (no sharing)\n");
  fprintf(stderr,"       -v                     (verbose)\n");
  fprintf(stderr,"       -a agt.cfg ..          (agent config)\n");
  fprintf(stderr,"       -p env.cfg ..          (environment config)\n");
  fprintf(stderr,"       -i <file.hrc>          (initial code)\n");
  fprintf(stderr,"       -c <num_cells>\n");
  fprintf(stderr,"       -f <cost_type>\n");
  fprintf(stderr,"       -g <max_eval>\n");
  fprintf(stderr,"       -t <target_cost>\n");
  fprintf(stderr,"       -z <interim_trim_comparisons>\n");
  exit(1);
}

/********************************************************//**
   Print the specified number of spaces
*/
void blank_spaces( FILE *fo, int num_spaces )
{
  int i;
  for( i=0; i < num_spaces; i++ ) {
    putc(' ',fo);
  }
}

