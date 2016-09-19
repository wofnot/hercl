/** \file hercl.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "step.h"
#include "scan_print.h"


void usage( char *argv_0 );

/********************************************************//**
   Interpret the specified hercl program,
   taking input from the specified file(s).
*/
int main( int argc, char *argv[] )
{
  Agent   *agt;
  FILE     *code_file=NULL;
  FILE    *input_file=NULL;
  Channel *in, *out, *feat=NULL;
  char line[MAX_LINE];
  int  freq[NUM_CODONS][NUM_CODONS]={{0}};
  int  op;
  char ch,key,command='s';
  long seed=0;
  int  c=NONE,k=NONE; // breakpoint cell and codon
  int  fp=NONE; // frame pointer
  int  a=1;

  in  = new_channel();
  out = new_channel();

  init_codons();

  if( argc < 3 ) {
    usage(argv[0]);
  }
  while( a < argc &&( code_file == NULL || input_file == NULL )) {

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
      else {
        usage(argv[0]);
      }
    }
    else if( argv[a][0] == '+' ) {
      FILE *reg_file;
      a++;
      if( feat != NULL || a >= argc ) {
        usage(argv[0]);
      }
      reg_file = fopen(argv[a],"r");
      if( reg_file == NULL ) {
        printf("File not found: %s\n",argv[a]);
        exit(1);
      }
      a++;
      fgets( line,MAX_LINE,reg_file );
      ch = getc( reg_file );
      if( ch != '\n' ) {
        ungetc( ch, reg_file );
        line[strlen(line)-1] = '\0';
      }
      fclose(reg_file);
      feat = new_channel();
      clear_message(feat,1);
      scan_next_input(feat,line);
    }
    else if( code_file == NULL ) {
      code_file = fopen(argv[a],"r");
      if( code_file == NULL ) {
        printf("File not found: %s\n",argv[a]);
        exit(1);
      }
      a++;
    }
    else {
      input_file = fopen(argv[a],"r");
      if( input_file == NULL ) {
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

  srandom( seed );

  reset_agent( agt );

  if( feat != NULL ) {
    int i=0;
    if( feat->length[1] > agt->cd->num_reg ) {
      agt->cd->num_reg = feat->length[1];
      agt->reg=(Floating *)realloc(agt->reg,
                           agt->cd->num_reg*sizeof(Floating));
      printf("num_reg = %d\n",agt->cd->num_reg);
      fflush(stdout);
    }
    // load features into registers
    while( i < agt->cd->num_reg && i < feat->length[1] ) {
      agt->reg[i] = feat->val[feat->index[1]+i];
      i++;
    }
    free_channel(feat);
  }

  while( agt->running ) {

    if( command != 'g' ) {
      Snode *sf = agt->call_stack + agt->fp;
      if(    command == 's' || command == 'h' || command == 'b'
         ||( command == 'c' && sf->c == c && sf->k == k )
         ||( command == 'n' && agt->fp <= fp )) {
        if( command == 'h' || command == 'b' ) {
          command = 's';
        }
        else {
          print_state( agt, in, out, stdout );
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
      op = step( agt, in, out );
      if( op == INP ) {
        agt->bstack[agt->bp] = TRUE;
        if( fgets(line,MAX_LINE,input_file) <= 0 ) {
          if( a >= argc ) {
            agt->bstack[agt->bp] = FALSE;
          }
          else {
            fclose( input_file );
            input_file = fopen(argv[a],"r");
            if( input_file == NULL ) {
              printf("File not found: %s\n",argv[a]);
              exit(1);
            }
            if( fgets(line,MAX_LINE,input_file) <= 0 ) {
              agt->bstack[agt->bp] = FALSE;
            }
            a++;
          }
        }
        if( agt->bstack[agt->bp] ) {
          ch = getc( input_file );
          if( ch != '\n' ) {
            ungetc( ch,input_file );
            line[strlen(line)-1] = '\0';
          }
          clear_message(in,1); // overwrite previous input
          scan_next_input( in,line );
          agt->bstack[agt->bp] = fetch_input( in );
        }
      }
      else if( op == OUT ) {
        print_message( out,out->om,stdout );
        clear_message( out,1 );
      }
    }
  }
  fclose( input_file );
  if( command != 'g' ) {
    print_state( agt,in,out,stdout );
  }
  printf("Total number of steps = %lld\n",agt->step );

  free_channel( in );
  free_channel( out );
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
  fprintf(stderr,"Usage: %s [-d <seed>] [+r reg.in] [-v] file.hrc file1.in [file2.in ..]\n",argv_0);
  exit(1);
}


