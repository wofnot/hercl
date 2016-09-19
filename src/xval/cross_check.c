/** \file cross_check.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../base/step.h"
#include "../base/scan_print.h"
#include "../base/cross_mutate.h"
#include "../base/eval.h"
#include "../base/super.h"

#define   MAX_MESSAGE   1024
#define   MAX_BUF      16386

/********************************************************//**
   Print a usage message and exit.
*/
void usage( char *argv_0 )
{
  fprintf( stderr,"Usage: %s <super.in> <num_fold> <task.fold> <num_run> <task.hrc>\n",argv_0);
  exit(1);
}

int * scan_folds( FILE *fp, int num_items )
{
  int *fold;
  int a,r,s;

  fold=(int *)malloc((num_items+1)*sizeof(int));
  check_null(fold);

  for( s=1; s <= num_items; s++ ) {
    fscanf( fp,"fold[%d] = %d\n", &r, &a );
    if( r != s ) {
      printf("HEY!\n");
      exit( 1 );
    }
    else {
      fold[r] = a;
    }
  }
  return( fold );
}

void print_code_by_miss(
                        int     num_fold,
                        int     num_run,
                        Code ***code,
                        int   **miss
                       )
{
  int prev_miss = -1;
  int min_miss = 0;
  int a,r;

  while( min_miss < 1000000 ) {
    min_miss = 1000000;
    for( a=0; a < num_fold; a++ ) {
      for( r=1; r <= num_run; r++ ) {
        if( miss[a][r] > prev_miss && miss[a][r] < min_miss ) {
          min_miss = miss[a][r];
        }
      }
    }
    for( a=0; a < num_fold; a++ ) {
      for( r=1; r <= num_run; r++ ) {
        if( miss[a][r] == min_miss ) {
          printf("miss = %d\n",min_miss);
          print_cell( code[a][r],0,stdout );
        }
      }
    }
    prev_miss = min_miss;
  }
}

int main( int argc, char *argv[] )
{
  char line[256];
  FILE      *fp;
  Super     *sup;
  Channel   *out;
  Candidate *can;
  Code    ***code;
  //int     ***item;
  //int      **num_items;
  int       *fold;
  int       *run;
  Floating **val;
  Floating min_val,max_val,mid_val;
  Floating target_val;
  Code *cd;
  double sum;
  int mismatch,pos_count,neg_count;
  int num_fold;
  int num_run;
  int rmax,r0;
  int a,j,n,r,s;

  if( argc < 6 ) {
    usage( argv[0] );
  }

  init_codons();

  fp = fopen_check( argv[1],"r" ); // super.in
  sup = scan_super( fp );
  fclose( fp );

  if(   sup->inputs_per_item == 1
     && sup->inputs_same_length && sup->cost_type != LED ){
    sup->register_features = TRUE;
    /*
    if( sup->input->length[1] > 10 ) {
      cd->num_reg = sup->input->length[1];
      printf("num_reg = %d\n",cd->num_reg);
      fflush(stdout);
    }
    */
  }

  min_val = 0.0;
  max_val = 0.0;
  for( s=1; s <= sup->total_items; s++ ) {
    target_val = sup->target->val[sup->target->index[s]];
    if( target_val < min_val ) {
      min_val = target_val;
    }
    if( target_val > max_val ) {
      max_val = target_val;
    }
  }
  mid_val = 0.5 *( min_val + max_val );

  num_fold = atoi( argv[2] );
  num_run  = atoi( argv[4] );

  val=(double **)malloc((sup->total_items+1)*sizeof(double *));
  check_null( val );
  for( s=1; s <= sup->total_items; s++ ) {
    val[s] =(double *)malloc((num_run+1)*sizeof(double));
  }
  for( s=1; s <= sup->total_items; s++ ) {
    val[s][0] = sup->target->val[sup->target->index[s]];
  }

  fp = fopen_check( argv[3],"r" );
  fold = scan_folds( fp,sup->total_items );
  fclose( fp );

  code = (Code ***)malloc( num_fold*sizeof(Code **));
  for( a=0; a < num_fold; a++ ) {
    code[a] = (Code **)malloc((num_run+1)*sizeof(Code *));
  }

  //printf("opening %s..\n", argv[5] );
  fp = fopen_check( argv[5],"r" );

  for( r=1; r <= num_run; r++ ) {
    for( a=0; a < num_fold; a++ ) {

      while( fgets( line,256,fp ) > 0 && line[0] != '=' )
        ;

      j = 0;
      while( line[j++] != 'A' )
        ;
      if( a != line[j]-'0' || r != line[j+1]+1-'a' ) {
        printf("Error! fold:%d=%d? run:%d=%d?\n",
               a,line[j]-'0',r,line[j+a]+1-'a');
        exit(1);
      }
      printf("%d %d\n",r,a);
      fflush(stdout);
      //fgets( line,256,fp );
      cd = scan_code( fp );
      //print_code( cd, stdout );
      code[a][r] = cd;
    }
  }
  fclose( fp );

  for( a=0; a < num_fold; a++ ) {
    for( r=1; r <= num_run; r++ ) {
      printf("fold=%d run=%d\n",a,r);
      print_code( code[a][r], stdout );
    }
  }

  can = new_candidate();
  out = new_channel( MAX_MESSAGE-2,MAX_BUF );

  sup->test_items = 1;

  /*
  for( s=1; s <= sup->total_items; s++ ) {
    a = fold[s];
    for( r=1; r <= num_run; r++ ) {
      can->agt->cd = code[a][r];
      sup->test_set[1] = s;
      reset_score( can );
      eval_super( can,sup,NULL,out,sup->test_set,1,1 );
      if( can->score.reject || out->om != 1 ) {
        val[s][r] = mid_val;
      }
      else {
        val[s][r] = out->val[out->index[1]];
      }
      if(  ( val[s][0] < mid_val && val[s][r] >= mid_val )
         ||( val[s][0] > mid_val && val[s][r] <= mid_val )) {
        miss[a][r]++;
      }
    }
    //printf("\n");
    fflush( stdout );
  }
  //print_code_by_miss( num_fold,num_run,code,miss );
  */

  for( s=1; s <= sup->total_items; s++ ) {
    sup->test_set[1] = s;
    for( r=1; r <= num_run; r++ ) {
      can->agt->cd = code[fold[s]][r];
      reset_score( can );
      eval_super( can,sup,NULL,out,sup->test_set,1,1 );
      if( can->score.reject || out->om != 1 ) {
        val[s][r] = mid_val;
      }
      else {
        val[s][r] = out->val[out->index[1]];
      }
    }
    //printf("\n");
    fflush( stdout );
  }

  for( s=1; s <= sup->total_items; s++ ) {
    if( val[s][0] < mid_val ) {
      printf("%3d %d ",s,(int)val[s][0]);
      for( r=1; r <= num_run; r++ ) {
        if( val[s][r] <= min_val ) {
          printf("v");
        }
        else if( val[s][r] >= max_val ) {
          printf("^");
        }
        else if( val[s][r] == mid_val ) {
          printf(" ");
        }
        else if( val[s][r] > mid_val ) {
          printf("'");
        }
        else {
          printf(",");
        }
      }
      printf("\n");
    }
  }
  for( s=1; s <= sup->total_items; s++ ) {
    if( val[s][0] > mid_val ) {
      printf("%3d %d ",s,(int)val[s][0]);
      for( r=1; r <= num_run; r++ ) {
        if( val[s][r] <= min_val ) {
          printf("v");
        }
        else if( val[s][r] >= max_val ) {
          printf("^");
        }
        else if( val[s][r] == mid_val ) {
          printf(" ");
        }
        else if( val[s][r] > mid_val ) {
          printf("'");
        }
        else {
          printf(",");
        }
      }
      printf("\n");
    }
  }

  run = (int *)malloc((sup->total_items+1)*sizeof(int));
  check_null( run );
  /*
  printf("HARD VOTE:\n");
  for( rmax=1; rmax <= num_run; rmax++ ) {
    mismatch = 0;
    for( s=1; s <= sup->total_items; s++ ) {
      pos_count = 0;
      neg_count = 0;
      for( r=1; r <= rmax; r++ ) {
        if( val[s][r] >= max_val ) pos_count++;
        if( val[s][r] <= min_val ) neg_count++;
      }
      if(  ( val[s][0] < mid_val && pos_count >= neg_count )
         ||( val[s][0] > mid_val && pos_count <= neg_count )) {
        //printf(" %d",s);
        mismatch++;
      }
    }
    //printf("\n");
    printf("%2d -> %3d\n",rmax,mismatch);
  }
  */

  //  printf("SOFT VOTE:\n");
  for( rmax=1; rmax <= num_run; rmax ++ ) {
    mismatch = 0;
    for( n=0; n < 10000; n++ ) {
      r = 0;
      while( r < rmax ) {
        run[r] = 1 + random()% num_run;
        for( r0=0; r0 < r && run[r0] != run[r]; r0++ )
          ;
        if( r0 == r ) {
          r++;
        }
      }
      for( s=1; s <= sup->total_items; s++ ) {
        sum = 0.0;
        pos_count = 0;
        neg_count = 0;
        for( r=0; r < rmax; r++ ) {
          sum += val[s][run[r]] - mid_val;
          if( val[s][run[r]] > mid_val ) pos_count++;
          if( val[s][run[r]] < mid_val ) neg_count++;
        }
        if( pos_count > neg_count ) {
          sum = 1.0;
        }
        else if( neg_count > pos_count ) {
          sum = -1.0;
        }
        if(  ( val[s][0] < mid_val && sum >= 0.0 )
           ||( val[s][0] > mid_val && sum <= 0.0 )) {
          //printf(" %d",s);
          mismatch++;
        }
      }
    }
    //printf("\n");
    printf(" %1.1f",mismatch/10000.0);
  }
  printf("\n");

  /*
  printf("AGGREGATING:\n");
  for( rmax=1; rmax <= num_run; rmax++ ) {
    mismatch = 0;
    for( s=1; s <= sup->total_items; s++ ) {
      val_sum = 0.0;
      for( r=1; r <= rmax; r++ ) {
        if( val[s][r] >= max_val ) {
          val_sum += ( max_val - mid_val );
        }
        if( val[s][r] <= min_val ) {
          val_sum += ( min_val - mid_val );
        }
        else {
          val_sum += ( val[s][r] - mid_val );
        }
      }
      if(  ( val[s][0] < mid_val && val_sum >= 0.0 )
         ||( val[s][0] > mid_val && val_sum <= 0.0 )) {
        //printf(" %d",s);
        mismatch++;
      }
    }
    //printf("\n");
    printf("%2d -> %3d\n",rmax,mismatch);
  }
  */
  /*
    if( out->length[1] != 1 ) {
      mismatch++;
    }
    else {
      output_val = out->val[out->index[1]];
      if(  ( target_val <= 0.5 && output_val >= 0.5 )
         ||( target_val >= 0.5 && output_val <= 0.5 )) {
        mismatch++;
      }
    }
  }
  printf("%d errors out of %d items\n",mismatch,sup->total_items);
  */

  return 0;
}
