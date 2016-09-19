#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <opencv2/opencv.hpp>

#include "../base/step.h"
#include "../base/eval.h"
#include "../base/scan_print.h"

#include "../metrics/metrics.h"

#include "art-school.h"

#include "../hashmap/hashmap.h"

//#define ALL_FRAMES

using namespace std;
using namespace cv;

double run_critic(Agent *critic, Channel *chan, metrics_t metrics) {
  // run the critic
  
  reset_agent( critic );
  
  // load metrics into critic registers
  
      
  int i;
  for (i = 0; i < (int) metrics.size(); i++) {
    assert(i < critic->cd->num_reg); // critic must have enough registers
    critic->reg[i] = metrics[i];
  }
  
  clear_message(chan, 1);
  
  int op = NONE;
  while( op != OUT && critic->running ){
    //Note that this will loop indefinitely if agent provided does not output/terminate
    op = step( critic, chan, chan );
  }
  
  if(op == OUT){
    cue_message(chan, chan->om);
  }
  
  // print_channel(es->chan, stdout);
  
  double score;
  
  //get response from agent-to-mimic and set target to this value
  if(!fetch_input(chan)){
    printf("Error: critic failed to output message\n");
    exit(1);
  }
  if(!scan_value(chan,&score)){
    printf("Error: critic failed to output item\n");
    exit(1);
  }
  
  return score;
}

double score_state(EnvState *es) {
  double cost = 0;
  
  metrics_t metrics;
  
  metrics = getMetrics(getMatrix(es->canvas), es->use_colour);

  double critic_score;
  int i;
  for (i = 0; es->critics[i] != NULL; i++) {
    critic_score = run_critic(es->critics[i], es->chan, metrics);
    critic_score = MAX(0, MIN(critic_score, 1.0)); //clamps score in [0,1]
    
   // printf("Critic %d gives %lf x %lf\n", i, 100*(1-critic_score), es->weights[i]);
    
    cost += es->weights[i]*100*(1-critic_score);
    
  }
  
  // 1 represents success
  
  
  return cost;
}

Agent *load_agent(char *fn) {

    Agent *agt = new_agent();
    FILE *fp = fopen_check( fn ,(char *)"r" );
    agt->cd = scan_code( fp );
    fclose(fp);
    return agt;
}

void load_settings(EnvState *es) {
  FILE *fp = fopen_check( (char *) SETTINGS_FILE ,(char *)"r" );
  
  int n;
  
  char fn[256];
  char c;
  
  // critic_fn
  n = fscanf(fp, "%c %s", &c, fn);
  assert(n == 2);
  
  if (c == 'l' || c == 'L') {
    // library of critics
    FILE *fp = fopen_check( fn ,(char *)"r" );
    
    int i = 0;
    // load critics one by one
    // format: <weight> <filename>
    while( fscanf( fp,"%lf %s", &(es->weights[i]), fn ) > 0 ) {
      es->critics[i] = load_agent(fn);
      i++;
    }
    fclose(fp);
    
    es->critics[i] = NULL;
  } else {
    // single critic
    es->critics[0] = load_agent(fn);
    es->weights[0] = 1.00;
    
    es->critics[1] = NULL;
  }
  
  // output_fn
  n = fscanf(fp, "%s", es->output_fn);
  assert(n == 1);
  if (es->output_fn[0] == '-' && es->output_fn[1] == 0) {
    es->output_fn[0] = 0;
  }
  
  // tolerance
  n = fscanf(fp, "%lf", &es->tolerance);
  assert(n == 1);
  
  // pushpop
  n = fscanf(fp, "%d", &es->pushpop);
  assert(n == 1);
  
  // no_loop
  n = fscanf(fp, "%d", &es->no_loop);
  assert(n == 1);

  // max_step
  n = fscanf(fp, "%d", &es->max_step);
  assert(n == 1);
  
  // max_movement
  n = fscanf(fp, "%lf", &es->max_movement);
  assert(n == 1);
  
  // use_colour
  n = fscanf(fp, "%d", &es->use_colour);
  assert(n == 1);
  
  // cout << "es->use_colour = " << es->use_colour << "\n";
  
  fclose(fp);
}

/********************************************************//**
   Allocate a new EnvState
*/
void * new_env_state()
{
  EnvState *es = (EnvState *)malloc(sizeof(EnvState));
  check_null( es );
  
  memset(es,0,sizeof(EnvState));
  
  es->sp = 0;
  
  load_settings(es);
  
  es->canvas = newCanvas(128,128, es->use_colour);
  
  if (es->no_loop) {
    is_allowed[(int)BRB] = FALSE;
  }
  
  // Allocate channels for this critic
  es->chan = new_channel();
  
  es->key_index = 0;
  
  es->cost_map = hashmap_new();

  return((void *)es);
}

int free_hashmap_data (char *key, any_t data) {
  //free(key);
  //free(data);
  return MAP_OK;
}

/********************************************************//**
   Free state occupied by env_state
*/
void free_env_state( void *void_es )
{
  EnvState *es = (EnvState *)void_es;
  
 //hashmap_iterate(es->cost_map, free_hashmap_data);
  
  hashmap_free(es->cost_map);
  
  // Free critic
  int i;
  for (i = 0; es->critics[i] != NULL; i++) {
    free_code( es->critics[i]->cd );
    compress_agent( es->critics[i] );
    free( es->critics[i] );
  }
  
  // Free channels
  free_channel(es->chan);
  
  freeCanvas(es->canvas);
  free( es );
}

/********************************************************//**
   Initialize the environment
*/
Boolean env_reset(
                  void    *void_es,
                  Channel *in,
                  Channel *out,
                  int r
                 )
{

  EnvState *es = (EnvState *)void_es;
 // printf(" ** env_reset **\n");
  es->score.penalty = FALSE;
  es->score.reject = FALSE;
  
  // reset history
  es->history[0] = 0;
  
  resetCanvas(es->canvas);
  
  es->step = 0;
  es->angle = 0;
  
  return( TRUE );
}

/********************************************************//**
   Continue running environment.
   Return TRUE if the trial is still running, FALSE otherwise
*/

Boolean env_continue(
                     void    *void_es,
                     Channel *in,
                     Channel *out
                    )
{
  double dmsg[100];
  int msg[100];
  int c;

  EnvState *es = (EnvState *)void_es;
 // printf("env_cont %d\n", es->step);
  
  
#ifdef DEBUG
  Scalar colour;
  getPenColour(es->canvas, &colour);

  cout << "current colour = " << colour[0] << " " << colour[1] << " " << colour[2] << "\n" ;
#endif
  
  Boolean cont = TRUE;
  
  es->score.successful = FALSE;
  
  double numPixels = getWidth(es->canvas)*getHeight(es->canvas); // upper bound on input messages
  
  if( !fetch_input(in)) {
    if (es->step == 0) { // only punish if it's the first step
      // *** No input ***
      es->score.cost = 1000.0;
      es->score.penalty = 1.0;
    }
    cont = FALSE;
  } else {
    c = 0;
    while (scan_value(in,&dmsg[c]) && c<5) {
    // FIX: not bounding message length causes segfault in step.c, should probs look into that
      if (!isfinite(dmsg[c])) {
        c = -1;
        break;
      }
      if (dmsg[c] > numPixels) {
        // clamp
        dmsg[c] = numPixels; // bit hacky, maybe just punish overflow instead of being so lenient?
      }
      msg[c] = ((int) dmsg[c])%500;
      
     // printf("msg[%d] = %d\n", c, msg[c]);
      
      c++;
    }
    if (c == -1) { // non finite
      es->score.cost = 800.0;
      es->score.penalty = 1.0;
      cont = FALSE;
    } else if(c == 0) {
      // *** empty message ***
      es->score.cost = 900.0;
      es->score.penalty = 1.0;
      cont = FALSE;
    } else {
      // parse message
      
      // restrict to valid command
      if (es->pushpop) {
        msg[0] = abs(msg[0]%CMD_COUNT);
      } else {
        msg[0] = abs(msg[0]%PUSH_TURTLE);
      }
      
      Command cmd = (Command) msg[0];
      
      if (cmd == TOGGLE) {
        togglePen(es->canvas);
        msg[1] = 0;
      } else if (cmd == FORWARD) {
        if (c < 2) {
          msg[1] = 1;
        }
        if (es->max_movement > 0 && isPenDown(es->canvas)) {
          if (abs(msg[1]) > es->max_movement) {
            msg[1] = ((msg[1] < 0) ? -1 : 1) * es->max_movement;
          }
        }
        movePen(es->canvas, msg[1]*sin(es->angle), msg[1]*cos(es->angle));
      } else if (cmd == TURN) {
        if (c < 2) {
          msg[1] = 10; // degrees
        }
        msg[1] = msg[1]%360;
        double delta = msg[1]*M_PI/180;
        es->angle = fmod(es->angle+delta + 2*M_PI, 2*M_PI);
      } else if (cmd == SET_SIZE) {
        if (c < 2) {
          msg[1] = 1;
        }
        msg[1] = abs(msg[1]%4);
        setPenSize(es->canvas, msg[1]); // cap size at 4.
        
      } else if (cmd == SET_COLOUR) {
        if (es->use_colour) {
          // set unspecified args to 0
          Scalar colour;
          
          getPenColour(es->canvas, &colour);
          
          
          // take arguments from artist in this order: l, h, s.
          // but the colour values are in the order: h, l, s.
          // corresponding permutation:
          int perm[3] = {1,0,2};
          
          int i;
          int ub = c - 1;
          if (c-1 > 3) {
            ub = 3;
          }
          for (i = 0; i < ub; i++) {
         //   cout << i << " " << c << "\n";
            colour[i] = abs(msg[perm[i] + 1]%256);
          }
          
   
#ifdef DEBUG
          cout << "set colour = " << colour[0] << " " << colour[1] << " " << colour[2] << "\n" ;
#endif
           
          setPenColour(es->canvas, colour);
        } else {
          if (c < 2) {
            msg[1] = 0;
          }
          msg[1] = abs(msg[1]%255);
          setPenColour(es->canvas, msg[1]);
        }
        
      } else if (cmd == PUSH_TURTLE) {
        if (es->sp < STACK_SIZE) {
          // write current settings to top of stack, increment sp
          turtleFrame *current = &(es->turtleStack[es->sp]);
          getPenPos(es->canvas, &(current->penX), &(current->penY));
          getPenColour(es->canvas, &(current->penColour));
          getPenSize(es->canvas, &(current->penSize));
          current->angle = es->angle;
          es->sp++;
        }
        
      } else if (cmd == POP_TURTLE) {
        // decrement sp, load settings from top of stack
        if (es->sp > 0) {
          es->sp--;
          
          turtleFrame current = es->turtleStack[es->sp];
         
          if (isPenDown(es->canvas)) {
            togglePen(es->canvas);
            setPenPos(es->canvas, current.penX, current.penY);
            togglePen(es->canvas);
          } else {
            setPenPos(es->canvas, current.penX, current.penY);
          }
          
          setPenColour(es->canvas, current.penColour);
          setPenSize(es->canvas, current.penSize);
          es->angle = current.angle;
        }
        
      } else {
        printf("Unknown command: %d\n", msg[0]);
        assert(0); // shouldn't happen
      }
      
      // write latest move
      char this_move[300];
      snprintf(this_move, MAX_HISTORY, "%d %d;",msg[0],msg[1]);
      
      int l = strlen(es->history);
      int l2= strlen(this_move);
      if (l+l2+1 > MAX_HISTORY) {
        fprintf(stderr, "error: MAX_HISTORY too small\n");
        fprintf(stderr, "%s%s\n",es->history, this_move);
      }
      
      strncat(es->history, this_move, l2);
      
     // printf("history: %s\n", es->history);
    }
  }
  
  es->step++;
  
  if (es->step > es->max_step) {
    cont = FALSE;
  }
  
  if (!cont && es->score.penalty < 0.01) {
  
    double *x;
  
    // check to see if this sequence has already been evaluated:
    
    if (hashmap_get(es->cost_map, es->history, (void **) &x) == MAP_OK)  {
      // already been evaluted, use previous score
      es->score.cost = *x;
      //printf("get\n");
    } else {
      // evaluate and add to es->cost_map

      //printf("put\n");
      strncpy(es->keys[es->key_index], es->history, MAX_HISTORY);

      es->score.cost = score_state(es);
      es->values[es->key_index] = es->score.cost;
      int e = hashmap_put(es->cost_map, es->keys[es->key_index], (void *) &(es->values[es->key_index]));
      if (e == MAP_OMEM) {
        fprintf(stderr,"Hashmap out of memory\n");
        exit(0);      
      }
      
      es->key_index = (es->key_index + 1)%NUM_KEYS;
      
      if(es->key_index == 0) {
        fprintf(stderr,"key array used up (not an error)\n");
      }
    }
     // es->score.cost = score_state(es);

    if (es->score.cost < es->tolerance*100) {
      es->score.successful = TRUE;
    }
    if (es->output_fn[0] != 0
#ifndef DEBUG
        && es->score.successful // if it's debug mode print it anyway
#endif
    ) {
        // save image 
      fprintf(stderr, "saving %s (after %d steps).\n", es->output_fn,es->step);
      savePNG(es->canvas, es->output_fn);
      
    }
    
#ifdef ALL_FRAMES
    //showCanvas(es->canvas);
    char fn[100];
    fprintf(stderr, "Step %d concluded\n", es->step);
    sprintf(fn, "pics/painter/%d.png", es->step);
    savePNG(es->canvas, fn);
#endif

  }



  return(cont); // terminate after second iteration
}

/********************************************************//**
   Return the score assigned to the agent by the environment
*/
Score *env_get_score( void *void_es )
{
  EnvState *es = (EnvState *)void_es;
  // printf("-----score returned-----\n");
  //printf("cost: %f \n", es->score.cost);
  //printf("penalty: %f \n", es->score.penalty);
  // printf("reject: %d \n", es->score.reject);
  return(&(es->score));
}
