/** \file message.c
*/
// TODO: Variable length messages?

#ifdef USE_MPI

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
#include "message.h"

#include <time.h>

#include <mpi.h>

Long cumul = 0;

int serialise_code(
               char *s,
               Code *cd,
               int   task_id
              )
{
  int k = string_from_int(s, task_id);
  k += string_from_code(s+k, cd);
  return k;
}

void deserialise_code(
               char *s,
               Code **cd,
               int  *task_id
              )
{  
  *task_id = int_from_codons(s,0);
  int k = skip_int(s);
  *cd = code_from_string( s+k ); 
}


// only includes entries after the code base.
int serialise_lib(
              char *s,
              Library *lib,
              Boolean which[MAX_LIB]
)
{
  int i,k;
  int q =0;
  k = 0;
  Code *cd;
  for (i = lib->code_base; i < lib->num_code; i++) {
    if (which[i]) {
      q++;
      cd = lib->code[i];
      if (cd != NULL) {
        // write this code (with index) as a line to s
        k += serialise_code(s+k, cd, i - lib->code_base);
        // replace null-terminator with newline char
        s[k] = '\n';
        k++;
      }
    }
  }
  //finish with null
  s[k] = '\0';
  return k;
}

// takes a string and updates the given library accordingly
void deserialise_lib(
              char *s,
              Library *lib
)
{
  // iterate through lines in string
  // deserialise each line and add to lib
  int i = 0;
  int k = 0; // start of current line
  int task_id;
  Code *cd;
  while (s[i] != '\0') {
    if (s[i] == '\n') {
      // deserialise this line
      s[i] = '\0'; // so that deserialise_code knows where to stop
      deserialise_code(s+k,&cd, &task_id);
#ifdef SHARING
      insert_library(lib,lib->code_base+task_id,cd);
#endif
      s[i] = '\n';
      k = i + 1; // start of next line
    }
    i++;
  }
}
              

#if NETWORK_TYPE == ALL_TO_ALL

/********************************************************//**
  check for waiting messages,
  update library with champs from
  other processes
  return true if a halt message has been received.
*/
Boolean update_library(Library *lib) {
  
  Boolean flag = FALSE;
  MPI_Status status;
  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
  
  while (flag) { // while messages are waiting
    // there is at least one message waiting, receive it and add the new champ
    char data[MAX_MSG_LEN];
    MPI_Recv(&data, MAX_MSG_LEN, MPI_CHAR,
              MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
              
    if (strncmp(data, "halt", 5) == 0) {
      return TRUE;
    } else {
      Code *cd;
      int task_id;
      if (!globals.no_share) {
        deserialise_code(data, &cd, &task_id);
        insert_library(lib,lib->code_base+task_id,cd);
      }
  
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);    
    }    
  }
  
  return FALSE;
}

/********************************************************//**
  check for waiting messages, while in trimming mode
  return true if a halt message has been received,
  and the sender of the halt has lower rank than this process
  (So that if multiple halts are sent, the lowest rank process continues)
*/
Boolean test_trim_halt(int sender_rank) {
    
  Boolean flag = FALSE;
  MPI_Status status;
  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
  
  while (flag) { // while messages are waiting
    // there is at least one message waiting, receive it and add the new champ
    char data[MAX_MSG_LEN];
    MPI_Recv(&data, MAX_MSG_LEN, MPI_CHAR,
              MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
              
    if ( (strncmp(data, "halt", 5) == 0)
      && (status.MPI_SOURCE < sender_rank) ) {
      // Only halt if the sender has lower rank than this process
      return TRUE;
    } else {
      // get the next message if there is one
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);    
    }    
  }
  
  return FALSE;
}

/********************************************************//**
  broadcast champ to all other processes
*/
void broadcast_code(Code *cd, int lib_index, int sender_rank, int world_size) {
   
  int dest;
  char message[MAX_MSG_LEN];
  serialise_code(message, cd, lib_index);
  for (dest = 0; dest < world_size; dest++) {
    if (dest == sender_rank) continue;
  
    MPI_Request request;
    MPI_Isend(message, MAX_MSG_LEN, MPI_CHAR, dest, 0, MPI_COMM_WORLD, &request);
    
    MPI_Wait(&request, MPI_STATUS_IGNORE);
  }
}

/********************************************************//**
  tell all processes to halt asap
*/
Boolean broadcast_halt(int sender_rank, int world_size) {
   
  int dest;
  char message[MAX_MSG_LEN];
  strncpy(message, "halt", 5);
  for (dest = 0; dest < world_size; dest++) {
    if (dest == sender_rank) continue;
    MPI_Request request;
    MPI_Isend(message, MAX_MSG_LEN, MPI_CHAR, dest, 0, MPI_COMM_WORLD, &request);
    
    MPI_Wait(&request, MPI_STATUS_IGNORE);
  }
  
  return TRUE;
}

#elif NETWORK_TYPE == STAR_NETWORK

// notify master process that worker is about to terminate (different from halting)
void signal_finished() {  
  MPI_Send(0, 0, MPI_CHAR, 0, TAG_FINISHED, MPI_COMM_WORLD);
}

Boolean update_library(Library *lib) {
  //struct timespec start,end;
  //printf("update library\n");
  // send request to master and update library according to result

  char data[MAX_MSG_LEN];

  MPI_Send(0, 0, MPI_CHAR, 0, TAG_REQUEST_UPDATE, MPI_COMM_WORLD);

  MPI_Status status;

  int msg_len;
  //clock_gettime( CLOCK_REALTIME, &start);
  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  //clock_gettime( CLOCK_REALTIME, &end);
  MPI_Get_count(&status, MPI_CHAR, &msg_len);
  //printf ("worker: msg_len = %d, tag = %d\n", msg_len, status.MPI_TAG);

  MPI_Recv(&data, msg_len, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  

  Boolean halt = FALSE;

  switch (status.MPI_TAG) {
  case (TAG_NO_UPDATE):
    // no action necessary
    break;
  case (TAG_HALT):
    halt = TRUE;
    break;
  case (TAG_UPDATE):
    // update library, will have no effect if data is empty
    
    if (!globals.no_share) {
      deserialise_lib(data, lib);
    }
    break;
  default:
    printf("update_library: Unexpected tag %d\n", status.MPI_TAG);
  }

  //Long t = ( end.tv_sec - start.tv_sec )*1000000000L + ( end.tv_nsec - start.tv_nsec );
  //cumul += t;
  //printf("t = %lld, cumul = %lld\n", t, cumul);
  
  return halt;
}

void broadcast_code(Code *cd, int lib_index, int sender_rank, int world_size) {
    
  char data[MAX_MSG_LEN];
  int k = serialise_code(data, cd, lib_index);
  
  MPI_Send(data, k+1, MPI_CHAR, 0, TAG_SUBMIT_CHAMP, MPI_COMM_WORLD);
  
}

Boolean broadcast_halt(int sender_rank, int world_size) {
   
  MPI_Status status;

  MPI_Send(0, 0, MPI_CHAR, 0, TAG_HALT, MPI_COMM_WORLD);
  MPI_Recv(NULL, 0, MPI_CHAR,0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  switch (status.MPI_TAG) {
  case (TAG_HALT):
    // someone else already halted
    return FALSE;
    break;
  case (TAG_HALT_ACK):
    break;
  default:
    printf("broadcast_halt: Unexpected tag %d\n", status.MPI_TAG);
  }
  return TRUE;
}

#endif

#endif
