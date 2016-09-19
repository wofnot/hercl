/** \File master.c
// experimental message passing system where this node is responsible
// for hosting the master copy of the shared library and communicating with worker processes
*/


#include <stdio.h>
#include <mpi.h>
#include <assert.h>

#include "step.h"
#include "scan_print.h"
#include "point.h"
#include "cross_mutate.h"
#include "eval.h"
#include "message.h"


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);


  char tags[6][100] = {"request", "submit", "halt", "update", "no update", "finished"};
  
  char in[MAX_MSG_LEN];
  char out[MAX_MSG_LEN];

  Code *champ;
  int task_id;

  int i,j;

  int num_processes;
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  // num_processes includes this one but it doesn't really matter

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  assert(my_rank == 0); // Master must be run on process 0

  Boolean outdated[num_processes][MAX_LIB];
  // outdated[i][j] means process i does not have the latest version of library entry j
  Boolean finished[num_processes];

  for (i = 1; i < num_processes; i++) {
    finished[i] = FALSE;
    for (j = 0; j < MAX_LIB; j++){
      outdated[i][j] = FALSE;
    }
  }

  Boolean halt = FALSE;
  
  Library *lib = new_library(MAX_LIB);
  
  MPI_Status status; 

  int msg_len;

  Boolean running = TRUE;

  while (running) {
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, &msg_len);
    //printf ("master: msg_len = %d, tag = %d\n", msg_len, status.MPI_TAG);
    MPI_Recv(&in, msg_len, MPI_CHAR,
             status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    //printf("Master received: %s from %d.\n", tags[status.MPI_TAG], status.MPI_SOURCE);
    switch (status.MPI_TAG) {
    case (TAG_REQUEST_UPDATE): // can ignore data
      // Received an update request from worker node
      if (halt) {
        // relay halt message
        MPI_Send(0, 0, MPI_CHAR, status.MPI_SOURCE, TAG_HALT, MPI_COMM_WORLD);
      } else {
        // send update (serialise_lib will return "\0" if no updates needed)
        j = serialise_lib(out, lib, outdated[status.MPI_SOURCE]);
        MPI_Send(out, j+1, MPI_CHAR, status.MPI_SOURCE, TAG_UPDATE, MPI_COMM_WORLD);
        for (j = 0; j < MAX_LIB; j++){
          outdated[status.MPI_SOURCE][j] = FALSE;
        }
      }
      break;

    case (TAG_FINISHED):
      finished[status.MPI_SOURCE] = TRUE;
      // if all the workers have finished, master can terminate.
      running = FALSE;
      for (i = 1; i < num_processes; i++) {
        if (!finished[i]){
          running = TRUE;
        }
      }
      break;

    case (TAG_SUBMIT_CHAMP):
      if (!halt) {
        // Update master copy of library.
        deserialise_code(in, &champ, &task_id);

        insert_library(lib,task_id,champ);

        // set other processes to outdated
        for (i = 0; i < num_processes; i++) {
          if (i != status.MPI_SOURCE) {
            outdated[i][task_id] = TRUE;
          }
        }
      }
      break;

    case (TAG_HALT): // ignore data
      if (halt) {
        // someone else has already halted, send HALT
        MPI_Send(0, 0, MPI_CHAR, status.MPI_SOURCE, TAG_HALT, MPI_COMM_WORLD);
      } else {
        halt = TRUE;

        // acknowledge halt
        MPI_Send(0, 0, MPI_CHAR, status.MPI_SOURCE, TAG_HALT_ACK, MPI_COMM_WORLD);
      }
      break;
    default:
      printf("Master: Unknown tag %d\n", status.MPI_TAG);
    }
    // printf("message dealt with\n");
  }
  printf("Master terminating\n");
  MPI_Finalize();
  return 0;
}
