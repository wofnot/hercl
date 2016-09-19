/** \File message.h
*/

#define SHARING  1

#define ALL_TO_ALL 0
#define STAR_NETWORK 1

#define NETWORK_TYPE STAR_NETWORK

#if NETWORK_TYPE == STAR_NETWORK
  #define MAX_MSG_LEN (100000)

  #define TAG_REQUEST_UPDATE 0
  #define TAG_SUBMIT_CHAMP 1
  #define TAG_HALT 2
  #define TAG_UPDATE 3
  #define TAG_NO_UPDATE 4
  #define TAG_FINISHED 5
  #define TAG_HALT_ACK 6

#else
  #define MAX_MSG_LEN (2000)
#endif



// serialise and deserialise functions
int serialise_code(char *s, Code *cd, int task_id);
void deserialise_code(char *s, Code **cd, int  *task_id);
int serialise_lib(char *s, Library *lib, Boolean which[MAX_LIB]);
void deserialise_lib(char *s, Library *lib);

// interface for interprocess communication
// implementation depends on the chosen topology
#if NETWORK_TYPE == STAR_NETWORK
void signal_finished();
#elif NETWORK_TYPE == ALL_TO_ALL
Boolean test_trim_halt(int sender_rank);
#endif

Boolean update_library(Library *lib);
void broadcast_code(Code *cd, int lib_index, int sender_rank, int world_size);
Boolean broadcast_halt(int sender_rank, int world_size);
