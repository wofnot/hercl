/** \File eval.h
*/

#define   RESHUFFLE_FAILS  6

//#define  RESERVE_CHAMP  1

#ifdef   RESERVE_CHAMP
#define  CHAMP          1
#else
#define  CHAMP          0
#endif

//#define  MAX_EVAL     900000000
//#define  MAX_EVAL     99900000000

//#define  MAX_TRACE        65536
//#define  PRINT_HISTOGRAM      1

//#define  ZERO_SEED            1

#define  MAX_SEEDS         1024
#define  MAX_LIB            256

#define  TRAIN         0
#define  TEST          1

#define  FAILED        0
#define  SUPERIOR      1
#define  FASTER        2
#define  SHORTER       3
#define  EQUAL         4
#define  LONGER        5
#define  SLOWER        6
#define  INFERIOR      7
#define  PENALTY       8
#define  REJECT        9

#define  COST_FACTOR  1000.0
#define  TIME_FACTOR     0.01
#define SPACE_FACTOR     0.01

#define  MAX_CHILD   65535

#define  SUPER         0
#define  INTERACT      1
#define  PREDICT       2

typedef struct search_param Search_Param;
typedef struct score     Score;
typedef struct candidate Candidate;
typedef struct library   Library;
typedef struct ladder    Ladder;

struct search_param {
  FILE *file_in;
  char *filename;
  char *output_dir;
  Boolean is_tune;  // tuning instead of evolving
  Boolean reactive; // only used for INTERACT
  Boolean is_multi; // multi-search
  Boolean silent;
  Boolean incremental;
  int min_trials;
  int max_trials;
  int num_instances;
  int task; // SUPER, INTERACT or PREDICT
  int max_epoch;
  int trim_interim;
  int trim_final;
  int num_cells;
  int cost_type;
  double item_cost;
};

struct score {
  int     num_trials;
  int     num_steps;
  int     mismatch;
  Boolean all_outputs_same;
  Boolean successful;
  Boolean reject;
  double  penalty;
  double  cost;
};

struct candidate {
  Agent *agt;
  int    status;
  Score  score;
};

struct library {
  Code  **code;
  int max_code;
  int num_code;
  int code_base;
  int index; // index of most recently inserted item
};

struct ladder {
  Candidate *can[MAX_LEVEL];
  Library  *bank[MAX_LEVEL];
  Boolean   adapting;
  // maximum children candidate can have before FAILED
  int       max_children[MAX_LEVEL];
  // number of children of candidate currently at step s
  int       num_children[MAX_LEVEL];
  // number of descendents of candidate currently at step s
  int    num_descendants[MAX_LEVEL];
  // total number of candidates at step s, level m
  Long  total_candidates[MAX_LEVEL][MAX_LEVEL][3];
  // number of descendants of candidates at step s, level m
  Long total_descendants[MAX_LEVEL][MAX_LEVEL][3];
  // number of superior or initially inferior candidates
  int      total_sup_inf[MAX_LEVEL][MAX_LEVEL][3];
  // number of (eventually) superior candidates at step s, level m
  int     total_superior[MAX_LEVEL][MAX_LEVEL][3];
#ifdef PRINT_HISTOGRAM
  int hist[MAX_LEVEL][MAX_LEVEL][3][10][17];
#endif
#ifdef MAX_TRACE
  char  ts[MAX_TRACE];
  char *ti[MAX_LEVEL];
#endif
  Long ncomp; // total number of fitness evaluations
  Long neval; // total number of items evaluated
  int max_bank;
  int max_child_per_epoch;
  int num_fails; // used for reshuffling data
  //int epoch;
  int s;
};

extern long eval_seed;

Candidate * new_candidate();
Boolean  better_candidate(Candidate *b,Candidate *a,
                          double target_cost,Boolean final);
void free_candidate(Candidate *can );
void    reset_score(Candidate *can);
void    print_score(Candidate *can,FILE *fp);
int   perfect_score(Candidate *can,double target_cost);
void   update_score(Candidate *can,Channel *out,Channel *tgt,int cost_type);

Library      *new_library(int max_code);
void scan_code_to_library(Library *lib,char *libname);
void       insert_library(Library *lib,int index,Code *cd);
void        print_library(Library *lib,FILE *fout);
void        clear_library(Library *lib );
void     move_to_codebank(Ladder *lad,int s);

Ladder     * new_ladder(int max_bank);
void       print_ladder(Ladder *lad,FILE *fp);
void        free_ladder(Ladder *lad);
void       reset_ladder(Ladder *lad,Code *cd);
void   adjust_max_child(Ladder *lad,int num_trials,int max_trials);
Boolean check_reshuffle(Ladder *lad,int min_trials,int prev_trials,
                        int num_trials,int max_trials);
void     reset_codebank(Ladder *lad);
void    random_codebank(Ladder *lad,int num_code);
void          procreate(Ladder *lad,Library *lib);
void              breed(Ladder *lad,Code *cd1,Level level);
Candidate          *top(Ladder *lad);
Candidate          *pop(Ladder *lad);
int            cull_top(Ladder *lad,FILE *fp);
void    top_replace_pop(Ladder *lad,FILE *fp);
#ifdef PRINT_HISTOGRAM
void         print_hist(Ladder *lad,FILE *fo);
#endif
void        print_stats(Ladder *lad,FILE *fo);

void     generate_seeds(int max_seeds,long seed[]);
void      shuffle_seeds(int max_seeds,long seed[]);
void     bring_to_front(int rank[],int n);

void  print_termination(Ladder *lad,Long ncomp,Long neval,
                        int epoch,int task,FILE *fo);

//void compress_ladder(Ladder *lad );
