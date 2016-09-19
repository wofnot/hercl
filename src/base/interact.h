/** \file interact.h
*/

#define  AGT_STEPS_PER_MOVE  10000

// These functions are provided in either base/inter_hercl.c or inter/xxx.c
void  *new_env_state();
Boolean    env_reset(void *void_env_state,Channel *in,Channel *out,int r);
Boolean env_continue(void *void_env_state,Channel *in,Channel *out);
Score *env_get_score(void *void_env_state);
void  free_env_state(void *void_env_state);

// These functions are provided in base/interact.c
void   tune_interact(Search_Param *par,long seed[],Channel *env_cfg,
                     Channel *agt_cfg,Ladder *lad);
void search_interact(Search_Param *par,long seed[],Channel *env_cfg,
                     Channel *agt_cfg,Ladder *lad,Library *lib);
void multi_instance_interact(Search_Param *par,Channel *env_cfg,Channel *agt_cfg,
                             Code *agt_code0,Library *lib);
void   test_interact(Candidate *can,void *env_state,Channel *env_cfg,Channel *agt_cfg,
                     Channel *env_out,Channel *agt_out,Channel *agt_out_prev,
                     int reactive,long seed[],int max_trials);

