/** super.h
*/

typedef struct super Super;

struct super {
  Channel *input;
  Channel *target;
  double   item_cost;
  int      cost_type;
  int      inputs_per_item;
  Boolean  inputs_same_length;
  Boolean  register_features;
  int *train_set;
  int *test_set;
  int  train_items;
  int  test_items;
  int  total_items;
};

void   tune_super(Search_Param *par,Super *sup,Channel *agt_cfg,
                  Channel *out,Ladder *lad);
void search_super(Search_Param *par,Super *sup,Channel *agt_cfg,
                  Channel *out,Ladder *lad,Library *lib);
void  multi_super(Search_Param *super_par,Code *agt_code0,Library *lib);
void   test_super(Candidate *can,Super *sup,Channel *agt_cfg,
                  Channel *out,int set);
void   eval_super(Ladder *lad,Candidate *can,Super *sup,
                  Channel *agt_cfg,Channel *out,int *tset,int r,int o);
Super  *new_super();
Super *scan_super(FILE  *fi);
void  print_super(Super *sup,FILE *fo);
void   free_super(Super *sup);
Boolean all_items_same(Channel *chan);

