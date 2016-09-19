/** cross_mutate.h */

#define  DIGIT     0
#define  REMOVE    1
#define  EXCHANGE  2
#define  INSERT    3

//int mlevel_prior();
Code * random_code(Code *cd0,int num_cells);
Code *cross_mutate(Code *cd0,Code *cd1,Level level);
