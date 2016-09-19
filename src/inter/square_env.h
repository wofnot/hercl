/** \File square_env.h
*/

typedef struct envstate EnvState;

struct envstate {
  Score  score;
  double expected;
};

double d_abs( double );
