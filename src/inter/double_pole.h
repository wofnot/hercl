/** \File env_double_pole.h
*/

typedef struct envstate EnvState;

struct envstate {
  int    running;
  Score  score;
  double xmin;
  double xmax;
  double Hmin;
  double Hmax;
  double Tmax;
  double tau;
  double mc;
  double m1;
  double m2;
  double l1;
  double l2;
  double g;
  double x,dx,ddx;
  double H1,dH1,ddH1;
  double H2,dH2,ddH2;
  double T;
};

