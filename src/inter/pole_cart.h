/** \File env_pole_cart.h
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
  double mp;
  double l;
  double g;
  double x,dx,ddx;
  double H,dH,ddH;
  double T;
};

