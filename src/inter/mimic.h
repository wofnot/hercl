/** \File mimic.h
*/

typedef struct _envstate {
  Score  score;
  
  double target;
  
  Agent  *agent;  // Agent-to-mimic
  Channel *chan;  // Separate channel between agent-to-mimic and environment
  
} EnvState;
