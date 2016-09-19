/** \File human.h
*/

// This environment allows a human to provide the input to an agent at each step,
// facilitating easy interactive testing of an agent.
// Compile in debug mode to use.

typedef struct envstate EnvState;

struct envstate {
  Score  score;
};
