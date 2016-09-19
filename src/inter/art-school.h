/** \File art-school.h

art-school implements the interactive environment for training artists,
within the context of artist-critic co-evolution (see artist-critic/README.txt).
The environment is not really easy to use standalone, it is intended to be used as a component of the larger system.

The interface for HERCL agents to draw pictures is based on the classic turtle-graphics environment. Via numerical messages the agent can control a virtual pen, adjusting its position, size, colour, and whether it is touching the canvas. Each message output by the agent is interpreted as a single command, based on the first item in the message. If the command requires arguments these are taken from the subsequent items in the message.

The agent has at most seven commands to choose from. Expected arguments, with default values, shown in brackets.

  TOGGLE      ()           -- toggles the state of the pen (on/off the canvas)
  FORWARD     (x = 10)     -- moves forward x pixels
  TURN        (x = 10)     -- turns x degrees clockwise
  SET_SIZE    (x =  1)     -- sets pen size to x. (Argument is taken modulo 4).
  SET_COLOUR  (l, [h, s])  -- sets pen colour (lightness, hue, saturation).
                              If colour mode is disabled the hue and saturation arguments are ignored
                              If less than 3 arguments are provided the remaining channels are left unchanged.
                              
[ (these are only permitted if pushpop is enabled in the settings)
  PUSH_TURTLE ()           -- Pushes the location, rotation, colour and size of the pen (turtle) to a stack.
  POP_TURTLE  ()           -- Pops pen features from the stack.
]
  
There are several parameters which determine properties of the drawing environment as well as the evaluation of the artist. The art-school executable expects a file art-school.txt to be present in the working directory, specifying each of these parameters. Information on these parameters and the file format is available in artist-critic/art-school.template. However, if you are using art-school as a component of artist-critic evaluation, art-school.txt will be automatically generated and you needn't worry about the file format.
*/


#include "OCVCanvas.h"
#include "../hashmap/hashmap.h"

#define SETTINGS_FILE "./art-school.txt"

#define STACK_SIZE 100

#define MAX_CRITICS 1000

typedef enum _command {
  TOGGLE, // []
  FORWARD, // [distance]
  TURN, // [delta_angle]
  SET_SIZE, // [size]
  SET_COLOUR, // [colour]
  PUSH_TURTLE,
  POP_TURTLE,
  CMD_COUNT
} Command;

#define MAX_HISTORY 10000
#define NUM_KEYS 10000

typedef struct _turtleFrame {
  double penX;
  double penY;
  double angle;
  colour_t penColour;
  int penSize;
} turtleFrame;

typedef struct _envstate {
// SETTINGS {
  // char critic_fn[100];
  char output_fn[100]; // empty for no save
  double tolerance;        // percent
  int pushpop;
  int no_loop;
  int max_step;
  double max_movement; // negative for no limit
  int use_colour;
  
// }
  
  Score  score;
  Canvas canvas;
  int step;
  
  double angle; // radians
  
  Agent  *critics[MAX_CRITICS]; // null terminated
  double weights[MAX_CRITICS];
  
  Channel *chan;  // Separate channel between critic and environment
  
  turtleFrame turtleStack[STACK_SIZE];
  int sp;
  
  char history[MAX_HISTORY];
  
  char keys[NUM_KEYS][MAX_HISTORY];
  double values[NUM_KEYS];
  int key_index;
  
  map_t cost_map;
} EnvState;


