/** File step.h
*/

//#define  DEBUG            1

#define  ASSERT           1

//#define  INTEGERS_ONLY    1

//#define  USE_MPI          1

#define  GRADIENT         1

#define  NONE            -1

#define  FALSE            0
#define  TRUE             1

#define  MAX_INT_PLUS_ONE  2147483648.0

#ifndef  M_PI
#define  M_PI  3.14159265358979323846264338328
#endif

#define  NUM_CODONS      39

#define  MAX_CODE      4096
#define  MAX_BAR        256
#define  MAX_CELL        64
#define  MAX_PUSH      1024

#define  MAX_LEVEL       24

// (items)
#define  MAX_CALL_STACK  65536

// values for mlevel (mutation level)

#define  COPY          0
#define  GRAD          1
#define  INTERP        2
#define  TUNE          3
#define  TRIM          4
#define  POINT         5
#define  BAR           6
#define  BRANCH        7
#define  CELL          8
#define  BLOCK         9

// values for level.g (transgenic level)

#define  ALIGNED       0
#define  NON_ALIGNED   1
#define  TRANSGENIC    2

// possible values for domain of a channel

#define   UNRESTRICTED   0
#define   FINITE         1
#define   INTEGRAL       2
#define   PRINTABLE      3
#define   ALPHABETIC     4
#define   UPPERCASE      5
#define   LOWERCASE      6

// cost types

#define  LED           5 /* Levenstein Edit Distance     */
#define  KLD2          4 /* KLD scaled between -1 and +1 */
#define  KLD           3 /* Kulback-Leibler Divergence   */
#define  SGN           2 /* L1 but capped at 1,-1        */
#define  SQR           1 /* Sum Squared Error            */
#define  LIN           0 /* L1-Norm (|x-y|)              */

#define  ADD  '+'
#define  AND  '&'
#define  ASN  'a'
#define  BLN  '|'
#define  BRB  ';'
#define  BRF  ':'
#define  CPY  'c'
#define  DEC  'v'
#define  DOT  '.'
#define  EQL  '='
#define  EXP  'e'
#define  GET  '<'
#define  GRT  'g'
#define  INC  '^'
#define  INP  'i'
#define  JSR  'j'
#define  LOD  '{'
#define  LOG  'n'
#define  MOD  '%'
#define  MLT  '*'
#define  NEG  '-'
#define  NOT  '~'
#define  OR   '/'
#define  OUT  'o'
#define  PAY  '$'
#define  PLR  'p'
#define  POP  '!'
#define  PSH  '#'
#define  PUT  '>'
#define  RAN  '?'
#define  RCP  'r'
#define  RND  'z'
#define  ROT  'y'
#define  SCN  's'
#define  SQT  'q'
#define  STO  '}'
#define  SWP  'x'
#define  TNH  'h'
#define  TRG  't'
#define  WRT  'w'

#ifdef SINGLE_FLOAT
typedef float  Floating;
#else
typedef double Floating;
#endif

typedef long long Long;

typedef int Boolean;

typedef unsigned short Index;

typedef struct global_param Global_Param;
typedef struct focus     Focus;
typedef struct level_mg  Level;
typedef struct template_ Template;
typedef struct code      Code;
typedef struct clip      Clip;
typedef struct gradient  Gradient;
typedef struct snode     Snode;
typedef struct agent     Agent;
typedef struct channel   Channel;

struct global_param {
  Boolean   verbose;
  Boolean   multi_process;
  Boolean   terminate_early;
  Boolean   no_share;
  int       base;
  double    tune_magnitude;
  double    interp_rate;
  double    learning_rate;
  long long max_eval;
};

struct level_mg {
  Index m;
  Index g;
};

struct focus {
  Index b;
  Index m;
};

struct template_ {
  char    codon[MAX_CODE];
  Index    ival[MAX_CODE];
  Index     bar[MAX_BAR];
  Index    cell[MAX_CELL];
  Index   param[MAX_PUSH];
  Floating fval[MAX_PUSH];
  Boolean could_align; // whether mutation could be aligned
  Boolean overflow;    // exceeded max codons or bars
  int   stack_size,mem_size;
  int   Lb, Rb;        // Left bar, Right bar
  Index base,num_reg,num_cells;
  Level level;         // mutation,transgenic level
  Index b, c;          // bar, cell
  Index k, k_stop;     // codon, stop_codon
  Index num_focus_bars;
  Focus Fb[16];        // focus bar(s)
  Index Fc,Nc;         // focus cell(s)
  Index d, p, v;       // digit, param, value codon
};

struct code {
  char  *codon;
  Index *ival;
  Index *index;
  Index *bar;
  Index *cell;
  Index *param;
  Floating *fval;
#ifdef GRADIENT
  Floating *delta;
#endif
  int   stack_size,mem_size;
  Index num_reg,num_cells,last_codon;
  Index num_value,num_param;
  Index base;
  Level level;
  Index num_focus_bars;
  Focus Fb[16];      // focus bar(s)
  Index Fc,Nc;       // focus cell(s)
};

struct clip {
  int   Lb, Rb;    // Left bar, Right bar
  Index b, c;      // bar, cell
  Index k, k_stop; // codon, stop_codon
  Index d;         // digit codon
};

struct gradient {
  char op;
  unsigned int arg0_grad;
  unsigned int arg1_grad;
  double value;
  double delta;
};

struct snode {
  Floating local_reg;
  Index c,b,k;
  int r;
#ifdef GRADIENT
  unsigned int local_grad;
#endif
  Index bp_reset;
};

struct agent {
  Code     *cd;
  Floating *reg;
  Floating *stack;
  Floating *mem;
#ifdef GRADIENT
  Gradient *grad;
  double   *delta;
  unsigned int *reg_grad;
  unsigned int *stack_grad;
  unsigned int *mem_grad;
  unsigned int *write_grad;
  unsigned int  grad_index;
  unsigned int  wgi;
#endif
  Snode    *call_stack;
  Boolean  *bstack;
  int stack_items; // number of items on stack
  int sp;          //  stack pointer
  int bp;          // bstack pointer
  int bstack_size;
  int call_stack_size;
  int fp;          // frame pointer;
  Boolean running;
  Long step;
};

struct channel {
  Floating *val;  //    val[i] is the value of item i
  int *index;     //  index[m] is (start) index of message m
  int *length;    // length[m] is length of message m (or -1, if non-message)
  int domain;     // range of acceptable values
  int min_length; // minimum acceptable length for a message
  int max_length; // maximum acceptable length for a message
  int im, max_im; //    input message, max  input message
  int om, max_om; //   output message, max output message
  int si;         //     scan-item
  int wi;         //    write-item
  int max_item;   //  size of array  val[]
  int max_message;//  size of arrays index[] and length[]
};

extern const int power2[17];
extern const char phen[NUM_CODONS];
extern int is_codon[256];
extern int op_uses_digits[256];
extern int is_allowed[256];

extern Global_Param globals;

void init_codons();
void  check_null(void *ptr);
void        fail(char *msg);
int       mysign(double x);

void    init_mutation_state(long seed);
void restore_mutation_state();
void  init_evaluation_state(long seed);

long   my_random();
int    random_bit();
double random_uniform();
double random_gaussian();

int  is_dot_or_digit(int ch);
int   codon_is_digit(int ch);
int codon_from_digit(int digit);
int digit_from_codon(int codon);

void      reset_template(Template *tp);
Code *code_from_template(Template *tp);
Code *empty_code(int num_reg,int num_cells,int stack_size,int mem_size);
void   free_code(Code *cd);

Channel     *new_channel();
Boolean matching_channel(Channel *chan0,Channel *chan1);
void        free_channel(Channel *cn);
Boolean      fetch_input(Channel *in);
Boolean  input_exhausted(Channel *in);
void       clear_message(Channel *chan,int im);
void         cue_message(Channel *chan,int om);
void      output_message(Channel *out);
void  output_non_message(Channel *out);
void       shift_message(Channel *chan,int new_m,int old_m);
int         same_message(Channel *chan0,int m0,Channel *chan1,int m1);
void        copy_message(Channel *copy,Channel *chan,int m);
void accommodate_message(Channel *chan,int max_message);
Boolean       scan_value(Channel *chan,double *p_val);
void         write_value(Channel *chan,double val);

Agent   * new_agent();
void    reset_agent(Agent *agt);
void compress_agent(Agent *agt);
void   next_boolean(Agent *agt);
int    stack_adjust(Agent *agt,int epsilon);
double top_of_stack(Agent *agt);
void     undo_input(Agent *agt,Channel *in);
void    undo_output(Agent *agt,Channel *out);


int    step(Agent *agt,Channel *in,Channel *out);
