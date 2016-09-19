/** \File scan_print.h
*/

#define  MAX_INPUT_BUF  65536
#define  MAX_INPUTS      1024
#define  MAX_LINE        4096

typedef struct fraction Fraction;
typedef struct tnode Tnode;

struct fraction {
  Long whole;   // value is whole+(num/denom)
  Long num;
  Long denom;
  int  is_negative;
  int  has_dot;
};

struct tnode {
  char  *symbol;
  char  *code;
  Tnode *left;
  Tnode *right;
};


FILE * fopen_check(char *filename,char *mode);

int  next_word(int j,char line[],char word[]);

Code   * scan_code(FILE *fp);
void    print_cell(Code *cd,int c,FILE *fp);
void    print_code(Code *cd,FILE *fp);
char * sprint_code(Code *cd,char *ts);

void   print_state(Agent *agt,Channel *in,Channel *out,FILE *fp);
void print_message(Channel *chan,int r,FILE *fp);

void scan_next_input(Channel *in,char line[]);
void     scan_inputs(Channel *in,FILE *fp);

Code *code_from_string(char s[]);
int           skip_int(char s[]);
int    string_from_int(char s[],int n);
int   string_from_code(char s[],Code *cd);

void         next_bar(Template *tp);
void       next_param(Template *tp);
void       next_value(Template *tp);
void transcribe_codon(Template *tp,int ch);

void         scan_fraction(Fraction *fr,char *codon,int d);
Floating val_from_fraction(Fraction *fr);
int        int_from_codons(char *codon, int d);

void aggregate_freq(int freq[NUM_CODONS][NUM_CODONS],Code *cd);
void     print_freq(int freq[NUM_CODONS][NUM_CODONS],FILE *fp);

void print_channel(Channel *in,FILE *fout);

