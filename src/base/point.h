/** \File point.h
*/

#define  WHOLE   0
#define  FRONT   1
#define  MIDDLE  2
#define  BACK    3
#define  FRINGE  4

#ifdef DEBUG
extern char ds[64*MAX_CODE];
extern int  di;
extern int  print_debug;

void         print_spaces(int num_spaces);
void           print_line();
void    print_parent_cell(Code *cd,int c);
void           print_clip(Code *cd,Clip *cp);
int   print_clip_template(Code *cd,Clip *cp,Template *tp);
void       print_template(Template *tp);
void print_completed_cell(Template *tp,int c);
void          flush_debug();
#endif

int      random_geometric(double mean,int max);
int     choose_jump_index(int c,int num_cells);
void      mutate_fraction(Fraction *fr,int mlevel,double scale, int singular);

void      reset_cell_bars(Template *tp,Code *cd,Clip *cp);
void       insert_integer(Template *tp,int n);
void           insert_dot(Template *tp);
void      insert_operator(Template *tp,int op);
void          insert_push(Template *tp,Fraction frac);
void    insert_comparison(Template *tp,int op);
void   insert_instruction(Template *tp);
void interpolate_fraction(Template *tp,Code *cd0,Code *cd1,int p,double alpha);
void       fix_bar_codons(Template *tp,char d_codon,int  must_change);
void           fix_digits(Template *tp,Code *cd,Clip *cp,int must_change);

/*

double     random_cauchy( int min );
//void   focus_weights( int w[],int c0,int i,int num_cells );

int      random_register( Template *tp );
//void     insert_register( Template *tp, int r );









void       mutate_digits( Template *tp, Code *cd, Clip *cp );
#ifdef DEBUG
void        print_spaces( int num_spaces );
void          print_line();
int  print_clip_template( Code *cd, Clip *cp, Template *tp );
void          print_clip( Code *cd, Clip *cp );
void      print_template( Template *tp );
void    print_parent_cell( Code *cd, int c );
void print_completed_cell(Template *tp,int c );
void         flush_debug();
#endif
*/
