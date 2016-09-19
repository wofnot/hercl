/** \file scan_print.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "step.h"
#include "scan_print.h"

void         scan_cell(Template *tp,Tnode *table,FILE *fp);
void       scan_symbol(char symbol[],int max_length,FILE *fp);
Tnode   *insert_symbol(char symbol[],char code[],Tnode *root);
void       find_symbol(char code[],Tnode *node,char *symbol);
void free_symbol_table(Tnode *root);
void       second_pass(Template *tp,char codon[MAX_CODE],Tnode *table);
void code_from_integer(char code[], int value);
int           scan_int(FILE *fp);
int             mygetc(FILE *fp);
void  print_call_stack(Agent *agt,Snode *sf,FILE *fout);
int        print_value(Floating value,FILE *fout);

const char *freq_phen = "#!cxy-+*rqenahz?%tp<>^v{}j|=g:;&/~iswo$";

/********************************************************//**
   Open the specified file, and check that it exists.
*/
FILE * fopen_check( char *filename, char *mode )
{
  FILE *fp;
  fp = fopen( filename, mode );
  if( fp == NULL ) {
    fprintf(stderr,"Failed to open file: %s\n",filename);
    exit(1);
  }
  return( fp );
}

/********************************************************//**
    Scan the next word on the line into word[].
    Return the first index after the end of the word,
    or -1 if no word was scanned.
*/
int next_word( int j, char line[], char word[] )
{
  int ch;
  int i=0;

  ch = line[j];
  while(  isspace(ch) && ch != '\n' && ch != '\0' ) {
    j++;
    ch = line[j];
  }
  while( !isspace(ch) && ch != '\n' && ch != '\0' ) {
    word[i++] = line[j++];
    ch = line[j];
  }
  word[i] = '\0';

  return(( i > 0 ) ? j : -1 );
}

/********************************************************//**
   Scan code from fp and return pointer to (new) Code.
*/
Code *scan_code( FILE *fp )
{
  Tnode *table=NULL;
  char codon[MAX_CODE];
  char symbol[256];
  char code[256];
  Template temp;
  int ch=' ';
  int j;

  while( ch != '[' && ch != '(' && ch != EOF ) {
    ch = getc( fp );
    if(ch == '`') {
      while( ch != '\n' && ch != EOF ) {
        ch = getc(fp);
      }
    }
  }
  if( ch == EOF ) {
    return( NULL );
  }
  while( ch != '[' ) {
    while( ch != '(' && ch != '[' ) {
      ch = mygetc(fp);
    }
    if( ch == '(' ) {
      // symbols for registers
      scan_symbol( symbol,256,fp );
      ch = getc(fp);
      while( isspace(ch)) {
        ch = getc(fp);
      }
      j=0;
      while( is_dot_or_digit(ch) && j < 255 ) {
        code[j++] = ch;
        ch = getc(fp);
      }
      code[j] = '\0';
      table = insert_symbol(symbol,code,table);
    }
  }
  temp.num_cells  = scan_int( fp );
  temp.num_reg    = scan_int( fp );
  temp.stack_size = scan_int( fp );
  temp.mem_size   = scan_int( fp );
  if( temp.num_cells  < 1 ) temp.num_cells  = 1;
  if( temp.stack_size < 1 ) temp.stack_size = 1;
  if( temp.num_reg    < 0 ) temp.num_reg    = 0;
  if( temp.mem_size   < 0 ) temp.mem_size   = 0;
  temp.level.m = BLOCK; // JUMP;
  while( power2[temp.level.m-BLOCK] <= temp.num_cells ) {
    temp.level.m++;
  }
  temp.level.g = ALIGNED;
  reset_template( &temp );
  temp.Fc = 0;
  temp.Nc = temp.num_cells;
  temp.num_focus_bars = 0;

  while( temp.c < temp.num_cells ) {
    scan_cell( &temp, table, fp );
  }

  if( temp.p > MAX_PUSH ) {
    fail((char *)"Error: Maximum number of tune values exceeded.\n");
  }
  if( temp.b >= MAX_BAR ) {
    fail((char *)"Error: Maximum number of bars exceeded.\n");
  }
  if( temp.k >= MAX_CODE-1 ) {
    fail((char *)"Error: Maximum number of codons exceeded.\n");
  }
  temp.codon[temp.k+1] = '\0';

  strcpy( codon,temp.codon );

  second_pass( &temp, codon, table );
  free_symbol_table( table );

  return( code_from_template(&temp));
}

/********************************************************//**
   Scan the code for a single cell from fp into tp->codon[],
   starting at tp->codon[tp->k+1]. Store new symbols in the
   symbol table, and replace recognized symbols with
   equivalent digits. Check max number of bars not exceeded.
   Also update k, n, c, bar[], cell[], param[], num_param.
*/
void scan_cell(
               Template *tp,
               Tnode *table,
               FILE  *fp
              )
{
  char symbol[256];
  char code[1024];
  int ch=' ';
  int j;

  while( ch != '[' && ch != '(' ) {
    ch = mygetc(fp);
  }
  if( ch == '(' ) { // symbol for this cell
    j=0;
    scan_symbol( symbol,256,fp );
    code_from_integer(code,tp->c);
    table = insert_symbol(symbol,code,table);
    while( ch != '[' ) {
      ch = mygetc(fp);
    }
  }
  while( ch != ']') {
    while(ch != '(' && ch != ']' && !is_codon[ch]) {
      ch = mygetc( fp );
    }
    if( ch == '(' ) { // embed symbol in codon[] buffer
      j = tp->k+1;
      transcribe_codon(tp,'(');
      while((ch = getc(fp)) != ')' && ch != EOF ) {
        transcribe_codon(tp,ch);
      }
      transcribe_codon(tp,'\0');
      find_symbol(code,table,&tp->codon[j+1]);
      if( code[0] == '\0' ) {   // forward reference
        tp->codon[tp->k] = ')'; // leave symbol for 2nd pass
      }
      else {
        tp->k = j-1;            // replace symbol with digits
        while( !is_codon[ch] ) {
          ch = mygetc(fp);
        }
        j=0;
        while( code[j] != '\0' ) {
          transcribe_codon(tp,code[j++]);
        }
        switch(ch) {
         case GET: case PUT: case LOD: case STO:
         case INC: case DEC: case EQL: case GRT:
         case JSR: case PSH: case DOT:
           break;
         default:
           fprintf(stderr,"Missing operator:%2d[",tp->c);
           fail((char *)&tp->codon[tp->bar[tp->cell[tp->c]]]);
           break;
        }
      }
    }
    else if( is_codon[ch] ) {
      tp->b += ( ch == BLN );
      //tp->p += ( ch == PSH ); ?????????????
      transcribe_codon(tp,ch);
      ch = mygetc( fp );
    }
  }
  transcribe_codon(tp,']');// will change to '|' in second_pass()
  next_bar( tp );
  tp->c++;
}

/********************************************************//**
   Scan a symbol from the specified FILE pointer.
*/
void scan_symbol( char symbol[],int max_length, FILE *fp )
{
  int ch;
  int j=0;
  while(( ch = getc(fp)) != ')' && j < max_length ) {
    symbol[j++] = ch;
  }
  if( j == max_length ) {
    symbol[j-1] = '\0';
    printf("Symbol too long: %s\n",symbol );
    exit( 1 );
  }
  symbol[j] = '\0';
}

/********************************************************//**
   Insert symbol into (BST) symbol table.
*/
Tnode *insert_symbol(
                     char symbol[],
                     char code[],
                     Tnode *root
                    )
{
  Tnode *parent,*child,*node;
  int test;
  node = (Tnode *)malloc(sizeof(Tnode));
  check_null(node);
  node->left  = NULL;
  node->right = NULL;
  node->symbol=(char *)malloc((strlen(symbol)+1)*sizeof(char));
  check_null(node->symbol);
  strcpy(node->symbol,symbol);
  node->code = (char *)malloc((strlen( code )+1)*sizeof(char));
  check_null(node->code);
  strcpy(node->code , code );
  if( root == NULL ) {
    root = node;
  }
  else {
    child = root;
    while( child != NULL ) {
      parent = child;
      test = strcmp(symbol,parent->symbol);
      if( test < 0 ) {
        child = parent->left;
      }
      else if( test > 0 ) {
        child = parent->right;
      }
      else {
        fprintf(stderr,"Duplicate symbol: ");
        fail((char *) symbol );
      }
    }
    if( test < 0 ) {
      parent->left  = node;
    }
    else {
      parent->right = node;
    }
  }
  return( root );
}

/********************************************************//**
   If symbol is a character, return ascii value;
   otherwise, look for symbol in symbol table.
   If not found, return -1.
*/
void find_symbol(
                 char  *code,
                 Tnode *node,
                 char  *symbol
                )
{
  int test;
  if( symbol[0] == '\'' && isprint(symbol[1])) {
    code_from_integer(code,symbol[1]); // ascii character
  }
  else {
    while(   node != NULL
          &&(test = strcmp(symbol,node->symbol)) != 0) {
      if( test < 0 ) {
        node = node->left;
      }
      else {
        node = node->right;
      }
    }
    if( node == NULL ) {
      code[0] = '\0';
    }
    else {
      strcpy( code, node->code );
    }
  }
}

/********************************************************//**
   Free the space occupied by symbol table (BST).
*/
void free_symbol_table( Tnode *node )
{
  if( node != NULL ) {
    free_symbol_table( node->left );
    free_symbol_table( node->right);
    free( node->symbol );
    free( node->code );
    free( node );
  }
}

/********************************************************//**
   Resolve any forward-referenced symbols in the code.
   Also compute ival, fval.
*/
void second_pass(
                 Template *tp,
                 char   codon[MAX_CODE],
                 Tnode *table
                )
{
  Fraction frac;
  char symbol[256];
  char code[1024];
  int k=1;
  int j;

  reset_template( tp );

  while( codon[k] != '\0') {
    tp->d = tp->k+1;
    if( codon[k] == '(') { // forward referenced symbol
      j=0;
      k++; // skip '('
      while( codon[k] != ')') {
        symbol[j++] = codon[k++];
      }
      symbol[j] = '\0';
      find_symbol( code,table,symbol );
      if( code[0] == '\0' ) {
        fprintf(stderr,"Undefined symbol: ");
        fail((char *) symbol );
      }
      k++;         // skip ')'
      j=0;
      while( code[j] != '\0' ) {
        transcribe_codon( tp,code[j++]) ;
      }
      switch( codon[k] ) { // operator or dot
       case GET: case PUT: case LOD: case STO:
       case INC: case DEC: case EQL: case GRT:
       case JSR: case PSH: case DOT:
         break;
       default:
         fail((char *)"Missing operator");
         break;
      }
    }
    while( is_dot_or_digit(codon[k])) {
      transcribe_codon(tp,codon[k++]);
    }
    transcribe_codon(tp,codon[k++]); // operator

    if( tp->k > tp->d ) {
      switch( tp->codon[tp->k] ) { // operator
       case GET: case PUT: case LOD: case STO:
       case INC: case DEC: case EQL: case GRT:
       case BRB: case BRF: case JSR:
         tp->ival[tp->d+1]=int_from_codons(tp->codon,tp->d);
         break;

       case PSH:
         if( codon[k] == NEG ) {
           transcribe_codon(tp,codon[k++]); // #-
         }
         tp->ival[tp->d+1] = tp->v;
         scan_fraction(&frac,tp->codon,tp->d);
         tp->fval[tp->v] = val_from_fraction(&frac);
         next_value( tp );
         if( frac.has_dot ) {
           tp->param[tp->p] = tp->d;
           next_param( tp );
         }
         break;

       case BLN: case ']': // no action
         break;

       default:
         fail((char *)"Missing operator");
         break;
      }
    }
    tp->ival[tp->d] = tp->k - tp->d;

    if( tp->codon[tp->k] == ']') {
      tp->codon[tp->k]  = BLN;
      next_bar( tp );
      tp->bar[tp->b]    = tp->k;
      tp->cell[++tp->c] = tp->b;
    }
    else if( tp->codon[tp->k] == BLN ) {
      next_bar( tp );
      tp->bar[tp->b] = tp->k;
    }
  }
  tp->codon[tp->k+1] = '\0';
}

/********************************************************//**
   store to code[] the digits encoding value, in globals.base.
*/
void code_from_integer( char code[], int value )
{
  int size=1,digit;
  int j=0;
  while( size <= value / globals.base ) {
    size *= globals.base;
  }
  while( size > 0 ) {
    digit = value / size;
    code[j++] = codon_from_digit( digit );
    value -= digit * size;
    size /= globals.base;
  }
  code[j] = '\0';
}

/********************************************************//**
   Scan an integer from the specified FILE pointer.
*/
int scan_int( FILE *fp )
{
  int ch;
  int n;

  ch = mygetc(fp);
  while( ch < '0' || ch > '9' ) {
    ch = mygetc(fp);
  }
  ungetc(ch,fp);
  fscanf(fp,"%d",&n);

  return(n);
}

/********************************************************//**
   Get the next character, skipping comments if necessary.
*/
int mygetc( FILE *fp )
{
  int ch;
  ch = getc(fp);
  if(ch == '`') {
    while( ch != '\n' && ch != EOF ) {
      ch = getc(fp);
    }
  }
  if( ch == EOF ) {
    fail((char *)"File ended unexpectedly\n");
  }
  return( ch );
}

/********************************************************//**
   Print one (specified) cell from code.
*/
void print_cell( Code *cd, int c, FILE *fout )
{
  int n,k;
  fprintf(fout,"%2d[",c);
  for( n = cd->cell[c]; n < cd->cell[c+1]; n++ ) {
    for( k = cd->bar[n]+1; k < cd->bar[n+1]; k++ ) {
      putc(cd->codon[k],fout);
    }
    if( cd->codon[k] != BLN ) {
      fail((char *)"WHOOPS! (bar)");
    }
    if( n < cd->cell[c+1]-1 ) {
      putc('|',fout);
    }
  }
  fprintf(fout,"]\n");
}

/********************************************************//**
   Print code to specified FILE pointer.
*/
void print_code( Code *cd, FILE *fout )
{
  int c;

  if( cd != NULL ) {
    fprintf(fout,"[c=%d,r=%d,s=%d,m=%d]\n",cd->num_cells,
            cd->num_reg,cd->stack_size,cd->mem_size);
    for( c=0; c < cd->num_cells; c++ ) {
      print_cell( cd, c, fout );
    }
  }
}

/********************************************************//**
*/
char * sprint_code( Code *cd, char *ts )
{
  return( ts + sprintf( ts," %s\n",cd->codon ));
}

/********************************************************//**
   Print the state of the agent.
*/
void print_state(Agent *agt,Channel *in,Channel *out,FILE *fout)
{
  Code  *cd = agt->cd;
  Snode *sf;
  int margin=3;
#ifdef DEBUG
  // char keych;
#endif
  int fp;
  int b,k;
  for( fp=0; fp <= agt->fp; fp++ ) {
    print_call_stack( agt, agt->call_stack + fp, fout );
  }
  sf = agt->call_stack + agt->fp;

  fprintf(fout,"\n[b=");
  for( b = sf->bp_reset; b <= agt->bp; b++ ) {
    //fprintf(fout, agt->bstack[agt->bp] ? "T" : "F" );
    fprintf(fout, agt->bstack[b] ? "T" : "F" );
  }
  fprintf(fout,"]\n");
  fprintf(fout,"i: ");
  if( in != NULL ) {
    for( k = in->si; k < in->index[in->im+1]; k++ ) {
      print_value( in->val[k], fout );
    }
  }
  fprintf(fout,"\n");
  fprintf(fout,"o: ");
  if( out != NULL ) {
    for( k = out->index[out->om+1]; k < out->wi; k++ ) {
      print_value( out->val[k], fout );
    }
  }
  fprintf(fout,"\n");
  fprintf(fout,"m: ");
  for( k=0; k < cd->mem_size; k++ ) {
    print_value( agt->mem[k], fout );
    if( k % 64 == 63 ) {
      fprintf(fout,"\n   ");
    }
  }
  fprintf(fout,"\nr: ");
  for( k=0; k < cd->num_reg; k++ ) {
    print_value( agt->reg[k], fout );
  }
  putc(' ', fout );
  print_value( sf->local_reg, fout );
  fprintf(fout,"\ns: ");
  for( k=0; k < cd->stack_size; k++ ) {
    int width = print_value(agt->stack[k], fout );
    if( k < agt->sp ) {
      margin += width;
    }
  }
  fprintf(fout,"\n");
  for( k=0; k < margin; k++ ) {
    putc(' ', fout );
  }
  fprintf(fout,"^\n");

  if( sf != NULL ) {
    fprintf(fout," %2d[",sf->c);
    for( k = cd->bar[cd->cell[sf->c]]+1;
         k < cd->bar[cd->cell[sf->c+1]]; k++ ) {
      putc(cd->codon[k], fout );
    }
    fprintf(fout,"]\n   ");
    for( k = cd->bar[cd->cell[sf->c]]; k < sf->k; k++ ) {
      putc(' ',fout);
    }
    fprintf(fout,"^\n");
  }
#ifdef DEBUG
  printf("? ");
  //scanf("%c", &keych );
#endif
}

/********************************************************//**
   Print the contents of the call stack.
*/
void print_call_stack( Agent *agt, Snode *sf, FILE *fout )
{
  if( sf != NULL ) {
    fprintf(fout,"[c=%d,i=%d,r=",sf->c,
            sf->k - agt->cd->bar[agt->cd->cell[sf->c]]);
    if( sf->r < 0 ) {
      fprintf(fout,".,r.=");
    }
    else {
      fprintf(fout,"%d,r.=",sf->r );
    }
    print_value( sf->local_reg, fout );
    fprintf(fout,"]");
  }
}

/********************************************************//**
   Print value, either as a character or integer.
   Return the number of characters printed.
*/
int print_value( Floating value, FILE *fout )
{
  int ch = (int) round( value );
  int margin=1;
  if(   value <   MAX_INT_PLUS_ONE
     && value >= -MAX_INT_PLUS_ONE && ch == value ) {
    if( ch > 0 && ch < 256 && isprint( ch )) {
      putc( ch, fout );
    }
    else if( ch == 0 ) {
      putc('.', fout );
    }
    else {
      fprintf(fout,"[%d]",ch);
      margin = 2;
      while( ch > 0 ) {
        margin++;
        ch /= 10;
      }
    }
  }
  else {
    fprintf(fout,"[%1.2f]",value);
    margin = 6;
    if( value >= 10.0 ) {
      margin +=   (int)floor(log( value)/log(10.0));
    }
    else if( value <= -10.0 ) {
      margin += 1+(int)floor(log(-value)/log(10.0));
    }
    else if( value < 0.0 ) {
      margin += 1;
    }
  }
  return( margin );
}

/********************************************************//**
   Print message m from specified channel.
*/
void print_message( Channel *chan, int m, FILE *fout )
{
  double value;
  int preceded_by_char = TRUE;
  int length;
  int ch;
  int i0,i;
  if( chan != NULL && m <= chan->om+1 ) {
    i0 = chan->index[m];
    if( m == chan->om+1 ) { // partial message not yet output
      length = chan->wi - i0;
    }
    else {
      length = chan->length[m];
    }
    if( length == -1 ) {
      fprintf(fout,"\\\n");
    }
    else if( length == 0 ) {
      fprintf(fout,",\n");
    }
    else {
      for( i=i0; i < i0+length; i++ ) {
        value = chan->val[i];
        ch = (int)(floor(0.5+value));
        if(   value <   2147483648.0
           && value >= -2147483648.0 && ch == value ) {
          if( ch >  0 && ch < 256 && isprint( ch )) {
            if( !preceded_by_char ) {
              putc(' ',fout);
            }
            preceded_by_char = TRUE;
            if( ch == ',' ) {
              fprintf(fout,"\\,");
            }
            else if( ch == '\\' ) {
              fprintf(fout,"\\\\");
            }
            else {
              putc( ch, fout );
            }
          }
          else { 
            fprintf(fout,",%d",ch);
            preceded_by_char = FALSE;
          }
        }
        else {
          fprintf(fout,",%1.2f",(double)value);
          preceded_by_char = FALSE;
        }
      }
      fprintf(fout,"\n");
    }
  }
}

/********************************************************//**
   Increment output message and scan items from line[]
   to form the content of the new message.
   Set length of new message and index of next message.
   Update write index to be the first index of the next message.
*/
void scan_next_input( Channel *in, char line[] )
{
  double value;
  int j=0;

  in->om++;
  accommodate_message(in,in->om);
  in->wi = in->index[in->om];

  if( line[0] == '\\' && line[1] == '\0' ) {
    in->length[in->om] = -1; // input fails
  }
  else if( line[0] == ',' && line[1] == '\0' ) {
    in->length[in->om] = 0;  // succeeds, but zero input items
  }
  else {
    while( line[j] != '\0' ) { // && in->wi < in->max_i ) {
      if( line[j] == ',' && line[j+1] == '?' ) {
        write_value(in,0.0);
        j = j+2;
        if( line[j] == ' ' ) {
          j++;
        }
      }
      else if(line[j] == ',' &&(line[j+1]=='-' || isdigit(line[j+1]))){
        int n;
        j++;
        sscanf(line+j,"%lf%n",&value,&n);
        write_value(in,value);
        j += n;
        //printf("%1.2f %d\n",value,n);
        /*if( line[j] == '-' ) {
          j++;
        }
        while( line[j] == '.' || isdigit(line[j])) {
          j++;
        }*/
        if( line[j] == ' ' ) {
          j++;
        }
      }
      else {
        if( line[j] == '\\' && line[j+1] != '\0' ) {
          j++;
        }
        write_value(in,(Floating)line[j++]);
      }
    }
    in->length[in->om] = in->wi - in->index[in->om];
  }
  in->index[in->om+1] = in->wi;
}

/********************************************************//**
   Scan several inputs from fp, until EOF, and store them
   successively into the specified buffer.
*/
void scan_inputs( Channel *in, FILE *fp )
{
  char line[MAX_LINE];
  int ch;
  while(  (in->max_om < 0 || in->max_om > in->om )
        && fgets(line,MAX_LINE,fp) > 0 ) {
    ch = getc( fp );
    if( ch != '\n' ) {
      ungetc( ch, fp );
      line[strlen(line)-1] = '\0';
    }
    scan_next_input( in, line );
  }
}

/********************************************************//**
   Scan code from string s[] and return pointer to (new) Code.
*/
Code *code_from_string( char s[] )
{
  Template temp;
  int k;
  temp.num_cells  = int_from_codons(s,0);
  k = skip_int(s);
  temp.num_reg    = int_from_codons(s,k);
  k = k + skip_int(s+k);
  temp.stack_size = int_from_codons(s,k);
  k = k + skip_int(s+k);
  temp.mem_size   = int_from_codons(s,k);
  k = k + skip_int(s+k);

  temp.level.m = BLOCK; // JUMP;
  while( power2[temp.level.m-BLOCK] < temp.num_cells ) {
    temp.level.m++;
  }
  temp.level.g = ALIGNED;
  reset_template( &temp );
  temp.Fc = 0;
  temp.Nc = temp.num_cells;
  temp.num_focus_bars = 0;

  second_pass( &temp,s+k,NULL );

  return( code_from_template(&temp));
}

/********************************************************//**
   Skip the next integer in string and return following index
*/
int skip_int( char s[] )
{
  int k=0;
  while( isdigit(s[k]) && s[k] != '\0' ) {
    k++;
  }
  while( isspace(s[k]) && s[k] != '\0' ) {
    k++;
  }
  return( k );
}

/********************************************************//**
*/
int string_from_int( char s[], int n )
{
  int digit;
  int size=1;
  int k=0;
  if( n == 0 ) {
    s[k++] = '0';
  }
  else {
    while( size <= n / 10 ) {
      size *= 10;
    }
    while( size > 0 ) {
      digit = n / size;
      s[k++] = '0'+ digit;
      n -= digit * size;
      size /= 10;
    }
  }
  s[k++] = ' ';
  return( k );
}

/********************************************************//**
   Convert code into a string s[] and return length of string.
*/
int string_from_code( char s[], Code *cd )
{
  int c,k;
  k =     string_from_int(s,  cd->num_cells);
  k = k + string_from_int(s+k,cd->num_reg);
  k = k + string_from_int(s+k,cd->stack_size);
  k = k + string_from_int(s+k,cd->mem_size);
  strcpy( s+k,cd->codon );
  s[k] = '[';
  for( c=1; c <= cd->num_cells; c++ ) {
    s[k+cd->bar[cd->cell[c]]] = ']';
  }
  k = k + cd->last_codon + 1;
  s[k] = '\0';
  return( k );
}

/********************************************************//**
   Increment tp->n.
   If MAX_BAR is exceeded, set tp->overflow to TRUE.
*/
void next_bar( Template *tp )
{
  if( tp->b < MAX_BAR-1 ) {
    tp->b++;
  }
  else {
    tp->overflow = TRUE;
  }
}

/********************************************************//**
   Increment tp->p.
   If MAX_PUSH is exceeded, set tp->overflow to TRUE.
*/
void next_param( Template *tp )
{
  if( tp->p < MAX_PUSH ) {
    tp->p++;
  }
  else {
    tp->overflow = TRUE;
  }
}

/********************************************************//**
   Increment tp->p.
   If MAX_PUSH is exceeded, set tp->overflow to TRUE.
*/
void next_value( Template *tp )
{
  if( tp->v < MAX_PUSH ) {
    tp->v++;
  }
  else {
    tp->overflow = TRUE;
  }
}

/********************************************************//**
   Increment tp->k and set codon[k] to the specified character.
   If MAX_CODE is exceeded, set tp->overflow to TRUE.
*/
void transcribe_codon( Template *tp, int ch )
{
  if( tp->k < MAX_CODE-1 ) {
    tp->k++;
  }
  else {
    tp->overflow = TRUE;
  }
  tp->codon[tp->k] = ch;
  //if( !is_dot_or_digit(tp->codon[tp->k-1])) {
  //  tp->ival[tp->k] = 0; // single-codon instruction
  //}
}

/********************************************************//**
   Scan fraction from string codon[], starting at codon[d],
   and store value of the fration into (*fr)
*/
void scan_fraction(
                   Fraction *fr,
                   char     *codon,
                   int       d
                  )
{
  fr->whole = 0;    // value is whole+(num/denom)
  fr->num   = 0;    // numerator
  fr->denom = 1;    // denominator
  fr->is_negative = FALSE;
  fr->has_dot     = FALSE;

  // scan whole number part, numerator and denominator
  while( codon_is_digit(codon[d])) {
    fr->whole = fr->whole*globals.base + digit_from_codon(codon[d]);
    d++;
  }
  if( codon[d] == DOT ){
    fr->has_dot = TRUE;
    d++; // skip the DOT
    if( codon_is_digit(codon[d])) {
      // scan fractional part ( num / denom )
      while( codon_is_digit(codon[d])) {
        fr->num = fr->num*globals.base + digit_from_codon(codon[d++]);
        fr->denom *= globals.base;
      }
    }
  }
  if( codon[d+1] == NEG ) { // value is negated
    fr->is_negative = TRUE;
  }
  //printf("%ld+%ld/%ld ->",whole,num,denom);
}

/********************************************************//**
   Return floating-point value of specified fraction
*/
Floating val_from_fraction( Fraction *fr )
{
  Floating val = fr->whole + fr->num/(double)fr->denom;
  if( fr->is_negative ) {
    val = -val;
  }
  return( val );
}

/********************************************************//**
   Convert sequence of digits into an integer, and return it.
*/
int int_from_codons( char *codon, int d )
{
  int value = 0;
  if( d >= 0 ) {
    while( codon_is_digit( codon[d] )) {
      // && value < ( MAX_INT/globals.base )) {
      value = globals.base*value + digit_from_codon(codon[d++]);
    }
  }
  return( value );
}

/********************************************************//**
   Print the parameters of the channel
*/
void print_channel( Channel *chan, FILE *fout )
{
  int i,m;
  printf("[");
  for( m=0; m <= chan->om; m++ ) {
    if( chan->length[m] > 0 ) {
      printf("%1.2f",chan->val[chan->index[m]]);
      for( i=1; i < chan->length[m]; i++ ) {
        printf(",%1.2f",chan->val[chan->index[m]+i]);
      }
    }
    printf("|");
  }
  printf("\n");
  printf("im=%d om=%d si=%d wi=%d\n",chan->im,chan->om,chan->si,chan->wi);
}

/********************************************************//**
   Add frequencies for this code to freq[][].
*/
void aggregate_freq(
                    int freq[NUM_CODONS][NUM_CODONS],
                    Code *cd
                   )
{
  int k,p,q;
  int op=BLN;
  for( p=0; p < NUM_CODONS && freq_phen[p] != op; p++ )
    ;
  for( k=1; k < cd->last_codon; k++ ) {
    op = cd->codon[k];
    for( q=0; q < NUM_CODONS && freq_phen[q] != op; q++ )
      ;
    if( q < NUM_CODONS ) {
      freq[p][q]++;
      p = q;
    }
  }
}

/********************************************************//**
   Print frequency table.
*/
void print_freq(
                int freq[NUM_CODONS][NUM_CODONS],
                FILE *fout
               )
{
  int ch;
  int p,q,n,sum;
  fprintf(fout,"  ");
  for( q=0; q < NUM_CODONS; q++ ) {
    putc(freq_phen[q],fout);
  }
  putc('\n',fout);
  for( p=0; p < NUM_CODONS; p++ ) {
    sum = 0;
    fprintf(fout,"%c ",freq_phen[p]);
    for( q=0; q < NUM_CODONS; q++ ) {
      n = freq[p][q];
      if( n == 0 ) {
        ch = '.';
      }
      else if( n < 10 ) {
        ch = '0' + n;
      }
      else if( n < 36) {
        ch = 'a' + n - 10;
      }
      else if( n < 62 ) {
        ch = 'A' + n - 36;
      }
      else {
        ch = '!';
      }
      putc(ch,fout);
      sum += n;
    }
    fprintf(fout,"%3d\n",sum);
  }
}




