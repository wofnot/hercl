/** mystring.h */

#define  TRUE         1
#define  FALSE        0

#define  NUM_TRAIN_TEST  1000

int   random_geometric(double mean,int max);
int         int_cauchy(int min,int max);
double   signed_cauchy();
double   pseudo_cauchy(double median);
double random_gaussian();
int   random_character();
void      print_string(FILE *fp,int length,int s[]);
void   print_character(FILE *fp,int ch);
void       divide_data(int train_items,int test_items,
                       int item_in[],int item_out[]);
