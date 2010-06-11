#define true 1;
#define false 0;

// inputtypes
#define NONE 0
#define RMAPPER_CUR 1
#define PROBCALC_CUR 2
#define RMAPPER_V09 3




static void usage(char *progname);
int file_iterator_n(char **paths, int npaths); 
int file_iterator(char *path);
static void usage(char *progname); 
int variant_transform(char *path);
void editstr_to_stats(char * str, long readloc, int is_forward);
int assert_editstring_char(char echar) ;
int complement(char ech);
