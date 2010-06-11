// $Id: probcalc_mp.h 321 2009-02-15 23:53:48Z dalcaadr $

// binary
#define true 1
#define false 0

#define MAX_READS 100 		// default for maximum read mappings to read
#define HIST_BINS 2000		// bins in the histogram

#define DEBUGLIMIT 1000000	// default matpair read limit for debug mode

#define MEAN_MPS 50000		// default for number of USEFUL mp to include in mean comp

#define PCHANCE_CUTOFF 0.05 // if 1.00, no pchance cutoff since all pchances will be less than or equal to 1
#define PGENOME_CUTOFF 0.95 // if 0.00, no pgenome cutoff since all pgenomes will be greater than or equal to 0

#define PRINT_MAX 20		// default for the number of matepair mappings to print

#define ALMOST_ZERO	0.000000001
#define ALMOST_ONE	0.999999999

// some math macros...
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define ABS(a) (MAX(a, -a))
#define FABS(a) (MAX(a, -a))


// pass types
#define MEAN_PASS 0
#define OUTPUT_PASS 1
#define COUNT_PASS 2 

// sorting
#define SORT_PGENOME 0
#define SORT_PCHANCE 1
#define SORT_NORMODDS 2

// time 
#define SYS_TIME 0
#define STATS_TIME 1
#define ANALYSIS_TIME 2
#define COUNT_TIME 3
#define TOTAL_TIME 4 // always the largest value!

// input file types 
#define ASCII 0
#define BINARY 1

#define MAXLINELEN 1023

// mate pair mapping struct
typedef struct {
	mapping_t 	*fwd_rs;
	mapping_t 	*rev_rs;
	uint64_t	dist;
	double		pchance;
	double		pgenome;
	double     	normodds;
} mate_pair_val ;



/*
 * usage of the program.
 */
static void usage(char *progname);

/*
 * The Mate Pair analysis. This does one or two passes over the combinations
 * of mate pairs...
 */
void mp_analysis(mapping_t *fwd_maps, mapping_t *rev_maps, int fwd_nr, 
		int rev_nr, int pass_type);

/*
 * test if the string readname is a forward or a reverse read
 */ 
int is_forward_test(char * readname, int readlen);

/*
 * compare mate pairs
 */
static int mate_pair_val_cmp(const void *a, const void *b);

/*
 * Pass through the readfile and mappingfile once through the method
 * specified by pass_type (MEAN_PASS or OUTPUT_PASS)
 */
uint64_t filepass(char * mappingfilename, int pass_type);

/*
 * Compute the cumsum histogram
 */
void compute_cumsum(); 

/*
 * Print histograms to log files.
 */
void print_hists();

/*
 * if the mate pair is "good" (see definition of good above), 
 * the method returns the distance between the two, including the reads
 * themselves. 
 * 
 * otherwise, returns 0
 */
uint64_t good_mp_dst(mapping_t * fwd_map, mapping_t * rev_map);

/* 
 * increment the mean, standard deviation, and histogram, using the
 * given distance
 */
void increments_stats(uint64_t good_mps_dist);

/*
 * Compute the probabilities, given a mapping. Add it to the mate pair array
 * if the probabilities pass the thresholds.
 */
inline int add_p_stats(mapping_t * fwd_map, mapping_t * rev_map, 
		mate_pair_val *mp_set, double * totnormodds, int * mp_set_index );

/*
 * Source of fuction: SHRiMP (original code probably by Steve Rumble)
 * 
 * Return a string on the stack corresponding to an unsigned integer that also
 * features commas. E.g.: int 1000 yields "1,000".
 *
 * There is a pool of sequentially allocated buffers returned, so this should
 * be safe to use multiple times in function arguments.
 */
char * comma_integer(uint64_t val);

/*
 * Print a dashed line of 80 chars, including the newline char at the end
 */
void print_dashed_line();

/*
 * Read a line from the mapping file, and assign the apropriate values to 
 * the mapping_t mapping.
 * The mappingfile can be binary or ascii, with input_file_type appropriately 
 * set.  
 */
int read_probcalc_line(FILE *mappingfile, mapping_t * mapping);
