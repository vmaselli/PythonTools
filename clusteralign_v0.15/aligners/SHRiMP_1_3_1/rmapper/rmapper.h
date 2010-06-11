/*	$Id: rmapper.h 397 2009-10-06 04:02:47Z matei $	*/

#ifndef _RMAPPER_H
#define _RMAPPER_H

#include "../common/bitmap.h"
#include "../common/anchors.h"

extern const bool use_colours;
extern const bool use_dag;

#define DEF_MODE_50BP 0
#define DEF_MODE_35BP 1
#define DEF_MODE_70BP 2

#define DEF_MODE_FAST 0
#define DEF_MODE_SENSITIVE 1

/* default parameters - optimised for human */
//#define DEF_SPACED_SEED_CS	"1111001111"	/* handle more adjacencies */
//#define DEF_SPACED_SEED_LS	"111111011111"	/* longer for solexa/454 reads*/
//#define DEF_SPACED_SEED_DAG	"11110111"	/* shorter for Helicos */ 
#define	DEF_WINDOW_LEN		135.0		/* 115% of read length */
#define DEF_WINDOW_OVERLAP	80.0
#define DEF_NUM_MATCHES		4
#define DEF_HIT_TABOO_LEN	4
#define DEF_SEED_TABOO_LEN	0
#define DEF_NUM_OUTPUTS		100
#define DEF_MAX_READ_LEN	1000		/* high sanity mark */
#define DEF_KMER_STDDEV_LIMIT	-1		/* disabled by default */

#define DEF_MAX_N_DEFAULT_SEEDS 4
typedef char const * const default_seed_array_t[DEF_MAX_N_DEFAULT_SEEDS];

static int const default_seed_weight_cs = 12;
static int const default_seed_min_weight_cs = 8;
static int const default_seed_max_weight_cs = 12;
static int const default_seeds_cnt_cs[5] = { 1, 1, 4, 4, 4 };
static default_seed_array_t const default_seeds_cs[5]=
  { { "1111001111" },
    { "11111001111" },
    { "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" },
    { "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" },
    { "111110001111111", "111100111001001111", "111001001000111001111", "1111001000010001001001111" } };

static int const default_seed_weight_ls = 12;
static int const default_seed_min_weight_ls = 8;
static int const default_seed_max_weight_ls = 12;
static int const default_seeds_cnt_ls[5] = { 1, 1, 4, 4, 4 };
static default_seed_array_t const default_seeds_ls[5] =
  { { "111101111" },
    { "1111101111" },
    { "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" },
    { "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" },
    { "111101101011111", "111010110011001111", "1110110001011010111", "11110010100000100110111" } };

static int const default_seed_weight_hs = 7;
static int const default_seed_min_weight_hs = 7;
static int const default_seed_max_weight_hs = 7;
static int const default_seeds_cnt_hs[1] = { 1 };
static default_seed_array_t const default_seeds_hs[1] =
  { { "11110111" } };

/* DAG Scores/Parameters */
#define DEF_DAG_EPSILON		  0
#define DEF_DAG_READ_MATCH	  4
#define DEF_DAG_READ_GAP	 -2
#define DEF_DAG_READ_MISMATCH	 -4

#define DEF_DAG_REF_MATCH		11
#define DEF_DAG_REF_MISMATCH	       -10
#define DEF_DAG_REF_HALF_MATCH		4
#define DEF_DAG_REF_NEITHER_MATCH      -5
#define DEF_DAG_REF_MATCH_DELETION	5
#define DEF_DAG_REF_MISMATCH_DELETION  -6
#define DEF_DAG_REF_ERROR_INSERTION    -6
#define DEF_DAG_REF_WEIGHTED_THRESHOLD  8.0

/* SW Scores */
#define DEF_MATCH_VALUE		10
#define DEF_MISMATCH_VALUE	-15
#define DEF_A_GAP_OPEN		-40
#define DEF_B_GAP_OPEN		 DEF_A_GAP_OPEN
#define DEF_A_GAP_EXTEND	-7
#define DEF_B_GAP_EXTEND	 DEF_A_GAP_EXTEND
#define DEF_XOVER_PENALTY	-14	/* CS only */

#define DEF_MATCH_VALUE_DAG	 5
#define DEF_MISMATCH_VALUE_DAG	-6
#define DEF_A_GAP_OPEN_DAG	 0
#define DEF_B_GAP_OPEN_DAG	 0
#define DEF_A_GAP_EXTEND_DAG	-6
#define DEF_B_GAP_EXTEND_DAG	-3

#define DEF_N_SCORE_SETS 2
typedef int const default_score_array_t[7];	/* match, mismatch, go (db), ge (db), xover, go (qr), ge (qr) */
static default_score_array_t const default_score_set[DEF_N_SCORE_SETS] =
  { { 10, -15, -40, -7, -14, -40, -7 },		// about 7% seq errors in CS
    { 8, -30, -35, -7, -25, -35, -7 } };	// about 3% seq errors in CS
static int const default_score_set_id = 0;


#define DEF_SW_VECT_THRESHOLD	60.0	/* == DEF_SW_FULL_THRESHOLD in lspace */
#define DEF_SW_VECT_THRESHOLD_DAG 50.0	/* smaller for Helicos DAG */
#define DEF_SW_FULL_THRESHOLD	68.0	/* read_length x match_value x .68 */

#define DEF_ANCHOR_WIDTH	8	/* width around anchors in full SW */
#define DEF_HASH_FILTER_CALLS	true

/*
 * The maximum seed weight (maximum number of 1's in the seed) sets an
 * upper limit on our lookup table allocation size. The memory usage of
 * rmapper corresponds strongly to 4^MAX_SEED_WEIGHT * (sizeof(void *) +
 * sizeof(uint32_t)). At 16, this is 32GB on 32-bit and 48GB on 64-bit
 * architectures.
 */
#ifndef MAX_SEED_WEIGHT
#define MAX_SEED_WEIGHT		14
#endif

/*
 * We hold seeds as bitmaps to reduce cache footprint.
 */
#define MAX_SEED_SPAN		64

/*
 * For larger seeds we'll just use a hash table. Presently, we're restricted to
 * 128 bytes in kmer_to_mapidx, but it's trivially extended.
 */
#define MAX_HASH_SEED_WEIGHT	128
#define MAX_HASH_SEED_SPAN	128
#define HASH_TABLE_POWER	12	/* 4^HASH_POWER entries in table */

/* Sanity check - pin reads to a fairly small size */
#define MAX_READ_LENGTH	15000

struct re_score {
  struct read_entry *	parent;		/* associated read_elem */
  struct re_score  *	next;		/* linked list (final pass) */
  struct sw_full_results * sfrp;	/* alignment results (final pass) */
  struct anchor *	anchors;
  u_int			contig_num;	/* contig index (for filename)*/
  union {
    int32_t		score;		/* doubles as heap cnt in [0] */
    uint32_t		heap_elems;
  };
  union {
    uint32_t		index;		/* doubles as heap alloc in [0]*/
    uint32_t		heap_alloc;
  };
  bool			revcmpl;	/* from contig's reverse cmpl */
};

struct read_hit {
  uint32_t	g_idx_start;		/* kmer index in genome */
  uint8_t	r_idx_end_first;	/* start index of first kmer hit in read */
  uint8_t	r_idx_end_last;		/* end index of last kmer hit in read */
  uint8_t	more_than_once;		/* set to 1 if kmer hits read more than once */
  uint8_t	sn;			/* seed where kmer originates from */
  // aligned to 8B
};

/* the following are used during genome scan */
struct read_entry_scan {
  uint32_t	last_swhit_idx;		/* index of last sw call */
  uint16_t	window_len;		/* per-read window length */
  uint8_t	read_len;
  uint8_t	last_hit;		/* index in 'hits'; goes around */
  struct read_hit hits[0];	/* size depends on num_matches*/
};

struct cache_entry; // in cache.h

/*
 * Main entry for every read.
 */
struct read_entry {
  struct re_score * scores;	/* top 'num_ouputs' scores */
  char *	name;
  uint32_t *	read1;		/* the read as a bitstring */
  uint32_t *	read2;		/* second of Helicos pair */
  dag_cookie_t	dag_cookie;	/* kmer graph cookie for glue */
  struct cache_entry * cache;

  int16_t	initbp;		/* colour space init letter */
  uint8_t	cache_sz;
  uint8_t	head;		/* this works like a queue */

  uint32_t	read1_len;
  uint32_t	read2_len;
  uint32_t	offset;		/* offset in read array */
  int		swhits;		/* num of hits with sw */
  uint32_t	final_matches;	/* num of final output matches*/
#ifdef EXTRA_STATS		/* computed only ifdef EXTRA_STATS */
  count32_t	es_filter_calls;
  count32_t	es_filter_calls_bypassed;
  count32_t	es_filter_passes;
#endif
};

/*
 * Each index of `readmap' points to an array of `readmap_entry' structures.
 * These structures in turn point to all reads that contain the indexing kmer
 * as well as track where in the read the kmer exists for colinearity checking.
 */
struct readmap_entry {
  uint32_t	offset;		/* offset to the read in array */
  uint16_t	r_idx_end_first;/* start index of first kmer hit in read (could use uint8_t) */
  uint16_t	r_idx_end_last;	/* end index of last kmer hit in read */
  // aligned to 8B
};

struct seed_type {
  bitmap_type	mask[1];	/* a bitmask, least significant bit = rightmost match */
  uint32_t	span;		/* max 64 (could be uint8_t) */
  uint32_t	weight;		/* max 64 */
  // aligned to 8B
};


#endif
