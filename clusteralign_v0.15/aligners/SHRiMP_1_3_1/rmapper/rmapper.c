/*	$Id: rmapper.c 448 2009-11-23 23:33:55Z matei $	*/

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include <xmmintrin.h>	// for _mm_prefetch

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "../common/fasta.h"
#include "../common/dag_glue.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/sw-full-ls.h"
#include "../common/sw-vector.h"
#include "../common/output.h"
#include "../common/util.h"
#include "../common/version.h"
#include "../common/anchors.h"
#include "../common/hash.h"
#include "../common/stats.h"
#include "../rmapper/rmapper.h"
#include "../rmapper/cache.h"

static int	mode_read_length	= DEF_MODE_50BP;
static int	mode_speed_tradeoff	= DEF_MODE_FAST;

/* Seed management */
static struct seed_type *seed = NULL;
static uint32_t **seed_hash_mask = NULL;
static uint	n_seeds = 0;
static uint	max_seed_span = 0;
static int	selected_seed_weight;

/* External parameters */
static double	window_len		= DEF_WINDOW_LEN;
static double	window_overlap		= DEF_WINDOW_OVERLAP;
static u_int	num_matches		= DEF_NUM_MATCHES;
static int	hit_taboo_len		= DEF_HIT_TABOO_LEN;
static u_int	seed_taboo_len		= DEF_SEED_TABOO_LEN;
static u_int	num_outputs		= DEF_NUM_OUTPUTS;
static int	kmer_stddev_limit  	= DEF_KMER_STDDEV_LIMIT;

static u_int	dag_epsilon		= DEF_DAG_EPSILON;
static int	dag_read_match		= DEF_DAG_READ_MATCH;
static int	dag_read_mismatch	= DEF_DAG_READ_MISMATCH;
static int	dag_read_gap		= DEF_DAG_READ_GAP;

static int	dag_ref_match		= DEF_DAG_REF_MATCH;
static int 	dag_ref_mismatch	= DEF_DAG_REF_MISMATCH;
static int	dag_ref_half_match	= DEF_DAG_REF_HALF_MATCH;
static int	dag_ref_neither_match	= DEF_DAG_REF_NEITHER_MATCH;
static int	dag_ref_match_deletion	= DEF_DAG_REF_MATCH_DELETION;
static int	dag_ref_mismatch_deletion=DEF_DAG_REF_MISMATCH_DELETION;
static int	dag_ref_error_insertion	= DEF_DAG_REF_ERROR_INSERTION;
static double	dag_ref_weighted_thresh = DEF_DAG_REF_WEIGHTED_THRESHOLD;

static int	match_score		= DEF_MATCH_VALUE;
static int	mismatch_score		= DEF_MISMATCH_VALUE;
static int	a_gap_open_score	= DEF_A_GAP_OPEN;
static int	a_gap_extend_score	= DEF_A_GAP_EXTEND;
static int	b_gap_open_score	= DEF_B_GAP_OPEN;
static int	b_gap_extend_score	= DEF_B_GAP_EXTEND;
static int	crossover_score		= DEF_XOVER_PENALTY;
//static int	selected_score_set_id	= default_score_set_id;

static double	sw_vect_threshold	= DEF_SW_VECT_THRESHOLD;
static double	sw_full_threshold	= DEF_SW_FULL_THRESHOLD;
static int	anchor_width		= DEF_ANCHOR_WIDTH;
static bool	hash_filter_calls	= DEF_HASH_FILTER_CALLS;

/*
 * If window_len, sw_vect_threshold, sw_full_threshold are absolute values,
 * we'll set them negative to distinguish.
 */
#define IS_ABSOLUTE(x)	((x) < 0)

/* If x is negative, return its absolute value; else return base*x% */
static inline double
abs_or_pct(double x, double base) {
  if (IS_ABSOLUTE(x))
    return -x;
  else
    return base * (x / 100.0);
}

/* Genomic sequence, stored in 32-bit words, first is in the LSB */
static uint32_t *genome;			/* genome -- always in letter*/
static uint32_t *genome_cs;			/* genome -- colourspace */
static uint32_t	 genome_len;
static bool      genome_is_rna = false;		/* is genome RNA (has uracil)?*/
static char     *contig_name;			/* name of contig in genome */
static uint32_t  ncontigs;			/* number of contigs read */

/* Reads */
static struct read_entry_scan * reads_scan;
static struct read_entry * reads;		/* list of reads */
static uint32_t		      reads_allocated;	/* slots allocated for reads */
static uint32_t		      read_scan_size;	/* size of each read struct */

static inline struct read_entry_scan *
get_re_scan(uint32_t idx) {
  return (struct read_entry_scan *)((char *)reads_scan + (read_scan_size*idx));
}

static u_int	nreads;				/* total reads loaded */
static u_int	nkmers;				/* total kmers of reads loaded*/
static u_int	longest_read_len;		/* longest read we've seen */
static u_int	nduphits;			/* number of duplicate hits */

/* Kmer to read index */
static struct readmap_entry ***readmap;		/* kmer to read map */
static uint32_t **readmap_len;

/* Flags */
static int Bflag = false;			/* print a progress bar */
static int Cflag = false;			/* do complement only */
static int Fflag = false;			/* do positive (forward) only */
static int Hflag = false;			/* use hash table, not lookup */
static int Pflag = false;			/* pretty print results */
static int Rflag = false;			/* add read sequence to output*/
static int Tflag = false;			/* reverse sw full tie breaks */
static int Uflag = false;			/* output unmapped reads, too */

/* Scan stats */
static uint64_t shortest_scanned_kmer_list;
static uint64_t longest_scanned_kmer_list;
static uint64_t kmer_lists_scanned;
static uint64_t kmer_list_entries_scanned;

static count_t mem_readmap;
static count_t mem_reads;
static count_t mem_scores;


/* Extra stats, computed ifdef EXTRA_STATS */
#ifdef EXTRA_STATS

static count_t es_colin_checks;
static count_t es_total_filter_passes;
static count_t es_reads_filtered;
static count_t es_total_filter_calls_bypassed;

char const * reads_profile_file_name;

#define ES_count_add(c, v) count_add((c), (v))
#define ES_count32_add(c, v) count32_add((c), (v))
#define ES_count_increment(c) count_increment((c))
#define ES_count32_increment(c) count32_increment((c))

#else

#define ES_count_add(c, v) ((void)0)
#define ES_count32_add(c, v) ((void)0)
#define ES_count_increment(c) ((void)0)
#define ES_count32_increment(c) ((void)0)

#endif


/* Misc stats */
static uint64_t scan_usecs;
static uint64_t revcmpl_usecs;

/* kmer_to_mapidx function */
static uint32_t (*kmer_to_mapidx)(uint32_t *, u_int) = NULL;

#define PROGRESS_BAR(_a, _b, _c, _d)	\
    if (Bflag) progress_bar((_a), (_b), (_c), (_d))

#define FETCH_AHEAD 8
#define USE_PREFETCH


static size_t
power(size_t base, size_t exp)
{
  size_t result = 1;

  while (exp > 0) {
    if ((exp % 2) == 1)
      result *= base;
    base *= base;
    exp /= 2;
  } 

  return (result);
}

/* percolate down in our min-heap */
static void
reheap(struct re_score *scores, int node)
{
  struct re_score tmp;
  int left, right, max;

  assert(node >= 1 && node <= (int)scores[0].index);

  left  = node * 2;
  right = left + 1;
  max   = node;
	
  if (left <= scores[0].score &&
      scores[left].score < scores[node].score)
    max = left;

  if (right <= scores[0].score &&
      scores[right].score < scores[max].score)
    max = right;

  if (max != node) {
    tmp = scores[node];
    scores[node] = scores[max];
    scores[max] = tmp;
    reheap(scores, max);
  }
}

/* percolate up in our min-heap */
static void
percolate_up(struct re_score *scores, int node)
{
  struct re_score tmp;
  int parent;

  assert(node >= 1 && node <= (int)scores[0].index);

  if (node == 1)
    return;

  parent = node / 2;

  if (scores[parent].score > scores[node].score) {
    tmp = scores[node];
    scores[node] = scores[parent];
    scores[parent] = tmp;
    percolate_up(scores, parent);
  }
}

static char *
seed_to_string(uint sn)
{
  static char buffer[100];
  bitmap_type tmp;
  int i;

  assert(sn < n_seeds);

  buffer[seed[sn].span] = 0;
  for (i = seed[sn].span - 1, tmp = seed[sn].mask[0];
       i >= 0;
       i--, tmp >>= 1) {
    if (bitmap_extract(&tmp, 1, 0) == 1)
      buffer[i] = '1';
    else
      buffer[i] = '0';
  }

  return buffer;
}


/* translate hits from scan phase into anchors for full sw */
static void
save_anchors(struct read_entry_scan * res, uint32_t g_index,
	     struct anchor * anchors) {
  uint i, crt;

  assert(res != NULL && anchors != NULL);

  for (i = 0, crt = (res->last_hit + 1) % num_matches;
       i < num_matches;
       i++, crt = (crt + 1) % num_matches) {
    anchors[i].x = res->hits[crt].g_idx_start - ((int32_t)g_index);
    anchors[i].y = res->hits[crt].r_idx_end_first - (seed[res->hits[crt].sn].span - 1);
    anchors[i].length = seed[res->hits[crt].sn].span;
    anchors[i].width = 1;
    anchors[i].more_than_once = res->hits[crt].more_than_once;
    anchors[i].y_alt = res->hits[crt].r_idx_end_last - (seed[res->hits[crt].sn].span - 1);

    assert(anchors[i].x < res->window_len); //&& anchors[i].x >= 0);
    assert(anchors[i].y < res->read_len);
    assert(!anchors[i].more_than_once || anchors[i].y < anchors[i].y_alt);
  }
}


/* update our score min-heap */
static void
save_score(uint32_t read_id, int score, uint index, int contig_num,
	   bool revcmpl)
{
  struct read_entry * re = &reads[read_id];
  struct read_entry_scan * res = get_re_scan(read_id);
  struct re_score * scores = re->scores;

#ifdef DEBUG_HITS
  fprintf(stderr, "\tsaving hit: contig_num=%u read=%s idx=%u revcmpl=%s score=%u\n",
	  contig_num, re->name, index, revcmpl ? "yes" : "no", score);
#endif

  if (scores[0].heap_elems == num_outputs) {
    if (score < scores[1].score)
      return;

    /* replace min score */
    scores[1].score = score;
    scores[1].index = index;
    scores[1].contig_num = contig_num;
    scores[1].revcmpl = revcmpl;
    save_anchors(res, index, scores[1].anchors);

    reheap(scores, 1);
  } else {
    int idx = 1 + scores[0].heap_elems++;

    /* We do the array doubling trick for O(1) amortised alloc. */
    if (scores[0].heap_elems > scores[0].heap_alloc) {
      unsigned int alloc_slots, old_alloc_slots = scores[0].heap_alloc;

      assert(scores[0].heap_elems > 0);

      if (scores[0].heap_elems * 2 > num_outputs)
	alloc_slots = num_outputs;
      else
	alloc_slots = scores[0].heap_elems * 2;

      assert(alloc_slots <= num_outputs);

      scores[0].heap_alloc = alloc_slots;
      scores = (struct re_score *)
	xrealloc_c(scores, sizeof(struct re_score) * (alloc_slots + 1),
		   sizeof(struct re_score) * (old_alloc_slots + 1),
		   &mem_scores);
      memset(&scores[idx], 0, sizeof(struct re_score) * (alloc_slots + 1 - idx));
      re->scores = scores;
    }

    scores[idx].score = score;
    scores[idx].index = index;
    scores[idx].contig_num = contig_num;
    scores[idx].revcmpl = revcmpl;
    scores[idx].anchors = (struct anchor *)xcalloc_c(num_matches * sizeof(struct anchor),
						     &mem_scores);
    save_anchors(res, index, scores[idx].anchors);
    
    percolate_up(scores, idx);
  }
}

/* pulled off the web; this may or may not be any good */
static uint32_t
hash(uint32_t a)
{
  a = (a+0x7ed55d16) + (a<<12);
  a = (a^0xc761c23c) ^ (a>>19);
  a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);
  a = (a+0xfd7046c5) + (a<<3);
  a = (a^0xb55a4f09) ^ (a>>16);
  return (a);
}

/* hash-based version or kmer -> map index function for larger seeds */
static uint32_t
kmer_to_mapidx_hash(uint32_t *kmerWindow, u_int sn)
{
  static uint32_t maxidx = ((uint32_t)1 << 2*HASH_TABLE_POWER) - 1;

  uint32_t mapidx = 0;
  uint i;

  assert(seed_hash_mask != NULL);

  for (i = 0; i < BPTO32BW(max_seed_span); i++)
    mapidx = hash((kmerWindow[i] & seed_hash_mask[sn][i]) ^ mapidx);

  return mapidx & maxidx;
}


/*
 * Compress the given kmer into an index in 'readmap' according to the seed.
 * While not optimal, this is only about 20% of the spaced seed scan time.
 *
 * This is the original version for smaller seeds.
 *
 * XXX- This algorithm only considers bases 0-3, which implies overlap
 *      when we have other bases (mainly uracil, but also wobble codes).
 *      This won't affect sensitivity, but may cause extra S-W calls. 
 */
static uint32_t
kmer_to_mapidx_orig(uint32_t *kmerWindow, u_int sn)
{
  bitmap_type a = seed[sn].mask[0];
  uint32_t mapidx = 0;
  int i = seed[sn].span - 1;

  do {
    if ((a & 0x1) == 0x1) {
      mapidx <<= 2;
      mapidx |= ((kmerWindow[i/8] >> (i%8)*4) & 0x3);
    }
    a >>= 1;
    i--;

  } while (a != 0x0);

  assert(mapidx < power(4, seed[sn].weight));

  return mapidx;
}

/*
 * If 'kmer_stddev_limit' is set to a reasonable value, prune all kmers, whose
 * frequencies are more than 'kmer_stddev_limit' standard deviations greater
 * than the mean.
 */
static void
readmap_prune()
{
  uint sn;

  if (kmer_stddev_limit < 0)
    return;

  /* this won't work with the hash table approach... */
  assert(!Hflag);

  for (sn = 0; sn < n_seeds; sn++) {
    double mean, stddev;
    uint64_t pruned = 0;
    uint32_t i, readmapSize;

    readmapSize = ((uint32_t)1) << 2*seed[sn].weight; // 4^(seed_weight)

    mean = 0;
    for (i = 0; i < readmapSize; i++)
      mean += readmap_len[sn][i];
    mean /= readmapSize;

    stddev = 0;
    for (i = 0; i < readmapSize; i++)
      stddev += pow((double)readmap_len[sn][i] - mean, 2);
    stddev = sqrt(stddev / readmapSize);

    fprintf(stderr, "- Pruning kmers; mu: %f, sigma: %f, sigma limit: +/- %f\n",
	    mean, stddev, kmer_stddev_limit * stddev);
    if (mean < 1.0)
      fprintf(stderr, "WARNING: low mean - are you sure you want to prune kmers?\n");

    for (i = 0; i < readmapSize; i++) {
      if (readmap_len[sn][i] > mean + (kmer_stddev_limit * stddev)) {
	free(readmap[sn][i]);
	readmap[sn][i] = NULL;
	readmap_len[sn][i] = 0;
	pruned++;
      }
    }

    fprintf(stderr, "  - Pruned %" PRIu64 " kmer(s) from seed %d\n", pruned, sn);
  }
}

/* reset the fields used by scan() for each read */
static void
reset_reads()
{
  uint32_t i, j;

  for (i = 0; i < nreads; i++) {
    struct read_entry_scan * res = get_re_scan(i);

    res->last_swhit_idx = UINT32_MAX;
    res->last_hit = 0;
    for (j = 0; j < num_matches; j++)
      res->hits[j].g_idx_start = UINT32_MAX;

    /* Perhaps reset freq_hits? */
  }
}

static inline int
get_sw_threshold(struct read_entry * re, double which, int readnum)
{
  assert(readnum == 1
	 || (readnum == 2 && shrimp_mode == MODE_HELICOS_SPACE));

  if (readnum == 1)
    return (int)abs_or_pct(which, match_score * re->read1_len);
  else
    return (int)abs_or_pct(which, match_score * re->read2_len);
}


/*
 * Compute a hash value for this genome region.
 */
static inline cache_key_t
hash_window(uint32_t * scan_genome, uint goff, uint glen) {
  /*
  uint i, j;
  uint8_t base;
  uint32_t hash = 0;
  bitmap_type buffer;
  static uint const bases_per_bitmap = (8*sizeof(bitmap_type))/2;

  for (i = 0; i < ((glen - 1) / bases_per_bitmap) + 1; i++) {
    buffer = 0;
    for (j = 0; j < bases_per_bitmap && i*bases_per_bitmap + j < glen; j++) {
      base = EXTRACT(scan_genome, goff + i*bases_per_bitmap + j);
      bitmap_prepend(&buffer, 2, base & 0x03);
    }
    hash = SuperFastHash((char const *)&buffer, sizeof(bitmap_type), hash);
  }
  */

  /*
  uint i;
  uint32_t hash = 0;
  uint8_t base;

  for (i = 0; i < glen; i++) {
    base = EXTRACT(scan_genome, goff + i);
    hash ^= (base & 0x3) << i%16;
  }
  */

  uint i, j;
  uint8_t base;
  uint32_t key = 0;
  uint32_t buffer;
  static uint const bases_per_buffer = 16;

  for (i = 0; i < ((glen - 1) / bases_per_buffer) + 1; i++) {
    buffer = 0;
    for (j = 0; j < bases_per_buffer && i * bases_per_buffer + j < glen; j++) {
      base = EXTRACT(scan_genome, goff + i * bases_per_buffer + j);
      buffer <<= 2;
      buffer |= (base & 0x03);
    }
    hash_accumulate(&key, buffer);
  }
  hash_finalize(&key);

  return (cache_key_t)key;
}


/*
 * Determine whether or not the kmer hits seen within this read
 * were colinear.
 */ 
static inline int
are_hits_colinear(struct read_entry_scan *res)
{
  uint i, prev, crt;

  ES_count_increment(&es_colin_checks);

  if (shrimp_mode == MODE_HELICOS_SPACE) // no checking in HELICOS mode
    return 1;

  for (i = 1, prev = (res->last_hit + 1) % num_matches, crt = (prev + 1) % num_matches;
       i < num_matches;
       i++, prev = crt, crt = (crt + 1) % num_matches) {
    assert(res->hits[prev].g_idx_start + seed[res->hits[prev].sn].span - 1
	   <= res->hits[crt].g_idx_start + seed[res->hits[crt].sn].span - 1);

    //if (re->hits[prev].r_idx_start + (seed[re->hits[prev].sn].span - 1) > re->hits[crt].r_idx_end_last)
    if (!res->hits[prev].more_than_once && res->hits[prev].r_idx_end_last > res->hits[crt].r_idx_end_last)
      return 0;
  }

  return 1;
}

/* scan the genome by kmers, running S-W as needed, and updating scores */
static void
scan(int contig_num, bool revcmpl)
{
  struct read_entry_scan *res;
  struct readmap_entry *rme1, *rme2;
  uint32_t *kmerWindow, *scan_genome, *mapidxs, *mapidxs_p1;
  uint32_t i, j, k, pf, idx, sn;
  uint32_t base, skip, score1 = 0, score2, base_p1 = 0;
  uint32_t goff, glen, prevhit;
  u_int thresh;
#ifdef DEBUG_HITS
  uint32_t _i, _k;
#endif

  if (shrimp_mode == MODE_COLOUR_SPACE)
    scan_genome = genome_cs;
  else
    scan_genome = genome;

  PROGRESS_BAR(stderr, 0, 0, 10);

  kmerWindow = (uint32_t *)xcalloc(sizeof(kmerWindow[0]) * BPTO32BW(max_seed_span));
  mapidxs = (uint32_t *)xcalloc(sizeof(mapidxs[0]) * n_seeds);
  mapidxs_p1 = (uint32_t *)xcalloc(sizeof(mapidxs_p1[0]) * n_seeds);

  skip = 0;
  for (i = 0 ; i < genome_len && i < max_seed_span; i++) {
    base_p1 = EXTRACT(scan_genome, i);
    bitfield_prepend(kmerWindow, max_seed_span, base_p1);

    /*
     * Note, with this approach, when scanning a contig, each time we hit
     * an N or X base, we wait for the next max_seed_span bps before
     * we continue generating kmers. As a result, we might miss
     * kmers from seeds with smaller span. We do this to avoid another test.
     */
    if (base_p1 == BASE_N || base_p1 == BASE_X)
      skip = max_seed_span;
    if (skip > 0)
      skip--;
  }
  i--;

  // 1 ahead for i=seed_span
  sn = 0;
  do {
    mapidxs_p1[sn] = kmer_to_mapidx(kmerWindow, sn);
#ifdef USE_PREFETCH
    _mm_prefetch((const char *)&readmap[sn][mapidxs_p1[sn]], _MM_HINT_NTA);
#endif
  } while (++sn < n_seeds);

  /*
   * Main loop over current contig
   */
  for ( ; i < genome_len; i++) {
    PROGRESS_BAR(stderr, i, genome_len, 10);

    /*
     * shift 1 ahead to 0 ahead
     */
    base = base_p1;
    sn = 0;
    do {
      mapidxs[sn] = mapidxs_p1[sn];
    } while (++sn < n_seeds);

    /*
     * Prefetch start of readmap_entry list for seed 0
     */
#ifdef USE_PREFETCH
    _mm_prefetch((const char *)(readmap[0][mapidxs[0]]), _MM_HINT_NTA);
#endif

    /*
     * 1 ahead for i+1
     */
    if (i + 1 < genome_len) {
      base_p1 = EXTRACT(scan_genome, i+1);
      bitfield_prepend(kmerWindow, max_seed_span, base_p1);

      sn = 0;
      do {
	mapidxs_p1[sn] = kmer_to_mapidx(kmerWindow, sn);
#ifdef USE_PREFETCH
	_mm_prefetch((const char *)&readmap[sn][mapidxs_p1[sn]], _MM_HINT_NTA);
#endif
      } while (++sn < n_seeds);

    }

    /*
     * 0 ahead, real work
     */
    if (base == BASE_N || base == BASE_X)
      skip = max_seed_span - 1; // see previous comment about size of skip

    if (skip > 0) {
      skip--;
      continue;
    }

    /*
     * Main loop over seeds
     */
    sn = 0;
    do {
      rme1 = rme2 = readmap[sn][mapidxs[sn]];
      if (rme1 == NULL) {
	kmer_lists_scanned++;
	continue;
      }
      idx = i - (seed[sn].span - 1);

      /* prefetch first FETCH_AHEAD struct read_elem's */
#ifdef USE_PREFETCH
      for (pf = 0; rme2->offset != UINT32_MAX && pf < FETCH_AHEAD; pf++)
	_mm_prefetch((const char *)get_re_scan((rme2++)->offset), _MM_HINT_T0);
#endif

      for (res = get_re_scan(rme1->offset), j = k = 0;
	   rme1->offset != UINT32_MAX;
	   res = get_re_scan((++rme1)->offset), k++) {

	/*
	 * Continue prefetching the struct read_elem's,
	 * or start of readmap_entry list for next seed.
	 */
#ifdef USE_PREFETCH
	if (rme2->offset != UINT32_MAX)
	  _mm_prefetch((const char *)get_re_scan((rme2++)->offset), _MM_HINT_T0);
	else if (sn + 1 < n_seeds)
	  _mm_prefetch((const char *)readmap[sn+1][mapidxs[sn+1]], _MM_HINT_NTA);
#endif

	prevhit = res->hits[res->last_hit].g_idx_start;
	if (hit_taboo_len >= 0
	    && prevhit != UINT32_MAX
	    && idx - prevhit <= (uint)hit_taboo_len
	    ) {
#ifdef DEBUG_HITS
	  fprintf(stderr, "found potential hit i=%u:\n", i);
	  fprintf(stderr, "\tcontig_num=%u read=%s idx=%u r_idx_end_last=%u seed=%u\n",
		  contig_num, reads[rme1->offset].name, idx, rme1->r_idx_end_last, sn);
	  fprintf(stderr, "\tskipping, last_hit.g_idx_start=%u hit_taboo_len=%u\n",
		  res->hits[res->last_hit].g_idx_start, hit_taboo_len);
#endif
	  continue;
	}

	res->last_hit = (res->last_hit + 1) % num_matches;
	res->hits[res->last_hit].g_idx_start = idx;
	res->hits[res->last_hit].r_idx_end_first = rme1->r_idx_end_first;
	res->hits[res->last_hit].r_idx_end_last = rme1->r_idx_end_last;
	res->hits[res->last_hit].more_than_once = (rme1->r_idx_end_first != rme1->r_idx_end_last? 1 : 0);
	res->hits[res->last_hit].sn = sn;

	/*
	 * Compute start of window in genome for potential SW call.
	 */
	if (i < rme1->r_idx_end_first + (uint)(res->window_len - res->read_len)/2)
	  goff = 0;
	else
	  goff = i - (rme1->r_idx_end_first + (res->window_len - res->read_len)/2);

#ifdef DEBUG_HITS
	fprintf(stderr, "found potential hit i=%u:\n", i);
	fprintf(stderr, "\tcontig_num=%u read=%s idx=%u goff=%u r_idx_end_last=%u seed=%u\n",
		contig_num, reads[rme1->offset].name, idx, goff, rme1->r_idx_end_last, sn);
	fprintf(stderr, "\tprevious hits:\n");
	for (_i = 0, _k = (res->last_hit + 1) % num_matches;
	     _i < num_matches;
	     _i++, _k = (_k + 1) % num_matches) {
	  fprintf(stderr, "\tg_idx_start=%u r_idx_end_last=%u more_than_once=%u sn=%u\n",
		  res->hits[_k].g_idx_start, res->hits[_k].r_idx_end_last, res->hits[_k].more_than_once, res->hits[_k].sn);
	}
	if (res->hits[(res->last_hit + 1) % num_matches].g_idx_start == UINT32_MAX) {
	  fprintf(stderr, "\tskipping: 1st hit\n");
	} else if (goff > res->hits[(res->last_hit + 1) % num_matches].g_idx_start) {
	  fprintf(stderr, "\tskipping: previous hit too far back\n");
	} else if (res->last_swhit_idx != UINT32_MAX
		   && goff < res->last_swhit_idx + res->window_len/4) {
	  fprintf(stderr, "\tskipping: previous sw invocation too close\n");
	} else {
	  if (!are_hits_colinear(res)) {
	    fprintf(stderr, "\tskipping: hits not colinear\n");
	  }
	  ES_count_add(&es_colin_checks, -1);
	}
#endif

	if (res->hits[(res->last_hit + 1) % num_matches].g_idx_start != UINT32_MAX
	    && goff <= res->hits[(res->last_hit + 1) % num_matches].g_idx_start
	    && (res->last_swhit_idx == UINT32_MAX
		|| goff >= res->last_swhit_idx + res->window_len - (int)abs_or_pct(window_overlap, res->read_len))
	    && are_hits_colinear(res)
	    )
	  {
	    struct read_entry * re = &reads[rme1->offset];
	    bool meets_thresh = false;
	    cache_key_t key;
	    uint rel_pos, abs_pos;

	    assert(re->offset == rme1->offset);

	    ES_count32_increment(&re->es_filter_calls);

	    glen = res->window_len;
	    if (goff + glen > genome_len)
	      glen = genome_len - goff;

	    if (shrimp_mode == MODE_HELICOS_SPACE) {
	      score1 = sw_vector(scan_genome, goff, glen, re->read1, re->read1_len,
				 NULL, -1, genome_is_rna);
	      thresh = get_sw_threshold(re, sw_vect_threshold, 1);
	      if (score1 >= thresh) {
		score2 = sw_vector(scan_genome, goff, glen, re->read2, re->read2_len,
				   NULL, -1, genome_is_rna);
		thresh = get_sw_threshold(re, sw_vect_threshold, 2);
		if (score2 >= thresh) {
		  score1 = score1 * re->read1_len;
		  score2 = score2 * re->read2_len;
		  score1 = (score1 + score2) / (re->read1_len + re->read2_len);
		  meets_thresh = true;
		}
	      }
	    } else {
	      /*
	       * CS and LS
	       */

	      if (hash_filter_calls) // Use cache
		{
		  if (re->cache_sz == 0) {
		    ES_count_increment(&es_reads_filtered);

		    cache_init(re, &mem_reads);
		  }

		  key = hash_window(scan_genome, goff, glen);
		  //fprintf(stderr, "%08llX\n", (long long unsigned int)key);
		  cache_lookup(re, key, &rel_pos, &abs_pos);

		  if (cache_key_found(re, &rel_pos, &abs_pos)) {
		    /*
		     * We have run vector SW on this region before; use previous score.
		     */
		    score1 = re->cache[abs_pos].score;
		    re->cache[abs_pos].count++;

		    ES_count_increment(&es_total_filter_calls_bypassed);
		    ES_count32_increment(&re->es_filter_calls_bypassed);
		  } else {
		    if (shrimp_mode == MODE_COLOUR_SPACE) {
		      score1 = sw_vector(scan_genome, goff, glen, re->read1, re->read1_len,
					 genome, re->initbp, genome_is_rna);
		    } else if (shrimp_mode == MODE_LETTER_SPACE) {
		      score1 = sw_vector(scan_genome, goff, glen, re->read1, re->read1_len,
					 NULL, -1, genome_is_rna);
		    }

		    /* save score */
		    if (cache_is_full(re, &rel_pos, &abs_pos))
		      cache_maybe_resize(re, &mem_reads);

		    cache_add_key(re, &rel_pos, &abs_pos);
		    re->cache[abs_pos].key = key;
		    re->cache[abs_pos].count = 1;
		    re->cache[abs_pos].score = score1;
		  }
		}
	      else  // Don't use cache
		{
		  if (shrimp_mode == MODE_COLOUR_SPACE) {
		    score1 = sw_vector(scan_genome, goff, glen, re->read1, re->read1_len,
				       genome, re->initbp, genome_is_rna);
		  } else if (shrimp_mode == MODE_LETTER_SPACE) {
		    score1 = sw_vector(scan_genome, goff, glen, re->read1, re->read1_len,
				       NULL, -1, genome_is_rna);
		  }
		}

	      thresh = get_sw_threshold(re, sw_vect_threshold, 1);
	      if (score1 >= thresh)
		meets_thresh = true;
	    } // CS and LS

	    if (meets_thresh) {
	      ES_count_increment(&es_total_filter_passes);
	      ES_count32_increment(&re->es_filter_passes);

	      save_score(re->offset, score1, goff, contig_num, revcmpl);
	      res->last_swhit_idx = goff;
	      re->swhits++;
#ifdef DEBUG_HITS
	    } else {
	      fprintf(stderr, "\thit doesn't meet SW vector threshold: %d<%d\n",
		      score1, thresh);
#endif
	    }
	  }
	j++;
      } // loop over read_elems

      assert(k == readmap_len[sn][mapidxs[sn]]);

      kmer_list_entries_scanned += j;
      shortest_scanned_kmer_list = MIN(shortest_scanned_kmer_list, j);
      longest_scanned_kmer_list = MAX(longest_scanned_kmer_list, j);
      kmer_lists_scanned++;

    } while (++sn < n_seeds);// loop over spaced seeds
  } // loop over genome

  if (Bflag) {
    PROGRESS_BAR(stderr, genome_len, genome_len, 10);
    putc('\n', stderr);
  }

  free(kmerWindow);
  free(mapidxs);
  free(mapidxs_p1);
}

static uint32_t
readalloc()
{
  struct read_entry_scan *res;
  struct read_entry *re;
  size_t bytes, bytes_scan, prevbytes, prevbytes_scan, i;

  assert(nreads <= reads_allocated);

  if (nreads == reads_allocated) {
    read_scan_size = sizeof(*res) +
      (sizeof(res->hits[0]) * num_matches);

    prevbytes_scan = reads_allocated * read_scan_size;
    prevbytes = reads_allocated * sizeof(struct read_entry);

    reads_allocated += 1000;

    bytes_scan = reads_allocated * read_scan_size;
    bytes = reads_allocated * sizeof(struct read_entry);

    reads_scan = (struct read_entry_scan *)xrealloc_c(reads_scan, bytes_scan,
						      prevbytes_scan,
						      &mem_reads);
    memset((char *)reads_scan + prevbytes_scan, 0, bytes_scan - prevbytes_scan);

    reads = (struct read_entry *)xrealloc_c(reads, bytes,
					    prevbytes,
					    &mem_reads);
    memset((char *)reads + prevbytes, 0, bytes - prevbytes);
  }

  res = get_re_scan(nreads);
  re = &reads[nreads];

  re->offset = nreads;

  for (i = 0; i < num_matches; i++)
    res->hits[i].g_idx_start = UINT32_MAX;

  bytes = sizeof(struct re_score) * 1;
  re->scores = (struct re_score *)xcalloc_c(bytes,
					    &mem_scores);

  nreads++;

  return nreads-1;
}

static bool
load_reads_lscs(const char *file)
{
  struct read_entry_scan *res;
  struct read_entry *re;
  uint32_t read_id;
  char *name, *seq;
  fasta_t fasta;
  size_t seqlen;
  ssize_t bases = 0;
  int space;

  if (shrimp_mode == MODE_LETTER_SPACE)
    space = LETTER_SPACE;
  else
    space = COLOUR_SPACE;
  fasta = fasta_open(file, space);
  if (fasta == NULL) {
    fprintf(stderr, "error: failed to open reads file [%s]\n", file);
    return (false);
  } else {
    fprintf(stderr, "- Processing reads file [%s]\n", file);
  }

  fprintf(stderr, "- Loading reads...");

  while (true) {
    u_int i, load;
    uint32_t *kmerWindow;

    if (Bflag && (nreads % 17) == 0) 
      fprintf(stderr, "\r- Loading reads... %u", nreads);

    if (!fasta_get_next(fasta, &name, &seq, NULL))
      break;

    if (strchr(name, '\t') != NULL || strchr(seq, '\t') != NULL) {
      fprintf(stderr, "error: tabs are not permitted in fasta names "
	      "or sequences. Tag: [%s].\n", name);
      exit(1);
    }

    seqlen = strlen(seq);
    if (seqlen == 0) {
      fprintf(stderr, "error: read [%s] had no sequence!\n",
	      name);
      exit(1);
    }
    if (seqlen > MAX_READ_LENGTH) {
      fprintf(stderr, "error: read [%s] had unreasonable "
	      "length!\n", name);
    }
    if (shrimp_mode == MODE_COLOUR_SPACE) {
      /* the sequence begins with the initial letter base */
      if (seqlen < 1) {
	fprintf(stderr, "error: read [%s] had sequence "
		"with no colours!\n", name);
	exit(1);
      }
      seqlen--;
    }
    bases += seqlen;

    read_id = readalloc();
    re = &reads[read_id];
    res = get_re_scan(read_id);

    res->last_swhit_idx = UINT32_MAX;
    re->name = name;

    re->read1 = fasta_sequence_to_bitfield(fasta, seq);
    if (re->read1 == NULL) {
      fprintf(stderr, "error: invalid sequence; tag: [%s]\n", name); 
      exit(1);
    }
    re->read1_len = seqlen;
    res->read_len = seqlen;
    longest_read_len = MAX(re->read1_len, longest_read_len);

    if (shrimp_mode == MODE_COLOUR_SPACE)
      re->initbp = fasta_get_initial_base(fasta, seq);

#ifdef DEBUG_KMERS
    fprintf(stderr, "processing read:\"%s\" length:%u sequence:\"%s\" bitfield:%s\n",
	    re->name, res->read_len, seq,
	    bitmap32v_string(re->read1, res->read_len, 4, true, false, false));
#endif

    kmerWindow = (uint32_t *)xcalloc(sizeof(kmerWindow[0])*BPTO32BW(max_seed_span));

    for (i = 0, load = 0; i < seqlen; i++) {
      uint base, sn;

      base = EXTRACT(re->read1, i);
      bitfield_prepend(kmerWindow, max_seed_span, base);

      /*
       * Ignore kmers that contain an 'N'.
       *
       * load = number of bases inside kmer
       */
      if (base == BASE_N || base == BASE_X)
	load = 0;
      else if (load < max_seed_span)
	load++;

      /*
       * Disabled seed_taboo; no skipping
       */
      //if (skip > 0) {
      //  skip--;
      //  continue;
      //}

      for (sn = 0; sn < n_seeds; sn++) {
	if (load < seed[sn].span)
	  continue;

	/*
	 * For simplicity we throw out the first kmer when in colour space. If
	 * we did not do so, we'd run into a ton of colour-letter space
	 * headaches. For instance, how should the first read kmer match
	 * against a kmer from the genome? The first colour of the genome kmer
	 * depends on the previous letter in the genome, so we may have a
	 * matching read, but the colour representation doesn't agree due to
	 * different initialising bases.
	 *
	 * If we wanted to be complete, we could compute the four permutations
	 * and add them, but I'm not so sure that'd be a good idea. Perhaps
	 * this should be investigated in the future.
	 */
	if (shrimp_mode == MODE_COLOUR_SPACE && i == seed[sn].span - 1)
	  continue;

	uint32_t mapidx = kmer_to_mapidx(kmerWindow, sn);

#ifdef DEBUG_KMERS
	fprintf(stderr, "\tfound kmer:%s",
		bitmap32v_string(&mapidx, seed[sn].weight, 2, false, false, false));
	fprintf(stderr, " in kmerWindow:%s using seed:%u\n",
		bitmap32v_string(kmerWindow, max_seed_span, 4, false, false, false),
		sn);
#endif

	/* permit only one kmer reference per read; but update r_idx_end_last */
	if (readmap[sn][mapidx] != NULL &&
	    readmap[sn][mapidx][readmap_len[sn][mapidx] - 1].offset == re->offset) {
	  readmap[sn][mapidx][readmap_len[sn][mapidx] - 1].r_idx_end_last = i;
	  continue;
	}
				
	readmap_len[sn][mapidx]++;
	readmap[sn][mapidx] = (struct readmap_entry *)
	  xrealloc_c(readmap[sn][mapidx], sizeof(struct readmap_entry)*(readmap_len[sn][mapidx] + 1),
		     sizeof(struct readmap_entry)*readmap_len[sn][mapidx],
		     &mem_readmap);
	readmap[sn][mapidx][readmap_len[sn][mapidx] - 1].offset = re->offset;
	readmap[sn][mapidx][readmap_len[sn][mapidx] - 1].r_idx_end_first = i;
	readmap[sn][mapidx][readmap_len[sn][mapidx] - 1].r_idx_end_last = i;
	readmap[sn][mapidx][readmap_len[sn][mapidx]].offset = UINT32_MAX;

	nkmers++;
      }

      /*
       * Taboo the seed generation as well, if requested.
       * Disabled; for multiple seeds it should be done on a per-seed basis.
       */
      //if (seed_taboo_len)
      //skip = seed_taboo_len;
    }

    free(seq);
    free(kmerWindow);
  }
  if (Bflag)
    fprintf(stderr, "\r- Loading reads... %u\n", nreads);
  else
    fprintf(stderr, "\n");

  fprintf(stderr, "- Loaded %s %s in %s reads (%s kmers)\n",
	  comma_integer(bases), (shrimp_mode == MODE_COLOUR_SPACE) ?
	  "colours" : "letters", comma_integer(nreads),comma_integer(nkmers));

  fasta_close(fasta);

  return (true);
}

static bool
load_reads_dag(const char *file)
{
  struct read_entry_scan *res;
  struct read_entry *re;
  uint32_t read_id;
  char *name1, *name2, *seq1, *seq2;
  fasta_t fasta;
  size_t seq1len, seq2len;
  ssize_t bases = 0;
  uint sn;

  fasta = fasta_open(file, LETTER_SPACE);
  if (fasta == NULL) {
    fprintf(stderr, "error: failed to open reads file [%s]\n",
	    file);
    return (false);
  }

  fprintf(stderr, "- Loading read pairs...");

  while (true) {
    int i;
    bool b1, b2;
    char ** kmers;

    if (Bflag && (nreads % 17) == 0) 
      fprintf(stderr, "\r- Loading read pairs... %u", nreads);

    b1 = fasta_get_next(fasta, &name1, &seq1, NULL);
    b2 = fasta_get_next(fasta, &name2, &seq2, NULL);

    if (b1 == false)
      break;

    if (b2 == false) {
      fprintf(stderr, "error: helicos fasta file [%s] does "
	      "not appear to have a pair for read [%s]\n",
	      file, name1);
      return (false);
    }

    if (strchr(name1, '\t') != NULL || strchr(seq1, '\t') != NULL) {
      fprintf(stderr, "error: tabs are not permitted in fasta names "
	      "or sequences. Tag: [%s].\n", name1);
      exit(1);
    }
    if (strchr(name2, '\t') != NULL || strchr(seq2, '\t') != NULL) {
      fprintf(stderr, "error: tabs are not permitted in fasta names "
	      "or sequences. Tag: [%s].\n", name2);
      exit(1);
    }

    seq1len = strlen(seq1);
    seq2len = strlen(seq2);
    bases += seq1len + seq2len;

    if (seq1len == 0 || seq2len == 0) {
      fprintf(stderr, "error: read [%s] had no sequence!\n",
	      name1);
      exit(1);
    }
    if (seq1len > MAX_READ_LENGTH || seq2len > MAX_READ_LENGTH) {
      fprintf(stderr, "error: read [%s] had unreasonable "
	      "length!\n", name1);
    }

    read_id = readalloc();
    re = &reads[read_id];
    res = get_re_scan(read_id);

    res->last_swhit_idx = UINT32_MAX;
    re->name = name1;
    re->read1 = fasta_sequence_to_bitfield(fasta, seq1);
    if (re->read1 == NULL) {
      fprintf(stderr, "error: invalid sequence; tag: [%s]\n", name1);
      exit(1);
    }
    re->read2 = fasta_sequence_to_bitfield(fasta, seq2);
    if (re->read2 == NULL) {
      fprintf(stderr, "error: invalid sequence; tag: [%s]\n", name2);
      exit(1);
    }
    re->read1_len = seq1len;
    res->read_len = seq1len;
    re->read2_len = seq2len;

    longest_read_len = MAX(re->read1_len, longest_read_len);
    longest_read_len = MAX(re->read2_len, longest_read_len);

    re->dag_cookie = dag_build_kmer_graph(seq1, seq2, dag_epsilon);

    /*
     * Should this be called for every seed???
     */
    for (sn = 0; sn < n_seeds; sn++) {
      kmers = dag_get_kmers(re->dag_cookie, seed[sn].span);
      assert(kmers != NULL);

      for (i = 0; kmers[i] != NULL; i++) {
	uint32_t *kmer;

	/* scan() builds kmers in reverse, so we must as well... */
	strrev(kmers[i]);

	kmer = fasta_sequence_to_bitfield(fasta, kmers[i]);
	if (kmer == NULL) {
	  fprintf(stderr, "error: invalid hs kmer: [%s]\n", kmers[i]);
	  exit(1);
	}

	uint32_t mapidx = kmer_to_mapidx(kmer, sn);

	/* permit only one kmer reference per read; but update r_idx_end_last */
	if (readmap[sn][mapidx] != NULL &&
	    readmap[sn][mapidx][readmap_len[sn][mapidx] - 1].offset == re->offset) {
	  readmap[sn][mapidx][readmap_len[sn][mapidx] - 1].r_idx_end_last = (res->read_len*3)/4; // ???
	  continue;
	}
				
	readmap_len[sn][mapidx]++;
	readmap[sn][mapidx] = (struct readmap_entry *)
	  xrealloc_c(readmap[sn][mapidx], sizeof(struct readmap_entry)*(readmap_len[sn][mapidx] + 1),
		     sizeof(struct readmap_entry)*readmap_len[sn][mapidx],
		     &mem_readmap);
	readmap[sn][mapidx][readmap_len[sn][mapidx] - 1].offset = re->offset;
	readmap[sn][mapidx][readmap_len[sn][mapidx] - 1].r_idx_end_first = (res->read_len*3)/4; // ???
	readmap[sn][mapidx][readmap_len[sn][mapidx] - 1].r_idx_end_last = (res->read_len*3)/4; // ???
	readmap[sn][mapidx][readmap_len[sn][mapidx]].offset = UINT32_MAX;

	nkmers++;
	free(kmer);
      }

      if (i == 0)
	fprintf(stderr, "WARNING: Got %d kmers for [%s] and [%s]\n", i, name1, name2);

      dag_free_kmers(kmers);
    }

    free(name2);
    free(seq1);
    free(seq2);
  }

  if (Bflag)
    fprintf(stderr, "\r- Loading read pairs... %u\n", nreads);
  else
    fprintf(stderr, "\n");

  fprintf(stderr, "- Loaded %s letters in %u read pairs (%u kmers)\n",
	  comma_integer(bases), nreads, nkmers);

  fasta_close(fasta);

  return true;
}

static bool
load_reads(const char *file)
{
  uint sn;
  ssize_t ret, map_size;

  // allocate as many readmaps as we have seeds
  readmap = (struct readmap_entry ***)xmalloc_c(n_seeds * sizeof(readmap[0]),
						&mem_readmap);
  readmap_len = (uint32_t **)xmalloc_c(n_seeds * sizeof(readmap_len[0]),
				       &mem_readmap);

  for (sn = 0; sn < n_seeds; sn++) {
    if (Hflag)
      map_size = power(4, HASH_TABLE_POWER);
    else
      map_size = power(4, seed[sn].weight);

    readmap[sn] = (struct readmap_entry **)xcalloc_c(map_size * sizeof(readmap[sn][0]),
						     &mem_readmap);
    readmap_len[sn] = (uint32_t *)xcalloc_c(map_size * sizeof(readmap_len[sn][0]),
					    &mem_readmap);
  }

  if (shrimp_mode == MODE_HELICOS_SPACE)
    ret = load_reads_dag(file);
  else
    ret = load_reads_lscs(file);

  return ret;
}

/*
 * XXX - this needs to be hooked into our regular output format.
 */
static void
generate_output_dag(struct re_score *rs_array, size_t rs_len, bool revcmpl)
{
  uint32_t goff, glen;
  size_t i;

  /* Obtain final results for each element of our array. */
  for (i = 0; i < rs_len; i++) {
    struct re_score * rs = &rs_array[i];
    struct read_entry * re = rs->parent;
    struct read_entry_scan * res = get_re_scan(re->offset);
    double score;
    bool meets_thresh = false;

    goff = rs->index;
    assert(goff < genome_len);

    glen = res->window_len;
    if (goff + glen > genome_len)
      glen = genome_len - goff;

    u_int j;
    char *refseq = (char *)xmalloc(glen + 1);

    for (j = 0; j < glen; j++) {
      refseq[j] = base_translate(EXTRACT(genome, goff + j), false);
    }
    refseq[j] = '\0';

    struct dag_alignment *dap = dag_build_alignment(refseq, re->dag_cookie);

    double avglen;

    avglen = (double)(re->read1_len + re->read2_len) / 2;
    if (((double)dap->score / avglen) >= dag_ref_weighted_thresh) {
      meets_thresh = true;
      score = ((double)dap->score / avglen);
    }

    if (meets_thresh) {
      uint32_t genome_start = goff + dap->start_index + 1;
      uint32_t genome_end = goff + dap->end_index;

      if (revcmpl) {
	uint32_t tmp = genome_start;
	genome_start = genome_len - genome_end - 1;
	genome_end = genome_len - tmp - 1;
      }

      printf("\n--------------------------------------------------------------------------------\n");
      printf("ALIGNMENT:\n");
      printf("  contig: [%s]\n", contig_name);
      printf("  read:   [%s]\n", re->name);
      printf("  strand: %c, score %.3f, genome_start: %u, genome_end: %u\n",
	     revcmpl ? '-' : '+', score, genome_start, genome_end);
      printf("\n");
      printf("  READ 1:     %s\n", dap->read1);
      printf("  READ 2:     %s\n", dap->read2);
      printf("  SEQUENCE:   %s\n", dap->sequence);
      printf("--------------------------------------------------------------------------------\n");
      printf("\n");
      re->final_matches++;
    }

    free(refseq);
    dag_free_alignment(dap);
  }
}

static int
score_cmp(const void *arg1, const void *arg2)
{
  const struct re_score *one = (const struct re_score *)arg1;
  const struct re_score *two = (const struct re_score *)arg2;

  assert(one->revcmpl == two->revcmpl);

  if (one->score > two->score)
    return (-1);
  if (one->score < two->score)
    return (1);

  if (one->sfrp->genome_start > two->sfrp->genome_start)
    return (1);
  if (one->sfrp->genome_start < two->sfrp->genome_start)
    return (-1);

  if (one->sfrp->matches > two->sfrp->matches)
    return (-1);
  if (one->sfrp->matches < two->sfrp->matches)
    return (1);

  return (0);
}

static void
generate_output_lscs(struct re_score *rs_array, size_t rs_len, bool revcmpl)
{
  static bool firstcall = true;

  struct sw_full_results last_sfr;
  uint32_t goff, glen;
  int thresh;
  u_int i;

  assert(shrimp_mode == MODE_LETTER_SPACE || shrimp_mode == MODE_COLOUR_SPACE);

  if (rs_len == 0)
    return;

  /* shut up, gcc */
  thresh = 0;

  /* Obtain final results for each element of our array. */
  for (i = 0; i < rs_len; i++) {
    struct re_score * rs = &rs_array[i];
    struct read_entry * re = rs->parent;
    struct read_entry_scan * res = get_re_scan(re->offset);

    goff = rs->index;
    assert(goff < genome_len);

    glen = res->window_len;
    if (goff + glen > genome_len)
      glen = genome_len - goff;

    thresh = (int)abs_or_pct(sw_full_threshold, match_score * re->read1_len);

    /*
     * Hacky Optimisation:
     *   In letter-space, our full alignment will be the same as
     *   what our vector aligner detected. For this reason, we
     *   can call the vector aligner a few more times to home in
     *   on the alignment we want and reduce our genome length
     *   in the much slower, full aligner, providing a net win.
     */
    if (shrimp_mode == MODE_LETTER_SPACE) {
      uint32_t tryoff, trylen;
      int score;

      trylen = glen / 2;
      tryoff = goff + trylen / 2;
      assert(tryoff < (goff + glen));
      score = sw_vector(genome, tryoff, trylen, re->read1, re->read1_len,
			NULL, -1, genome_is_rna);
      if (score == rs->score) {
	goff = tryoff;
	glen = trylen;
      }
    }

    rs->sfrp = (struct sw_full_results *)xmalloc(sizeof(*rs->sfrp));

    if (shrimp_mode == MODE_COLOUR_SPACE) {
      sw_full_cs(genome, goff, glen, re->read1, re->read1_len,
		 re->initbp, thresh, rs->sfrp, revcmpl && Tflag, genome_is_rna,
		 rs->anchors, num_matches);
    } else {
      /*
       * The full SW in letter space assumes it's given the correct max score.
       * This might not be true just yet if we're using hashing&caching because
       * of possible hash collosions.
       */
      rs->score = sw_vector(genome, goff, glen, re->read1, re->read1_len,
			    NULL, -1, genome_is_rna);
      if (rs->score >= thresh) {
	sw_full_ls(genome, goff, glen, re->read1, re->read1_len,
		   thresh, rs->score, rs->sfrp, revcmpl && Tflag,
		   rs->anchors, num_matches);
	assert(rs->sfrp->score == rs->score);
      } else { // this wouldn't have passed the filter; eliminated in loop below
	rs->sfrp->score = rs->score;
	rs->sfrp->dbalign = NULL;
	rs->sfrp->qralign = NULL;
      }
    }
    rs->score = rs->sfrp->score;
  }

  /* Sort the scores. */
  qsort(rs_array, rs_len, sizeof(rs_array[0]), score_cmp);

  /* Output sorted list, removing any duplicates. */
  memset(&last_sfr, 0, sizeof(last_sfr));
  for (i = 0; i < rs_len; i++) {
    struct re_score * rs = &rs_array[i];
    struct read_entry * re = rs->parent;
    bool dup;

    dup = sw_full_results_equal(&last_sfr, rs->sfrp);
    if (dup)
      nduphits++;

    if (rs->score >= thresh && dup == false) {
      char *output;

      re->final_matches++;

      if (firstcall) {
	output = output_format_line(Rflag);
	puts(output);
	free(output);
      }
      firstcall = false;

      output = output_normal(re->name, contig_name, rs->sfrp,
			     genome_len, (shrimp_mode == MODE_COLOUR_SPACE), re->read1,
			     re->read1_len, re->initbp, revcmpl, Rflag);
      puts(output);
      free(output);

      if (Pflag) {
	output = output_pretty(re->name, contig_name, rs->sfrp,
			       genome, genome_len,
			       (shrimp_mode == MODE_COLOUR_SPACE), re->read1,
			       re->read1_len, re->initbp, revcmpl);
	puts("");
	puts(output);
	free(output);
      }
    }

    last_sfr = *rs->sfrp;
    if (rs->sfrp->dbalign != NULL)
      free(rs->sfrp->dbalign);
    if (rs->sfrp->qralign != NULL)
      free(rs->sfrp->qralign);
    free(rs->sfrp);
  }
}

/*
 * Split our linked list of scores into segments based on read and send an
 * array of them to the appropriate space-specific output function.
 */
static void
generate_output(struct re_score *rs, bool revcmpl)
{
  struct re_score *byread;
  struct read_entry *re, *last_re;
  u_int i;

  byread = (struct re_score *)xcalloc(num_outputs * sizeof(*byread));

  i = 0;
  last_re = NULL;
  while (true) {
    re = (rs == NULL) ? NULL : rs->parent;

    if ((rs == NULL && byread[0].score > 0) ||
	(re != NULL && last_re != NULL && re != last_re)) {
      if (shrimp_mode == MODE_HELICOS_SPACE)
	generate_output_dag(byread, i, revcmpl);
      else
	generate_output_lscs(byread, i, revcmpl);

      memset(byread, 0, sizeof(*byread) * num_outputs);
      i = 0;
    }

    if (rs == NULL)
      break;

    assert(re != NULL);
    assert(i < num_outputs);
    byread[i++] = *rs;
    last_re = re;

    rs = rs->next;
  }

  free(byread);
}

static void
final_pass_contig(struct re_score *scores, struct re_score *scores_revcmpl)
{

  if (Fflag)
    generate_output(scores, false);

  if (Cflag) {
    uint64_t before;

    before = gettimeinusecs();
    reverse_complement(genome, NULL, genome_len, genome_is_rna);
    revcmpl_usecs += (gettimeinusecs() - before);
    generate_output(scores_revcmpl, true);
  }
}

static int
final_pass(char **files, u_int nfiles)
{
  struct re_score **scores, **scores_revcmpl;
  u_int i, j, hits = 0;
  fasta_t fasta;
  char *name, *sequence;
  bool is_rna = false;

  scores = (struct re_score **)xcalloc(ncontigs * sizeof(*scores));
  scores_revcmpl = (struct re_score**)xcalloc(ncontigs * sizeof(*scores));

  /* For each contig, generate a linked list of hits from it. */
  for (i = 0; i < nreads; i++) {
    struct read_entry * re = &reads[i];

    if (re->swhits == 0)
      continue;

    hits++;

    for (j = 1; j <= (u_int)re->scores[0].score; j++) {
      struct re_score *rs = &re->scores[j];
      u_int cn = rs->contig_num;

      assert(rs->score >= 0);
      assert(cn < ncontigs);
      assert(rs->parent == NULL);
      assert(rs->next == NULL);

      /* extra sanity */
#ifndef NDEBUG
      if (shrimp_mode != MODE_HELICOS_SPACE) {
	int thresh = get_sw_threshold(re, sw_vect_threshold, 1);
	assert(rs->score >= thresh);
      }
#endif

      /* prepend to linked list */
      if (rs->revcmpl) {
	rs->next = scores_revcmpl[cn];
	scores_revcmpl[cn] = &re->scores[j];
      } else {
	rs->next = scores[cn];
	scores[cn] = &re->scores[j];
      }

      rs->parent = re;
    }
  }

  /* Now, do a final s-w for all reads of each contig and output them. */
  for (i = j = 0; i < nfiles; i++) {
    fprintf(stderr, "- Processing contig file [%s] (%u of %u)\n",
	    files[i], i + 1, nfiles);

    fasta = fasta_open(files[i], LETTER_SPACE);
    if (fasta == NULL) {
      fprintf(stderr, "error: failed to reload genome file "
	      "[%s] for final pass\n", files[i]);
      exit(1);
    }

    while (fasta_get_next(fasta, &name, &sequence, &is_rna)) {
      assert(j < ncontigs);

      if (strchr(name, '\t') != NULL || strchr(sequence, '\t') != NULL) {
	fprintf(stderr, "error: tabs are not permitted in fasta names "
		"or sequences. Tag: [%s].\n", name);
	exit(1);
      }

      contig_name = name;
      genome_len = strlen(sequence);
      genome_is_rna = is_rna;

      genome = fasta_sequence_to_bitfield(fasta, sequence);
      if (genome == NULL) {
	fprintf(stderr, "error: invalid contig sequence; contig: [%s]\n", name);
	exit(1);
      }
      genome_cs = NULL;

      fprintf(stderr, "  - Loaded %s letters from contig \"%s\"\n",
	      comma_integer(genome_len), name);

      free(sequence);
      final_pass_contig(scores[j], scores_revcmpl[j]);

      free(genome);
      free(contig_name);
      genome = genome_cs = NULL;
      contig_name = NULL;
      j++;
    }

    fasta_close(fasta);
  }

  free(scores);
  free(scores_revcmpl);

  return (hits);
}

static bool
valid_spaced_seeds()
{
  uint sn;

  for (sn = 0; sn < n_seeds; sn++) {
    if (seed[sn].weight > MAX_SEED_WEIGHT && !Hflag)
      return false;

    if (Hflag && (seed[sn].span > MAX_HASH_SEED_SPAN ||
		  seed[sn].weight > MAX_HASH_SEED_WEIGHT))
      return false;
  }

  return true;
}

static void
scan_genomes_contig()
{
  uint64_t before;

  if (Fflag) {
    /*
     * Do it forwards.
     */
    before = gettimeinusecs();
    scan(ncontigs, false);
    scan_usecs += (gettimeinusecs() - before);
    reset_reads();
  }

  if (Cflag) {
    /*
     * Do it reverse-complementwards.
     */
    fprintf(stderr, "    - Processing reverse complement\n");

    before = gettimeinusecs();
    if (shrimp_mode == MODE_COLOUR_SPACE)
      reverse_complement(genome, genome_cs, genome_len, genome_is_rna);
    else
      reverse_complement(genome, NULL, genome_len, genome_is_rna);
    revcmpl_usecs += (gettimeinusecs() - before);

    before = gettimeinusecs();
    scan(ncontigs, true);
    scan_usecs += (gettimeinusecs() - before);
    reset_reads();
  }

  ncontigs++;
}

static int
scan_genomes(char **files, int nfiles)
{
  fasta_t fasta;
  char *name, *sequence;
  bool is_rna = false;
  int i;

  for (i = 0; i < nfiles; i++) {
    fprintf(stderr, "- Processing contig file [%s] (%d of %d)\n",
	    files[i], i + 1, nfiles);

    fasta = fasta_open(files[i], LETTER_SPACE);
    if (fasta == NULL) {
      fprintf(stderr, "error: failed to reload genome file "
	      "[%s] for final pass\n", files[i]);
      return (1);
    }

    while (fasta_get_next(fasta, &name, &sequence, &is_rna)) { 
      contig_name = name;
      genome_len = strlen(sequence);
      genome_is_rna = is_rna;
 
      if (strchr(name, '\t') != NULL || strchr(sequence, '\t') != NULL) {
	fprintf(stderr, "error: tabs are not permitted in fasta names "
		"or sequences. Tag: [%s].\n", name);
	exit(1);
      }

      genome = fasta_sequence_to_bitfield(fasta, sequence);
      if (genome == NULL) {
	fprintf(stderr, "error: invalid contig sequence; contig: [%s]\n", name);
	exit(1);
      }
      if (shrimp_mode == MODE_COLOUR_SPACE)
	genome_cs = fasta_bitfield_to_colourspace(fasta, genome, genome_len, genome_is_rna);

      fprintf(stderr, "  - Loaded %s letters from contig \"%s\"\n",
	      comma_integer(genome_len), name);

      free(sequence);
      scan_genomes_contig();

      free(genome);
      if (shrimp_mode == MODE_COLOUR_SPACE)
	free(genome_cs);
      free(contig_name);
      genome = genome_cs = NULL;
      contig_name = NULL;
    }

    fasta_close(fasta);
  }

  return (0);
}

static u_int
set_window_lengths()
{
  uint32_t i;
  u_int max_window_len = 0;

  if (!IS_ABSOLUTE(window_len) && window_len < 100.0) {
    fprintf(stderr, "WARNING: window length is < 100%% of read "
	    "length - is this what you want?\n");
  }

  for (i = 0; i < nreads; i++) {
    struct read_entry * re = &reads[i];
    struct read_entry_scan * res = get_re_scan(i);

    res->window_len = (uint16_t)abs_or_pct(window_len, MAX(re->read1_len, re->read2_len));

    max_window_len = MAX(max_window_len, res->window_len);
  }

  return max_window_len;
}

static void
print_statistics()
{
  double cellspersec, hz, seedscantime;
  uint64_t invocs, cells, reads_matched, total_matches, fasta_load_ticks, vticks;
  fasta_stats_t fs;
  uint32_t i;

  fs = fasta_stats();
  fasta_load_ticks = fs->total_ticks;
  free(fs);
  hz = cpuhz();

  reads_matched = total_matches = 0;
  for (i = 0; i < nreads; i++) {
    struct read_entry * re = &reads[i];

    reads_matched += (re->final_matches == 0) ? 0 : 1;
    total_matches += re->final_matches;
  }

  sw_vector_stats(&invocs, &cells, &vticks, &cellspersec);

  seedscantime = ((double)scan_usecs / 1000000.0) - (vticks / hz);
  seedscantime = MAX(0, seedscantime);

  fprintf(stderr, "\nStatistics:\n");
  fprintf(stderr, "    Spaced Seed Scan:\n");
  fprintf(stderr, "        Run-time:               %.2f seconds\n",
	  seedscantime); 
  fprintf(stderr, "        Total Kmers:            %s\n",
	  comma_integer(nkmers));
  fprintf(stderr, "        Minimal Reads/Kmer:     %s\n",
	  comma_integer(shortest_scanned_kmer_list));
  fprintf(stderr, "        Maximal Reads/Kmer:     %s\n",
	  comma_integer(longest_scanned_kmer_list));
  fprintf(stderr, "        Kmer List Entries:      %s\n",
	  comma_integer(kmer_list_entries_scanned));
  fprintf(stderr, "        Average Reads/Kmer:     %.2f\n",
	  (kmer_lists_scanned == 0) ? 0 :
	  (double)kmer_list_entries_scanned / (double)kmer_lists_scanned);
#ifdef EXTRA_STATS
  fprintf(stderr, "        Colinearity Checks:     %s\n",
	  comma_integer(count_get_count(&es_colin_checks)));
#endif
  fprintf(stderr, "\n");

  fprintf(stderr, "    Vector Smith-Waterman:\n");
  fprintf(stderr, "        Run-time:               %.2f seconds\n",
	  (cellspersec == 0) ? 0 : cells / cellspersec);
  fprintf(stderr, "        Invocations:            %s\n",
	  comma_integer(invocs));
#ifdef EXTRA_STATS
  fprintf(stderr, "        Passed threshold:       %s\n",
          comma_integer(count_get_count(&es_total_filter_passes)));
  fprintf(stderr, "        Reads filtered:         %s\n",
	  comma_integer(count_get_count(&es_reads_filtered)));
  fprintf(stderr, "        Bypassed calls:         %s\n",
	  comma_integer(count_get_count(&es_total_filter_calls_bypassed)));
#endif
  fprintf(stderr, "        Cells Computed:         %.2f million\n",
	  (double)cells / 1.0e6);
  fprintf(stderr, "        Cells per Second:       %.2f million\n",
	  cellspersec / 1.0e6);
  fprintf(stderr, "\n");

  if (shrimp_mode != MODE_HELICOS_SPACE) {
    if (shrimp_mode == MODE_COLOUR_SPACE)
      sw_full_cs_stats(&invocs, &cells, NULL, &cellspersec);
    else
      sw_full_ls_stats(&invocs, &cells, NULL, &cellspersec);

    fprintf(stderr, "    Scalar Smith-Waterman:\n");
    fprintf(stderr, "        Run-time:               %.2f seconds\n",
	    (cellspersec == 0) ? 0 : cells / cellspersec);
    fprintf(stderr, "        Invocations:            %s\n",
	    comma_integer(invocs));
    fprintf(stderr, "        Cells Computed:         %.2f million\n",
	    (double)cells / 1.0e6);
    fprintf(stderr, "        Cells per Second:       %.2f million\n",
	    cellspersec / 1.0e6);
  } else {
    struct dag_statistics *dsp = dag_get_statistics();

    fprintf(stderr, "    DAG Read-Read Alignment:\n");
    fprintf(stderr, "        Run-time:               %.2f seconds\n",
	    dsp->kmers_seconds);
    fprintf(stderr, "        Avg Kmers/Read Pair:    %.2f\n",
	    (double)dsp->kmers_total_kmers / (double)dsp->kmers_invocations);
    fasta_load_ticks = (uint64_t)MAX(0, (float)fasta_load_ticks -
				     (float)dsp->kmers_seconds * hz);
    fprintf(stderr, "\n");

    fprintf(stderr, "    DAG Pair-Genome Alignment:\n");
    fprintf(stderr, "        Run-time:               %.2f seconds\n",
	    dsp->aligner_seconds);
    fprintf(stderr, "        Invocations:            %s\n",
	    comma_integer(dsp->aligner_invocations));
    free(dsp);
  }
  fprintf(stderr, "\n");

  fprintf(stderr, "    Miscellaneous:\n");
  fprintf(stderr, "        Load Time:              %.2f seconds\n",
	  (double)fasta_load_ticks / hz);
  fprintf(stderr, "        Revcmpl. Time:          %.2f seconds\n",
	  (double)revcmpl_usecs / 1000000.0);
  fprintf(stderr, "\n");

  fprintf(stderr, "    General:\n");
  fprintf(stderr, "        Reads Matched:          %s    "
	  "(%.4f%%)\n", comma_integer(reads_matched),
	  (nreads == 0) ? 0 : ((double)reads_matched / (double)nreads) * 100);
  fprintf(stderr, "        Total Matches:          %s\n",
	  comma_integer(total_matches));
  fprintf(stderr, "        Avg Hits/Matched Read:  %.2f\n",
	  (total_matches == 0) ? 0 : ((double)total_matches /
				      (double)reads_matched));
  fprintf(stderr, "        Duplicate Hits Pruned:  %s\n",
	  comma_integer(nduphits));

#ifdef EXTRA_STATS
  FILE * reads_profile_file;
  if (reads_profile_file_name != NULL
      && (reads_profile_file = fopen(reads_profile_file_name, "w")) != NULL) {

    for (i = 0; i < nreads; i++) {
      uint l = 0;

      if (reads[i].cache_sz != 0) {
	uint crt = reads[i].head;
	while (l < reads[i].cache_sz
	       && reads[i].cache[crt].count != 0) {
	  l++;
	  crt = (crt + reads[i].cache_sz - 1) % reads[i].cache_sz;
	}
      }
      fprintf(reads_profile_file, "%-20s\t%u\t%u\t%u\t%u\t%u\n",
	      reads[i].name,
	      reads[i].es_filter_calls,
	      reads[i].es_filter_calls_bypassed,
	      l,
	      reads[i].es_filter_passes,
	      reads[i].final_matches);
    }
    fclose(reads_profile_file);
  }
#endif

  fprintf(stderr, "\n    Memory usage:\n");
  fprintf(stderr, "        Readmap:                %s\n",
	  comma_integer(count_get_count(&mem_readmap)));
  fprintf(stderr, "        Reads:                  %s\n",
	  comma_integer(count_get_count(&mem_reads)));
  fprintf(stderr, "        Scores:                 %s\n",
	  comma_integer(count_get_count(&mem_scores)));
}

static bool
add_spaced_seed(const char *seedStr)
{
  uint i;

  seed = (struct seed_type *)xrealloc(seed, sizeof(struct seed_type) * (n_seeds + 1));
  seed[n_seeds].mask[0] = 0x0;
  seed[n_seeds].span = strlen(seedStr);
  seed[n_seeds].weight = strchrcnt(seedStr, '1');

  if (seed[n_seeds].span < 1
      || seed[n_seeds].span > MAX_SEED_SPAN
      || seed[n_seeds].weight < 1
      || strchrcnt(seedStr, '0') != seed[n_seeds].span - seed[n_seeds].weight)
    return false;

  for (i = 0; i < seed[n_seeds].span; i++)
    bitmap_prepend(seed[n_seeds].mask, 1, (seedStr[i] == '1' ? 1 : 0));

  if (seed[n_seeds].span > max_seed_span)
    max_seed_span = seed[n_seeds].span;

  n_seeds++;
  return true;
}


/*
 * Load default seeds of a given weight. If passed -1, load default seeds for default weight.
 * Returns 1 if there are no default seeds of given weight, 0 otherwise.
 */
static int
load_default_seeds(int requested_weight) {
  int i;
  int default_seed_min_weight = 0, default_seed_max_weight = 0;
  int const * default_seeds_cnt = NULL;
  default_seed_array_t * default_seeds = NULL;

  switch(shrimp_mode) {
  case MODE_COLOUR_SPACE:
    default_seed_min_weight = default_seed_min_weight_cs;
    default_seed_max_weight = default_seed_max_weight_cs;
    default_seeds_cnt = default_seeds_cnt_cs;
    default_seeds = default_seeds_cs;
    break;
  case MODE_LETTER_SPACE:
    default_seed_min_weight = default_seed_min_weight_ls;
    default_seed_max_weight = default_seed_max_weight_ls;
    default_seeds_cnt = default_seeds_cnt_ls;
    default_seeds = default_seeds_ls;
    break;
  case MODE_HELICOS_SPACE:
    default_seed_min_weight = default_seed_min_weight_hs;
    default_seed_max_weight = default_seed_max_weight_hs;
    default_seeds_cnt = default_seeds_cnt_hs;
    default_seeds = default_seeds_hs;
    break;
  }
  
  if (requested_weight == -1)
    requested_weight = selected_seed_weight;
  if (requested_weight < default_seed_min_weight
      || requested_weight > default_seed_max_weight)
    return 1;
  for (i = 0; i < default_seeds_cnt[requested_weight - default_seed_min_weight]; i++)
    add_spaced_seed(default_seeds[requested_weight - default_seed_min_weight][i]);

  return 0;
}


static void
init_seed_hash_mask(void)
{
  uint sn;
  int i;

  seed_hash_mask = (uint32_t **)xmalloc(sizeof(seed_hash_mask[0])*n_seeds);
  for (sn = 0; sn < n_seeds; sn++) {
    seed_hash_mask[sn] = (uint32_t *)xcalloc(sizeof(seed_hash_mask[sn][0])*BPTO32BW(max_seed_span));

    for (i = seed[sn].span - 1; i >= 0; i--)
      bitfield_prepend(seed_hash_mask[sn], max_seed_span,
		       bitmap_extract(seed[sn].mask, 1, i) == 1? 0xf : 0x0);
  }
}


static void
usage(char * progname, bool full_usage)
{
  char *slash;
  double vect_sw_threshold = -1;
  uint sn;

  switch (shrimp_mode) {
  case MODE_LETTER_SPACE:
    vect_sw_threshold = DEF_SW_VECT_THRESHOLD;
    break;
  case MODE_COLOUR_SPACE:
    vect_sw_threshold = DEF_SW_VECT_THRESHOLD;
    break;
  case MODE_HELICOS_SPACE:
    vect_sw_threshold = DEF_SW_VECT_THRESHOLD_DAG;
    break;
  default:
    assert(0);
  }

  if (n_seeds == 0)
    load_default_seeds(-1);

  slash = strrchr(progname, '/');
  if (slash != NULL)
    progname = slash + 1;

  fprintf(stderr, "usage: %s [parameters] [options] "
	  "reads_file genome_file1 genome_file2...\n", progname);

  fprintf(stderr, "Parameters:\n");

  fprintf(stderr,
	  "    -s    Spaced Seed(s) (weight/span)            (default: ");
  for (sn = 0; sn < n_seeds; sn++) {
    if (sn > 0)
      fprintf(stderr, "                                                            ");
    fprintf(stderr, "%s (%u/%u)%s\n", seed_to_string(sn), seed[sn].weight, seed[sn].span,
	    (sn == n_seeds - 1? ")" : ","));
  }

  fprintf(stderr,
	  "    -n    Seed Matches per Window                 (default: %d)\n",
	  DEF_NUM_MATCHES);

  fprintf(stderr,
	  "    -t    Seed Hit Taboo Length                   (default: %d)\n",
	  DEF_HIT_TABOO_LEN);

  if (shrimp_mode != MODE_HELICOS_SPACE) {
    fprintf(stderr,
	    "    -9    Seed Generation Taboo Length            (default: %d)\n",
	    DEF_SEED_TABOO_LEN);
  }

  fprintf(stderr,
	  "    -w    Seed Window Length                      (default: %.02f%%)\n",
	  DEF_WINDOW_LEN);

  fprintf(stderr,
	  "    -W    Seed Window Overlap Length              (default: %.02f%%)\n",
	  DEF_WINDOW_OVERLAP);

  fprintf(stderr,
	  "    -o    Maximum Hits per Read                   (default: %d)\n",
	  DEF_NUM_OUTPUTS);

  /*
  fprintf(stderr,
	  "    -r    Maximum Read Length                     (default: %d)\n",
	  DEF_MAX_READ_LEN);
  */

  fprintf(stderr,
	  "    -d    Kmer Std. Deviation Limit               (default: %d%s)"
	  "\n", DEF_KMER_STDDEV_LIMIT,
	  (DEF_KMER_STDDEV_LIMIT < 0) ? " [None]" : "");

  if (shrimp_mode == MODE_HELICOS_SPACE) {
    fprintf(stderr, "\n");
    fprintf(stderr,
	    "    -p    DAG Epsilon                             (default: %d)\n",
	    DEF_DAG_EPSILON);

    fprintf(stderr,
	    "    -1    DAG Read Match                          (default: %d)\n",
	    DEF_DAG_READ_MATCH);

    fprintf(stderr,
	    "    -y    DAG Read Mismatch                       (default: %d)\n",
	    DEF_DAG_READ_MISMATCH);

    fprintf(stderr,
	    "    -z    DAG Read Gap                            (default: %d)\n",
	    DEF_DAG_READ_GAP);

    fprintf(stderr,
	    "    -a    DAG Reference Match                     (default: %d)\n",
	    DEF_DAG_REF_MATCH);

    fprintf(stderr,
	    "    -b    DAG Reference Mismatch                  (default: %d)\n",
	    DEF_DAG_REF_MISMATCH);

    fprintf(stderr,
	    "    -c    DAG Reference Half Match                (default: %d)\n",
	    DEF_DAG_REF_HALF_MATCH);

    fprintf(stderr,
	    "    -j    DAG Reference Neither Match             (default: %d)\n",
	    DEF_DAG_REF_NEITHER_MATCH);

    fprintf(stderr,
	    "    -k    DAG Reference Match,Deletion            (default: %d)\n",
	    DEF_DAG_REF_MATCH_DELETION);

    fprintf(stderr,
	    "    -l    DAG Reference Mismatch,Deletion         (default: %d)\n",
	    DEF_DAG_REF_MISMATCH_DELETION);

    fprintf(stderr,
	    "    -u    DAG Reference Error,Insertion           (default: %d)\n",
	    DEF_DAG_REF_ERROR_INSERTION);

    fprintf(stderr,
	    "    -2    DAG Reference Weighted Threshold        (default: %.3f)\n",
	    DEF_DAG_REF_WEIGHTED_THRESHOLD);
  }

  fprintf(stderr, "\n");
  fprintf(stderr,
	  "    -m    S-W Match Value                         (default: %d)\n",
	  (shrimp_mode == MODE_HELICOS_SPACE) ? DEF_MATCH_VALUE_DAG : DEF_MATCH_VALUE);

  fprintf(stderr,
	  "    -i    S-W Mismatch Value                      (default: %d)\n",
	  (shrimp_mode == MODE_HELICOS_SPACE) ? DEF_MISMATCH_VALUE_DAG : DEF_MISMATCH_VALUE);

  fprintf(stderr,
	  "    -g    S-W Gap Open Penalty (Reference)        (default: %d)\n",
	  (shrimp_mode == MODE_HELICOS_SPACE) ? DEF_A_GAP_OPEN_DAG : DEF_A_GAP_OPEN);

  fprintf(stderr,
	  "    -q    S-W Gap Open Penalty (Query)            (default: %d)\n",
	  (shrimp_mode == MODE_HELICOS_SPACE) ? DEF_B_GAP_OPEN_DAG : DEF_B_GAP_OPEN);

  fprintf(stderr,
	  "    -e    S-W Gap Extend Penalty (Reference)      (default: %d)\n",
	  (shrimp_mode == MODE_HELICOS_SPACE) ? DEF_A_GAP_EXTEND_DAG : DEF_A_GAP_EXTEND);

  fprintf(stderr,
	  "    -f    S-W Gap Extend Penalty (Query)          (default: %d)\n",
	  (shrimp_mode == MODE_HELICOS_SPACE) ? DEF_B_GAP_EXTEND_DAG : DEF_B_GAP_EXTEND);

  if (shrimp_mode == MODE_COLOUR_SPACE) {
    fprintf(stderr,
	    "    -x    S-W Crossover Penalty                   ("
	    "default: %d)\n", DEF_XOVER_PENALTY);
  }

  if (shrimp_mode == MODE_COLOUR_SPACE) {
    fprintf(stderr,
	    "    -h    S-W Full Hit Threshold                  "
	    "(default: %.02f%%)\n", DEF_SW_FULL_THRESHOLD);
  }
  if (shrimp_mode == MODE_COLOUR_SPACE || shrimp_mode == MODE_HELICOS_SPACE) {
    fprintf(stderr,
	    "    -v    S-W Vector Hit Threshold                "
	    "(default: %.02f%%)\n", vect_sw_threshold);
  } else {
    fprintf(stderr,
	    "    -h    S-W Hit Threshold                       "
	    "(default: %.02f%%)\n", DEF_SW_FULL_THRESHOLD);
  }

  if (full_usage) {
    fprintf(stderr, "\n");
    if (shrimp_mode == MODE_COLOUR_SPACE || shrimp_mode == MODE_LETTER_SPACE) {
      fprintf(stderr,
	      "    -A    Anchor width limiting full SW           (default: %d; disable: -1)\n",
	      DEF_ANCHOR_WIDTH);
    }
#ifdef EXTRA_STATS
    //fprintf(stderr,
	//    "    -D    Dump per-read scan statistics in file\n");
#endif
    fprintf(stderr,
	    "    -Y    Cache sizes for genome scan             (default: %d,%d)\n",
	    cache_min_size, cache_max_size);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");

  fprintf(stderr,
	  "    -B    Print Scan Progress Bar                       (default: "
	  "disabled)\n"); 

  fprintf(stderr,
	  "    -C    Only Process Negative Strand (Rev. Compl.)    (default: "
	  "disabled)\n");

  fprintf(stderr,
	  "    -F    Only Process Positive Strand                  (default: "
	  "disabled)\n");

  fprintf(stderr,
	  "    -P    Pretty Print Alignments                       (default: "
	  "disabled)\n"); 

  fprintf(stderr,
	  "    -R    Print Reads in Output                         (default: "
	  "disabled)\n");

  if (shrimp_mode != MODE_HELICOS_SPACE) {
    fprintf(stderr,
	    "    -T    Reverse tie-break choice on negative strand   (default: "
	    "disabled)\n");
  }

  fprintf(stderr,
	  "    -U    Print Unmapped Read Names in Output           (default: "
	  "disabled)\n");
  fprintf(stderr,
	  "    -?    Full list of parameters and options\n");

  if (full_usage) {
    fprintf(stderr, "\n");
    fprintf(stderr,
	    "    -Z    Toggle hashing of genome windows              (default: %s)\n",
	    hash_filter_calls? "enabled" : "disabled");
  }

  exit(1);
}


static int
set_mode_from_string(char const * s) {
  if (!strcmp(s, "fast")) {
    mode_speed_tradeoff = DEF_MODE_FAST;
  } else if (!strcmp(s, "sensitive")) {
    mode_speed_tradeoff = DEF_MODE_SENSITIVE;
  } else if (!strcmp(s, "35bp")) {
    mode_read_length = DEF_MODE_35BP;
  } else if (!strcmp(s, "50bp")) {
    mode_read_length = DEF_MODE_50BP;
  } else if (!strcmp(s, "70bp")) {
    mode_read_length = DEF_MODE_70BP;
  } else
    return 0;
  return 1;
}


static void
set_params_from_mode(bool num_matches_set, bool window_len_set, bool hit_taboo_len_set) {
  switch (mode_speed_tradeoff) {
  case DEF_MODE_FAST:

    if (n_seeds == 0)
      load_default_seeds(12);
    switch (mode_read_length) {
    case DEF_MODE_35BP:

      if (!num_matches_set)
	num_matches = 3;
      if (!window_len_set)
	window_len = 140.0;
      if (!hit_taboo_len_set)
	hit_taboo_len = 2;
      break;

    case DEF_MODE_50BP:

      if (!num_matches_set)
	num_matches = 4;
      if (!window_len_set)
	window_len = 135.0;
      if (!hit_taboo_len_set)
	hit_taboo_len = 3;
      break;

    case DEF_MODE_70BP:

      if (!num_matches_set)
	num_matches = 6;
      if (!window_len_set)
	window_len = 130.0;
      if (!hit_taboo_len_set)
	hit_taboo_len = 4;
      break;

    }
    break;

  case DEF_MODE_SENSITIVE:
    
    if (n_seeds == 0)
      load_default_seeds(10);
    anchor_width = -1;
    switch (mode_read_length) {
    case DEF_MODE_35BP:

      if (!num_matches_set)
	num_matches = 2;
      if (!window_len_set)
	window_len = 160.0;
      if (!hit_taboo_len_set)
	hit_taboo_len = -1;
      break;

    case DEF_MODE_50BP:

      if (!num_matches_set)
	num_matches = 2;
      if (!window_len_set)
	window_len = 150.0;
      if (!hit_taboo_len_set)
	hit_taboo_len = -1;
      break;

    case DEF_MODE_70BP:

      if (!num_matches_set)
	num_matches = 3;
      if (!window_len_set)
	window_len = 140.0;
      if (!hit_taboo_len_set)
	hit_taboo_len = -1;
      break;

    }
    break;
  }
}


int
main(int argc, char **argv)
{
  char *reads_file, * progname = argv[0];
  char const * optstr = "";
  char **genome_files;
  int i, ch, ret, ngenome_files;
  u_int max_window_len;
  bool a_gap_open_score_set = false, b_gap_open_score_set = false;
  bool a_gap_extend_score_set = false, b_gap_extend_score_set = false;
  bool match_score_set = false, mismatch_score_set = false, crossover_score_set = false;
  bool num_matches_set = false, window_len_set = false, hit_taboo_len_set = false;
  char *c;
  int requested_weight;

  set_mode_from_argv(argv);

  /* set the appropriate defaults based on mode */
  switch (shrimp_mode) {
  case MODE_COLOUR_SPACE:
    optstr = "?s:n:t:9:w:o:r:d:m:i:g:q:e:f:x:h:v:BCFHPRTUA:ZY:W:M:";
    selected_seed_weight = default_seed_weight_cs;
    break;
  case MODE_LETTER_SPACE:
    optstr = "?s:n:t:9:w:o:r:d:m:i:g:q:e:f:h:X:BCFHPRTUA:ZY:W:M:";
    selected_seed_weight = default_seed_weight_ls;
    break;
  case MODE_HELICOS_SPACE:
    optstr = "?s:n:t:w:o:r:d:p:1:y:z:a:b:c:j:k:l:u:2:m:i:g:q:e:f:v:BCFHPRUZY:W:";
    selected_seed_weight = default_seed_weight_hs;
    match_score = DEF_MATCH_VALUE_DAG;
    mismatch_score = DEF_MISMATCH_VALUE_DAG;
    a_gap_open_score = DEF_A_GAP_OPEN_DAG;
    b_gap_open_score = DEF_B_GAP_OPEN_DAG;
    a_gap_extend_score = DEF_A_GAP_EXTEND_DAG;
    b_gap_extend_score = DEF_B_GAP_EXTEND_DAG;
    sw_vect_threshold = DEF_SW_VECT_THRESHOLD_DAG;
    break;
  default:
    assert(0);
  }

  fprintf(stderr, "--------------------------------------------------"
	  "------------------------------\n");
  fprintf(stderr, "rmapper: %s.\nSHRiMP %s\n[%s]\n", get_mode_string(),
	  SHRIMP_VERSION_STRING, get_compiler());
  fprintf(stderr, "--------------------------------------------------"
	  "------------------------------\n");

  while ((ch = getopt(argc, argv, optstr)) != -1) {
    switch (ch) {
    case 's':
      if (optarg[0] == 'w') { // request default seeds of given weight
	requested_weight = atoi(&optarg[1]);
	if (load_default_seeds(requested_weight)) {
	  fprintf(stderr, "error: no default seeds of given weight (%d)\n", requested_weight);
	  exit(1);
	}
      } else {
	c = strtok(optarg, ",");
	do {
	  if (!add_spaced_seed(c)) {
	    fprintf(stderr, "error: invalid spaced seed \"%s\"\n", c);
	    exit (1);
	  }
	  c = strtok(NULL, ",");
	} while (c != NULL);
      }
      break;
    case 'n':
      num_matches_set = true;
      num_matches = atoi(optarg);
      break;
    case 't':
      hit_taboo_len_set = true;
      hit_taboo_len = atoi(optarg);
      break;
    case '9':
      assert(shrimp_mode != MODE_HELICOS_SPACE);
      seed_taboo_len = atoi(optarg);
      break;
    case 'w':
      window_len_set = true;
      window_len = atof(optarg);
      if (window_len <= 0.0) {
	fprintf(stderr, "error: invalid window "
		"length\n");
	exit(1);
      }
      if (strcspn(optarg, "%.") == strlen(optarg))
	window_len = -window_len;		//absol.
      break;
    case 'o':
      num_outputs = atoi(optarg);
      break;
    case 'r':
      /* max_read_len is gone - silently ignore */
      break;
    case 'd':
      kmer_stddev_limit = atoi(optarg);
      break;
    case 'p':
      dag_epsilon = atoi(optarg);
      break;
    case '1':
      dag_read_match = atoi(optarg);
      break;
    case 'y':
      dag_read_mismatch = atoi(optarg);
      break;
    case 'z':
      dag_read_gap = atoi(optarg);
      break;
    case 'a':
      dag_ref_match = atoi(optarg);
      break;
    case 'b':
      dag_ref_mismatch = atoi(optarg);
      break;
    case 'c':
      dag_ref_half_match = atoi(optarg);
      break;
    case 'j':
      dag_ref_neither_match = atoi(optarg);
      break;
    case 'k':
      dag_ref_match_deletion = atoi(optarg);
      break;
    case 'l':
      dag_ref_mismatch_deletion = atoi(optarg);
      break;
    case 'u':
      dag_ref_error_insertion = atoi(optarg);
      break;
    case '2':
      dag_ref_weighted_thresh = atof(optarg);
      break;
    case 'm':
      match_score_set = true;
      match_score = atoi(optarg);
      break;
    case 'i':
      mismatch_score_set = true;
      mismatch_score = atoi(optarg);
      break;
    case 'g':
      a_gap_open_score_set = true;
      a_gap_open_score = atoi(optarg);
      break;
    case 'q':
      b_gap_open_score_set = true;
      b_gap_open_score = atoi(optarg);
      break;
    case 'e':
      a_gap_extend_score_set = true;
      a_gap_extend_score = atoi(optarg);
      break;
    case 'f':
      b_gap_extend_score_set = true;
      b_gap_extend_score = atoi(optarg);
      break;
    case 'x':
      assert(shrimp_mode == MODE_COLOUR_SPACE);
      crossover_score_set = true;
      crossover_score = atoi(optarg);
      break;
    case 'h':
      assert(shrimp_mode != MODE_HELICOS_SPACE);
      sw_full_threshold = atof(optarg);
      if (sw_full_threshold < 0.0) {
	fprintf(stderr, "error: invalid s-w full "
		"hit threshold\n");
	exit(1);
      }
      if (strcspn(optarg, "%.") == strlen(optarg))
	sw_full_threshold = -sw_full_threshold;	//absol.
      break;
    case 'v':
      assert(shrimp_mode == MODE_COLOUR_SPACE ||
	     shrimp_mode == MODE_HELICOS_SPACE);
      sw_vect_threshold = atof(optarg);
      if (sw_vect_threshold < 0.0) {
	fprintf(stderr, "error: invalid s-w vector "
		"hit threshold\n");
	exit(1);
      }
      if (strcspn(optarg, "%.") == strlen(optarg))
	sw_vect_threshold = -sw_vect_threshold;	//absol.
      break;
    case 'B':
      Bflag = true;
      break;
    case 'C':
      if (Fflag) {
	fprintf(stderr, "error: -C and -F are mutually "
		"exclusive\n");
	exit(1);
      }
      Cflag = true;
      break;
    case 'F':
      if (Cflag) {
	fprintf(stderr, "error: -C and -F are mutually "
		"exclusive\n");
	exit(1);
      }
      Fflag = true;
      break;
    case 'H':
      Hflag = true;
      break;
    case 'P':
      Pflag = true;
      break;
    case 'R':
      Rflag = true;
      break;
    case 'T':
      Tflag = true;
      break;
    case 'U':
      Uflag = true;
      break;
    /*
     * New options/parameters since SHRiMP 1.2.1
     */
    case 'A':
      anchor_width = atoi(optarg);
      if (anchor_width < -1 || anchor_width >= 100) {
	fprintf(stderr, "error: anchor_width requested is invalid (%s)\n",
		optarg);
	exit(1);
      }
      break;
    case 'D': // disabled
#ifdef EXTRA_STATS
      reads_profile_file_name = optarg;
#else
      fprintf(stderr, "warning: option -D is used only if -DEXTRA_STATS is defined at compile time; ignoring it\n");
#endif
      break;
    case 'Z':
      hash_filter_calls = !hash_filter_calls;
      break;
    case 'Y':
      c = strtok(optarg, ",");
      if (c == NULL) {
	fprintf(stderr, "error: format for cache sizes is \"-Y 2,16\"\n");
	exit(1);
      }
      cache_min_size = (uint)atoi(c);
      c = strtok(NULL, ",");
      if (c == NULL) {
	fprintf(stderr, "error: format for cache sizes is \"-Y 2,16\"\n");
	exit(1);
      }
      cache_max_size = (uint)atoi(c);
      if (cache_min_size > DEF_CACHE_MAX_SIZE || cache_max_size > DEF_CACHE_MAX_SIZE) {
	fprintf(stderr, "error: cache size too big (%d,%d)\n",
		cache_min_size, cache_max_size);
	exit(1);
      }
      break;
    case 'W':
      window_overlap = atof(optarg);
      if (window_overlap <= 0.0) {
	fprintf(stderr, "error: invalid window overlap\n");
	exit(1);
      }
      if (strcspn(optarg, "%.") == strlen(optarg))
	window_overlap = -window_overlap;		//absol.
      break;
    case '?':
      usage(progname, true);
      break;
    case 'M':
      assert(shrimp_mode == MODE_COLOUR_SPACE || shrimp_mode == MODE_LETTER_SPACE);

      c = strtok(optarg, ",");
      do {
	if (!set_mode_from_string(c)) {
	  fprintf(stderr, "error: unrecognized mode (%s)\n", c);
	  exit(1);
	}
	c = strtok(NULL, ",");
      } while (c != NULL);
      break;
    default:
      usage(progname, false);
    }
  }
  argc -= optind;
  argv += optind;

  /* Parameters from mode */
  if (shrimp_mode == MODE_COLOUR_SPACE || shrimp_mode == MODE_LETTER_SPACE) {
    set_params_from_mode(num_matches_set, window_len_set, hit_taboo_len_set);
  } else { // helicos
    if (n_seeds == 0)
      load_default_seeds(-1);
    anchor_width = -1;
  }

  /*
  if (!match_score_set)
    match_score = default_score_set[selected_score_set_id][0];
  if (!mismatch_score_set)
    mismatch_score = default_score_set[selected_score_set_id][1];
  if (!a_gap_open_score_set)
    a_gap_open_score = default_score_set[selected_score_set_id][2];
  if (!a_gap_extend_score_set)
    a_gap_extend_score = default_score_set[selected_score_set_id][3];
  if (!crossover_score_set)
    crossover_score = default_score_set[selected_score_set_id][4];
  if (!b_gap_open_score_set)
    b_gap_open_score = default_score_set[selected_score_set_id][5];
  if (!b_gap_extend_score_set)
    b_gap_extend_score = default_score_set[selected_score_set_id][6];
  */

  kmer_to_mapidx = kmer_to_mapidx_orig;
  if (Hflag)
    kmer_to_mapidx = kmer_to_mapidx_hash;

  if (Hflag)
    init_seed_hash_mask();

  if (Hflag && kmer_stddev_limit >= 0) {
    fprintf(stderr, "error: the -H flag may not be used with -d\n");
    exit(1);
  }

  if (argc < 2) {
    fprintf(stderr, "error: %sgenome_file(s) not specified\n",
	    (argc == 0) ? "reads_file, " : "");
    usage(progname, false);
  }

  reads_file    = argv[0];
  genome_files  = &argv[1];
  ngenome_files = argc - 1;

  if (!Cflag && !Fflag)
    Cflag = Fflag = true;

  if (shrimp_mode == MODE_LETTER_SPACE)
    sw_vect_threshold = sw_full_threshold;

  if (!valid_spaced_seeds()) {
    fprintf(stderr, "error: invalid spaced seed\n");
    if (!Hflag)
      fprintf(stderr, "       for longer seeds, try using the -H flag\n"); 
    exit(1);
  }

  if (!IS_ABSOLUTE(window_len) && window_len < 100.0) {
    fprintf(stderr, "error: window length < 100%% "
	    "of read length\n");
    exit(1);
  }

  if (num_matches < 1) {
    fprintf(stderr, "error: invalid number of matches\n");
    exit(1);
  }

  if (num_outputs < 1) {
    fprintf(stderr, "error: invalid maximum hits per read\n");
    exit(1);
  }

  if (a_gap_open_score > 0 || b_gap_open_score > 0) {
    fprintf(stderr, "error: invalid gap open penalty\n");
    exit(1);
  }

  if (a_gap_extend_score > 0 || b_gap_extend_score > 0) {
    fprintf(stderr, "error: invalid gap extend penalty\n");
    exit(1);
  }

  if (!IS_ABSOLUTE(sw_full_threshold) && sw_full_threshold > 100.0) {
    fprintf(stderr, "error: invalid s-w full hit threshold\n");
    exit(1);
  }

  if (shrimp_mode == MODE_COLOUR_SPACE && !IS_ABSOLUTE(sw_vect_threshold) &&
      sw_vect_threshold > 100.0) {
    fprintf(stderr, "error: invalid s-w full hit threshold\n");
    exit(1);
  }

  if ((a_gap_open_score_set && !b_gap_open_score_set)
      || (a_gap_extend_score_set && !b_gap_extend_score_set))
    fputc('\n', stderr);
  if (a_gap_open_score_set && !b_gap_open_score_set) {
    fprintf(stderr, "Notice: Gap open penalty set for reference but not query; assuming symmetry.\n");
    b_gap_open_score = a_gap_open_score;
  }
  if (a_gap_extend_score_set && !b_gap_extend_score_set) {
    fprintf(stderr, "Notice: Gap extend penalty set for reference but not query; assuming symmetry.\n");
    b_gap_extend_score = a_gap_extend_score;
  }
  if ((a_gap_open_score_set && !b_gap_open_score_set)
      || (a_gap_extend_score_set && !b_gap_extend_score_set))
    fputc('\n', stderr);

  fprintf(stderr, "Settings:\n");
  fprintf(stderr, "    Spaced Seed%s (weight/span)            %s%s (%u/%u)\n",
	  (n_seeds == 1) ? "" : "s", (n_seeds == 1) ? " " : "",
	  seed_to_string(0), seed[0].weight, seed[0].span);
  for (i = 1; i < (int)n_seeds; i++) {
    fprintf(stderr, "                                          %s (%u/%u)\n",
	    seed_to_string(i), seed[i].weight, seed[i].span);
  }
  fprintf(stderr, "    Seed Matches per Window:              %u\n",
	  num_matches);
  fprintf(stderr, "    Seed Hit Taboo Length:                %d%s\n",
	  hit_taboo_len, hit_taboo_len == -1 ? " (disabled)" : "");

  if (shrimp_mode != MODE_HELICOS_SPACE) {
    fprintf(stderr, "    Seed Generation Taboo Length:         %u\n",
	    seed_taboo_len);
  }

  if (IS_ABSOLUTE(window_len)) {
    fprintf(stderr, "    Seed Window Length:                   "
	    "%u\n", (u_int)-window_len);
  } else {
    fprintf(stderr, "    Seed Window Length:                   "
	    "%.02f%%\n", window_len);
  }

  if (IS_ABSOLUTE(window_overlap)) {
    fprintf(stderr, "    Seed Window Overlap Length:           %u\n",
	    (u_int)-window_overlap);
  } else {
    fprintf(stderr, "    Seed Window Overlap Length:           %.02f%%\n",
	    window_overlap);
  }

  fprintf(stderr, "    Maximum Hits per Read:                %u\n",
	  num_outputs);
  fprintf(stderr, "    Kmer Std. Deviation Limit:            %d%s\n",
	  kmer_stddev_limit, (kmer_stddev_limit < 0) ? " (None)" : "");

  if (shrimp_mode == MODE_HELICOS_SPACE) {
    fprintf(stderr, "\n");
    fprintf(stderr, "    DAG Epsilon:                          "
	    "%u\n", dag_epsilon);
    fprintf(stderr, "    DAG Read Match:                       "
	    "%d\n", dag_read_match);
    fprintf(stderr, "    DAG Read Mismatch:                    "
	    "%d\n", dag_read_mismatch);
    fprintf(stderr, "    DAG Read Gap:                         "
	    "%d\n", dag_read_gap);
    fprintf(stderr, "    DAG Reference Match:                  "
	    "%d\n", dag_ref_match);
    fprintf(stderr, "    DAG Reference Mismatch:               "
	    "%d\n", dag_ref_mismatch);
    fprintf(stderr, "    DAG Reference Half Match:             "
	    "%d\n", dag_ref_half_match);
    fprintf(stderr, "    DAG Reference Neither Match:          "
	    "%d\n", dag_ref_neither_match);
    fprintf(stderr, "    DAG Reference Match,Deletion:         "
	    "%d\n", dag_ref_match_deletion);
    fprintf(stderr, "    DAG Reference Mismatch,Deletion:      "
	    "%d\n", dag_ref_mismatch_deletion);
    fprintf(stderr, "    DAG Reference Error,Insertion:        "
	    "%d\n", dag_ref_error_insertion);
    fprintf(stderr, "    DAG Reference Weighted Threshold:     "
	    "%.03f\n", dag_ref_weighted_thresh);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "    S-W Match Value:                      "
	  "%d\n", match_score);
  fprintf(stderr, "    S-W Mismatch Value:                   "
	  "%d\n", mismatch_score);
  fprintf(stderr, "    S-W Gap Open Penalty (Ref):           "
	  "%d\n", a_gap_open_score);
  fprintf(stderr, "    S-W Gap Open Penalty (Qry):           "
	  "%d\n", b_gap_open_score);
  fprintf(stderr, "    S-W Gap Extend Penalty (Ref):         "
	  "%d\n", a_gap_extend_score);
  fprintf(stderr, "    S-W Gap Extend Penalty (Qry):         "
	  "%d\n", b_gap_extend_score);

  if (shrimp_mode == MODE_COLOUR_SPACE) {
    fprintf(stderr, "    S-W Crossover Penalty:                "
	    "%d\n", crossover_score);
  }

  if (shrimp_mode == MODE_COLOUR_SPACE || shrimp_mode == MODE_HELICOS_SPACE) {
    if (IS_ABSOLUTE(sw_vect_threshold)) {
      fprintf(stderr, "    S-W Vector Hit Threshold:     "
	      "        %u\n", (u_int)-sw_vect_threshold);
    } else {
      fprintf(stderr, "    S-W Vector Hit Threshold:     "
	      "        %.02f%%\n", sw_vect_threshold);
    }
  }
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    if (IS_ABSOLUTE(sw_full_threshold)) {
      fprintf(stderr, "    S-W Full Hit Threshold:       "
	      "        %u\n", (u_int)-sw_full_threshold);
    } else {
      fprintf(stderr, "    S-W Full Hit Threshold:       "
	      "        %.02f%%\n", sw_full_threshold);
    }
  }
  if (shrimp_mode == MODE_LETTER_SPACE) {
    if (IS_ABSOLUTE(sw_full_threshold)) {
      fprintf(stderr, "    S-W Hit Threshold:            "
	      "        %u\n", (u_int)-sw_full_threshold);
    } else {
      fprintf(stderr, "    S-W Hit Threshold:            "
	      "        %.02f%%\n", sw_full_threshold);
    }
  }

  fprintf(stderr, "\n");
  if (shrimp_mode == MODE_COLOUR_SPACE || shrimp_mode == MODE_LETTER_SPACE) {
    fprintf(stderr, "    Anchor width:                         %d%s\n",
	    anchor_width, anchor_width == -1 ? " (disabled)" : "");
  }
  fprintf(stderr, "    Hash filter SW calls:                 %s\n",
	  hash_filter_calls ? "yes" : "no");
  fprintf(stderr, "    Cache settings:                       (%d,%d)\n",
	  cache_min_size, cache_max_size);
  fprintf(stderr, "\n");

  dag_setup(dag_read_match, dag_read_mismatch, dag_read_gap, dag_ref_match,
	    dag_ref_mismatch, dag_ref_half_match, dag_ref_neither_match,
	    dag_ref_match_deletion, dag_ref_mismatch_deletion,
	    dag_ref_error_insertion);

  if (!load_reads(reads_file))
    exit(1);

  fprintf(stderr, "  - Configuring window lengths...\n");
  max_window_len = set_window_lengths();
  fprintf(stderr, "  - Maximum window length: %u\n", max_window_len);
		

  if (sw_vector_setup(max_window_len, longest_read_len, a_gap_open_score,
		      a_gap_extend_score, b_gap_open_score, b_gap_extend_score, match_score, mismatch_score,
		      shrimp_mode == MODE_COLOUR_SPACE, false)) {
    fprintf(stderr, "failed to initialise vector "
	    "Smith-Waterman (%s)\n", strerror(errno));
    exit(1);
  }

  if (shrimp_mode == MODE_COLOUR_SPACE) {
    /* XXX - a vs. b gap */
    ret = sw_full_cs_setup(max_window_len, longest_read_len,
			   a_gap_open_score, a_gap_extend_score, match_score, mismatch_score,
			   crossover_score, false, anchor_width);
  } else {
    ret = sw_full_ls_setup(max_window_len, longest_read_len,
			   a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,
			   match_score, mismatch_score, false, anchor_width);
  }
  if (ret) {
    fprintf(stderr, "failed to initialise scalar "
	    "Smith-Waterman (%s)\n", strerror(errno));
    exit(1);
  }

  readmap_prune();

  if (scan_genomes(genome_files, ngenome_files))
    exit(1);

  fprintf(stderr, "\nGenerating output...\n");
  final_pass(genome_files, ngenome_files);

  if (Uflag) {
    unsigned int j;

    fprintf(stderr, "\nListing unmapped reads...\n");
    printf("#\n#UNMAPPED READS:\n#\n");

    for (j = 0; j < nreads; j++) {
      if (reads[j].final_matches == 0)
	printf("%s\n", reads[j].name);
    }
  }

  print_statistics();

  return (0);
}
