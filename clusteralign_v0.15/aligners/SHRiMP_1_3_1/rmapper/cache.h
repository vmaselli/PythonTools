#ifndef _CACHE_H
#define _CACHE_H


#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>
#include <string.h>

#include "../common/util.h"
#include "../common/stats.h"

struct read_entry; // in rmapper.h


typedef uint32_t cache_key_t;
static uint cache_min_size = 2;
static uint cache_max_size = 32;

#define DEF_CACHE_MAX_SIZE 127 /* cannot be set to more than this value */

struct cache_entry {
  cache_key_t	key;
  uint16_t	score;		/* scores <= 2^15-1 (cf., sw-vector.c) */
  uint16_t	count;
};

/*
 * Cache manipulation routines.
 */

/* Initialize cache */
static inline void
cache_init(struct read_entry * re, count_t * c) {
  assert(re != NULL && re->cache == NULL && re->cache_sz == 0);

  re->cache = (struct cache_entry *)
    xcalloc_c(cache_min_size * sizeof(struct cache_entry),
	      c);
  re->cache_sz = cache_min_size;
}

/* Possibly double cache size; moving the tail of the queue */
static inline void
cache_maybe_resize(struct read_entry * re, count_t * c) {
  assert(re != NULL && re->cache != NULL && re->cache_sz > 0);

  if (re->cache_sz < cache_max_size
      //&& (re->cache_sz < 16 || 3 * re->filter_calls_bypassed > re->filter_calls)
      // only get 32 if more than .33 calls bypassed
      ) {
    re->cache = (struct cache_entry *)
      xrealloc_c(re->cache, 2 * re->cache_sz * sizeof(struct cache_entry),
		 re->cache_sz * sizeof(struct cache_entry),
		 c);
    if (re->head != re->cache_sz - 1) {
      memcpy(&re->cache[(re->head + 1) + re->cache_sz],
	     &re->cache[re->head + 1],
	     ((re->cache_sz - 1) - re->head) * sizeof(struct cache_entry));
    }
    memset(&re->cache[re->head + 1], 0,
	   re->cache_sz * sizeof(struct cache_entry));
    re->cache_sz *= 2;
  }
}

/* Lookup key */
static inline void
cache_lookup(struct read_entry * re, cache_key_t key, uint * rel_pos, uint * abs_pos) {
  assert(re != NULL && re->cache != NULL && rel_pos != NULL && abs_pos != NULL);
  assert(re->head < re->cache_sz);

  (*abs_pos) = re->head;
  (*rel_pos) = 0;
  while (*rel_pos < re->cache_sz
	 && re->cache[*abs_pos].count != 0
	 && re->cache[*abs_pos].key != key) {
    (*rel_pos)++;
    (*abs_pos) = ((*abs_pos) + re->cache_sz - 1) % re->cache_sz; // -= 1 mod ..
  }
  /*
   * Outcomes:
   *   *rel_pos == re->cache_sz:
   *     key not found; cache full
   *   *rel_pos < re->cache_sz && re->cache[*abs_pos].count == 0:
   *     key not found; cache not full
   *   *rel_pos < re->cache_sz && re->cache[*abs_pos].count != 0:
   *     key found at *abs_pos
   */
}

/* Predicates valid only after lookup call */
static inline bool
cache_key_found(struct read_entry * re, uint * rel_pos, uint * abs_pos) {
  assert (re != NULL && re->cache != NULL && rel_pos != NULL && abs_pos != NULL);

  return *rel_pos < re->cache_sz && re->cache[*abs_pos].count != 0;
}
static inline bool
cache_is_full(struct read_entry * re, uint * rel_pos, uint * abs_pos) {
  assert (re != NULL && re->cache != NULL && rel_pos != NULL && abs_pos != NULL);

  return *rel_pos == re->cache_sz;
}

/* Add key */
static inline void
cache_add_key(struct read_entry * re, uint * rel_pos, uint * abs_pos) {
  assert (re != NULL && re->cache != NULL && rel_pos != NULL && abs_pos != NULL);

  *abs_pos = (re->head + 1) % re->cache_sz;
  re->head = *abs_pos;
  /*
   * The caller must fill in the fields of freq_hits[*abs_pos]
   */
}


#endif
