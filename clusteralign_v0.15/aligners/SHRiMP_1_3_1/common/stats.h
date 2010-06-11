#ifndef _STATS_H
#define _STATS_H

#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>
#include <math.h>

#include "../common/util.h"


/*
 * Statistics for which we can measure mean, std dev.
 */

typedef struct {
  int64_t	sum;
  int64_t	sum_sq;
  int64_t	count;
} stat_t;

/*
 * Methods
 */

/* Constructor */
static inline stat_t *
stat_init(stat_t * a) {
  if (a == NULL)
    a = (stat_t *)calloc(sizeof(stat_t), 1);
  else {
    a->sum = 0;
    a->sum_sq = 0;
    a->count = 0;
  }

  return a;
}

/* Add measurement */
static inline void
stat_add(stat_t * a, int64_t val) {
  assert(a != NULL);

  a->sum += val;
  a->sum_sq += val*val;
  (a->count)++;
}

/* Get sum */
static inline int64_t
stat_get_sum(stat_t * a) {
  assert(a != NULL);

  return a->sum;
}

/* Get count */
static inline int64_t
stat_get_count(stat_t * a) {
  assert(a != NULL);

  return a->count;
}

/* Get mean */
static inline double
stat_get_mean(stat_t * a) {
  assert(a != NULL);

  return (double)a->sum / (double)a->count;
}

/* Get sample std dev */
static inline double
stat_get_sample_stddev(stat_t * a) {
  assert(a != NULL);

  return sqrt( (double)a->sum_sq / ((double)a->count - 1)
	       - ((double)a->sum)*((double)a->sum) / (((double)a->count) * ((double)a->count - 1)) );
}


/*
 * Simple count
 */

typedef int64_t count_t;

static inline count_t *
count_init(count_t * c) {
  if (c == NULL)
    c = (count_t *)calloc(sizeof(count_t), 1);
  else
    *c = 0;

  return c;
}

static inline void
count_increment(count_t * c) {
  assert(c != NULL);

  (*c)++;
}

static inline void
count_add(count_t * c, int64_t val) {
  assert(c != NULL);

  *c += val;
}

static inline int64_t
count_get_count(count_t * c) {
  assert(c != NULL);

  return *c;
}


typedef int32_t count32_t;

static inline count32_t *
count32_init(count32_t * c) {
  if (c == NULL)
    c = (count32_t *)calloc(sizeof(count32_t), 1);
  else
    *c = 0;

  return c;
}

static inline void
count32_increment(count32_t * c) {
  assert(c != NULL);

  (*c)++;
}

static inline void
count32_add(count32_t * c, int32_t val) {
  assert(c != NULL);

  *c += val;
}

static inline int32_t
count32_get_count(count32_t * c) {
  assert(c != NULL);

  return *c;
}


#endif
