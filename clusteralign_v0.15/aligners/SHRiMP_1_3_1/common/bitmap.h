#ifndef _BITMAP_H
#define _BITMAP_H

#include <assert.h>
#include <inttypes.h>
#include <sys/types.h>
#include <stdbool.h>


/*
 * Public declarations
 */

typedef uint64_t bitmap_type;
//#define bitmap_type uint64_t

static inline void bitmap_insert_clean(bitmap_type * a, uint bits_per_field, uint index, uint value);
static inline void bitmap_clear(bitmap_type * a, uint bits_per_field, uint index);
static inline void bitmap_prepend(bitmap_type * a, uint bits_per_field, uint value);
static inline uint bitmap_extract(bitmap_type * a, uint bits_per_field, uint index);

static inline uint bitmap32_extract(uint32_t * a, uint bits_per_field, uint index);

char const * bitmap32v_string(uint32_t * a, uint n_fields, uint bits_per_field,
			      bool reverse, bool use_bits, bool use_commas);

/*
 * Inline definitions
 */

static inline bitmap_type
bitmap_all1_field(uint bits_per_field) {
  return ((bitmap_type)1 << bits_per_field) - 1;
}


static inline uint32_t
bitmap32_all1_field(uint bits_per_field) {
  return ((uint32_t)1 << bits_per_field) - 1;
}


static inline void
bitmap_insert_clean(bitmap_type * a, uint bits_per_field, uint index, uint value) {
  assert(index < (8 * sizeof(bitmap_type)) / bits_per_field);
  assert((value & ~bitmap_all1_field(bits_per_field)) == 0);

  *a |= ((bitmap_type)value << index * bits_per_field);
}


static inline void
bitmap_clear(bitmap_type * a, uint bits_per_field, uint index) {
  assert(index < (8 * sizeof(bitmap_type)) / bits_per_field);

  *a &= ~(bitmap_all1_field(bits_per_field) << index * bits_per_field);
}


static inline void
bitmap_prepend(bitmap_type * a, uint bits_per_field, uint value) {
  assert((value & ~bitmap_all1_field(bits_per_field)) == 0);

  *a <<= bits_per_field;
  *a |= value;
}


static inline uint
bitmap_extract(bitmap_type * a, uint bits_per_field, uint index) {
  assert(index < (8 * sizeof(bitmap_type)) / bits_per_field);

  return (*a >> bits_per_field * index) & bitmap_all1_field(bits_per_field);
}


static inline uint
bitmap32_extract(uint32_t * a, uint bits_per_field, uint index) {
  assert(index < 32/bits_per_field);

  return (*a >> bits_per_field * index) & bitmap32_all1_field(bits_per_field);
}


#endif
