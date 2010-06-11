#ifndef _HASH_H
#define _HASH_H

#include <stdint.h>
#include <assert.h>

#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((uint32_t)(((const uint8_t *)(d))[1])) << 8)\
                       +(uint32_t)(((const uint8_t *)(d))[0]) )
#endif


static inline uint32_t
SuperFastHash (const char * data, int len, uint32_t hash) {
  uint32_t tmp;
  int rem;

  assert(data != NULL);

  rem = len & 3;
  len >>= 2;

  /* Main loop */
  for ( ; len > 0; len--) {
    hash  += get16bits (data);
    tmp    = (get16bits (data+2) << 11) ^ hash;
    hash   = (hash << 16) ^ tmp;
    data  += 2*sizeof (uint16_t);
    hash  += hash >> 11;
  }

  /* Handle end cases */
  switch (rem) {
  case 3:
    hash += get16bits (data);
    hash ^= hash << 16;
    hash ^= data[sizeof (uint16_t)] << 18;
    hash += hash >> 11;
    break;
  case 2:
    hash += get16bits (data);
    hash ^= hash << 11;
    hash += hash >> 17;
    break;
  case 1:
    hash += *data;
    hash ^= hash << 10;
    hash += hash >> 1;
  }

  /* Force "avalanching" of final 127 bits */
  hash ^= hash << 3;
  hash += hash >> 5;
  hash ^= hash << 4;
  hash += hash >> 17;
  hash ^= hash << 25;
  hash += hash >> 6;

  return hash;
}


static inline void
hash_accumulate(uint32_t * key, uint32_t val) {
  assert(key != NULL);

  uint32_t tmp;

  *key += (val >> 16); // high 16 bits
  tmp = ((val & 0xFFFF) << 11) ^ *key; // high 16 bits
  *key = (*key << 16) ^ tmp;
  *key += *key >> 11;
}


static inline void
hash_finalize(uint32_t * key) {
  assert(key != NULL);

  *key ^= *key << 3;
  *key += *key >> 5;
  *key ^= *key << 4;
  *key += *key >> 17;
  *key ^= *key << 25;
  *key += *key >> 6;
}


#endif
