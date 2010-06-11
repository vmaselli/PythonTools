#include <stdio.h>
#include <assert.h>
#include "bitmap.h"


char const *
bitmap32v_string(uint32_t *a, uint n_fields, uint bits_per_field,
		 bool reverse, bool use_bits, bool use_commas) {
  static char buffer[1001];

  int i, j, k;
  uint32_t val;
  uint fields_per_byte = 32 / bits_per_field;

  int first, last, increment;

  if (!reverse) {
    first = n_fields - 1;
    last = 0;
    increment = -1;
  } else {
    first = 0;
    last = n_fields - 1;
    increment = 1;
  }

  for (i = first, j = 0;
       i != last + increment && j < 990;
       i += increment) {
    val = ((a[i / fields_per_byte] >> bits_per_field * (i % fields_per_byte))
	   & bitmap32_all1_field(bits_per_field));
    if (use_bits) {
      for (k = bits_per_field - 1; k >= 0; k--) {
	buffer[j++] = (bitmap32_extract(&val, 1, k) == 1? '1' : '0');
      }
    } else {
      j += sprintf(&buffer[j], "%u", val); // unsafe
    }
    if (use_commas && i != last)
      buffer[j++] = ',';
  }
  buffer[j] = 0;

  return buffer;
}
