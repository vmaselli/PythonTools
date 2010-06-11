#ifndef _ANCHORS_H
#define _ANCHORS_H

#include <stdint.h>
#include <sys/types.h>

struct anchor {
  int16_t	x;
  int16_t	y;
  uint8_t	y_alt;
  uint8_t	length;
  uint8_t	width;
  uint8_t	more_than_once;
};


void join_anchors(struct anchor *, uint, struct anchor *);
void widen_anchor(struct anchor *, uint);
void get_x_range(struct anchor *, uint, uint, int, int *, int *);


#endif
