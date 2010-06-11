#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <stdbool.h>
#include "anchors.h"

#ifdef DEBUG_ANCHORS
#include <stdio.h>
#endif


void join_anchors(struct anchor * anchors, uint anchors_cnt,
		  struct anchor * dest) {
  int border_nw_min, border_sw_min, border_ne_max, border_se_max;
  uint i;
  bool once_more;

  assert(anchors != NULL && dest != NULL);

  border_nw_min = INT_MAX;
  border_sw_min = INT_MAX;
  border_ne_max = INT_MIN;
  border_se_max = INT_MIN;

  for (i = 0, once_more = false; i < anchors_cnt; i += (once_more? 0 : 1)) {
    int border_nw, border_sw, border_ne, border_se, y;

    if (!once_more) {
      y = anchors[i].y;
    } else {
      y = anchors[i].y_alt;
    }

    border_nw = anchors[i].x + y;
    border_sw = anchors[i].x - y;
    border_ne = border_sw + 2*(anchors[i].width - 1);
    border_se = border_nw + 2*(anchors[i].length - 1);

    if (border_nw < border_nw_min) border_nw_min = border_nw;
    if (border_sw < border_sw_min) border_sw_min = border_sw;
    if (border_ne > border_ne_max) border_ne_max = border_ne;
    if (border_se > border_se_max) border_se_max = border_se;

#ifdef DEBUG_ANCHORS
    fprintf(stderr, "i:%d once_more:%s border_nw_min=%d border_sw_min=%d border_ne_max=%d border_se_max=%d\n",
	    i, (once_more? "true" : "false"), border_nw_min, border_sw_min, border_ne_max, border_se_max);
#endif

    if (once_more) {
      once_more = false;
    } else {
      if (anchors[i].more_than_once != 0) {
	once_more = true;
      }
    }
  }

  if ((border_nw_min + border_sw_min) % 2 != 0) border_nw_min--;
  dest->x = (border_nw_min + border_sw_min)/2;
  dest->y = border_nw_min - dest->x;

  if ((border_ne_max - border_sw_min) % 2 != 0) border_ne_max++;
  dest->width = (border_ne_max - border_sw_min)/2 + 1;

  if ((border_se_max - border_nw_min) % 2 != 0) border_se_max++;
  dest->length = (border_se_max - border_nw_min)/2 + 1;

  dest->more_than_once = 0;
}


void widen_anchor(struct anchor * anchor, uint width) {
  assert(anchor != NULL);

  anchor->x -= width/2;
  anchor->y += width/2;
  anchor->width += width;
}


void get_x_range(struct anchor * anchor, uint x_len, uint y_len, int y,
		 int * x_min, int * x_max) {
  assert(anchor != NULL && x_min != NULL && x_max != NULL);

  if (y < anchor->y) {
    *x_min = 0;
  } else if (y <= anchor->y + (anchor->length - 1)) {
    *x_min = anchor->x + (y - anchor->y);
  } else {
    *x_min = anchor->x + anchor->length;
  }

  if (*x_min < 0) *x_min = 0;
  if ((uint)*x_min >= x_len) *x_min = x_len - 1;

  if (y < anchor->y - (anchor->width - 1)) {
    *x_max = anchor->x + (anchor->width - 1) - 1;
  } else if (y <= anchor->y - (anchor->width - 1) + (anchor->length - 1)) {
    *x_max = anchor->x + (anchor->width - 1) + (y - (anchor->y - (anchor->width - 1)));
  } else {
    *x_max = x_len - 1;
  }
  
  if (*x_max < 0) *x_max = 0;
  if ((uint)*x_max >= x_len) *x_max = x_len - 1;
}
