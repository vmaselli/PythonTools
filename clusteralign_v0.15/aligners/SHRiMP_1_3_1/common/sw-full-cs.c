/*	$Id: sw-full-cs.c 398 2009-10-06 04:20:35Z matei $	*/

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "../common/fasta.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/util.h"

struct swcell {
  struct {
    int	score_n;
    int	score_w;
    int	score_nw;

    int8_t	back_n;
    int8_t	back_w;
    int8_t	back_nw;
  } from[4];
};

#define FROM_A	0x00
#define FROM_B	0x01 
#define FROM_C	0x02
#define FROM_D	0x03

#define FROM_NORTH_NORTH		0x01
#define FROM_NORTH_NORTHWEST		0x02
#define FROM_WEST_NORTHWEST		0x03
#define FROM_WEST_WEST			0x04
#define FROM_NORTHWEST_NORTH		0x05
#define FROM_NORTHWEST_NORTHWEST	0x06
#define FROM_NORTHWEST_WEST		0x07

#define FROM_x(_mat, _dir)		(int8_t)(((_dir) << 2) | (_mat))

#define FROM_A_NORTH_NORTH		FROM_x(FROM_A, FROM_NORTH_NORTH)
#define FROM_A_NORTH_NORTHWEST		FROM_x(FROM_A, FROM_NORTH_NORTHWEST)
#define	FROM_A_WEST_NORTHWEST		FROM_x(FROM_A, FROM_WEST_NORTHWEST)
#define FROM_A_WEST_WEST		FROM_x(FROM_A, FROM_WEST_WEST)
#define FROM_A_NORTHWEST_NORTH		FROM_x(FROM_A, FROM_NORTHWEST_NORTH)
#define FROM_A_NORTHWEST_NORTHWEST	FROM_x(FROM_A, FROM_NORTHWEST_NORTHWEST)
#define FROM_A_NORTHWEST_WEST		FROM_x(FROM_A, FROM_NORTHWEST_WEST)

#define FROM_B_NORTH_NORTH		FROM_x(FROM_B, FROM_NORTH_NORTH)
#define FROM_B_NORTH_NORTHWEST		FROM_x(FROM_B, FROM_NORTH_NORTHWEST)
#define	FROM_B_WEST_NORTHWEST		FROM_x(FROM_B, FROM_WEST_NORTHWEST)
#define FROM_B_WEST_WEST		FROM_x(FROM_B, FROM_WEST_WEST)
#define FROM_B_NORTHWEST_NORTH		FROM_x(FROM_B, FROM_NORTHWEST_NORTH)
#define FROM_B_NORTHWEST_NORTHWEST	FROM_x(FROM_B, FROM_NORTHWEST_NORTHWEST)
#define FROM_B_NORTHWEST_WEST		FROM_x(FROM_B, FROM_NORTHWEST_WEST)

#define FROM_C_NORTH_NORTH		FROM_x(FROM_C, FROM_NORTH_NORTH)
#define FROM_C_NORTH_NORTHWEST		FROM_x(FROM_C, FROM_NORTH_NORTHWEST)
#define	FROM_C_WEST_NORTHWEST		FROM_x(FROM_C, FROM_WEST_NORTHWEST)
#define FROM_C_WEST_WEST		FROM_x(FROM_C, FROM_WEST_WEST)
#define FROM_C_NORTHWEST_NORTH		FROM_x(FROM_C, FROM_NORTHWEST_NORTH)
#define FROM_C_NORTHWEST_NORTHWEST	FROM_x(FROM_C, FROM_NORTHWEST_NORTHWEST)
#define FROM_C_NORTHWEST_WEST		FROM_x(FROM_C, FROM_NORTHWEST_WEST)

#define FROM_D_NORTH_NORTH		FROM_x(FROM_D, FROM_NORTH_NORTH)
#define FROM_D_NORTH_NORTHWEST		FROM_x(FROM_D, FROM_NORTH_NORTHWEST)
#define	FROM_D_WEST_NORTHWEST		FROM_x(FROM_D, FROM_WEST_NORTHWEST)
#define FROM_D_WEST_WEST		FROM_x(FROM_D, FROM_WEST_WEST)
#define FROM_D_NORTHWEST_NORTH		FROM_x(FROM_D, FROM_NORTHWEST_NORTH)
#define FROM_D_NORTHWEST_NORTHWEST	FROM_x(FROM_D, FROM_NORTHWEST_NORTHWEST)
#define FROM_D_NORTHWEST_WEST		FROM_x(FROM_D, FROM_NORTHWEST_WEST)

enum {
  BACK_INSERTION = 1,
  BACK_A_DELETION,
  BACK_B_DELETION,
  BACK_C_DELETION,
  BACK_D_DELETION,
  BACK_A_MATCH_MISMATCH,
  BACK_B_MATCH_MISMATCH,
  BACK_C_MATCH_MISMATCH,
  BACK_D_MATCH_MISMATCH
};

static int		initialised;
static int8_t	       *db, *qr[4];
static int		dblen, qrlen;
static int		gap_open, gap_ext;
static int		match, mismatch;
static int		xover_penalty;
static struct swcell   *swmatrix;
static uint8_t	       *backtrace;
static char	       *dbalign, *qralign;
static int		anchor_width;

/* statistics */
static uint64_t		swticks, swcells, swinvocs;

#define BT_CROSSOVER		0x80
#define BT_ISCROSSOVER(_x)	((_x) & BT_CROSSOVER)
#define BT_TYPE(_x)		((_x) & 0x0f)

inline static void
init_cell(int idx) {
  swmatrix[idx].from[0].score_nw = 0;
  swmatrix[idx].from[0].score_n  = -gap_open;
  swmatrix[idx].from[0].score_w  = -gap_open;
  swmatrix[idx].from[1].score_nw = xover_penalty;
  swmatrix[idx].from[1].score_n  = -gap_open + xover_penalty;
  swmatrix[idx].from[1].score_w  = -gap_open + xover_penalty;
  swmatrix[idx].from[2].score_nw = xover_penalty;
  swmatrix[idx].from[2].score_n  = -gap_open + xover_penalty;
  swmatrix[idx].from[2].score_w  = -gap_open + xover_penalty;
  swmatrix[idx].from[3].score_nw = xover_penalty;
  swmatrix[idx].from[3].score_n  = -gap_open + xover_penalty;
  swmatrix[idx].from[3].score_w  = -gap_open + xover_penalty;

  swmatrix[idx].from[0].back_nw = 0;
  swmatrix[idx].from[0].back_n  = 0;
  swmatrix[idx].from[0].back_w  = 0;
  swmatrix[idx].from[1].back_nw = 0;
  swmatrix[idx].from[1].back_n  = 0;
  swmatrix[idx].from[1].back_w  = 0;
  swmatrix[idx].from[2].back_nw = 0;
  swmatrix[idx].from[2].back_n  = 0;
  swmatrix[idx].from[2].back_w  = 0;
  swmatrix[idx].from[3].back_nw = 0;
  swmatrix[idx].from[3].back_n  = 0;
  swmatrix[idx].from[3].back_w  = 0;
}

/*
 * Perform a full Smith-Waterman alignment. For the colour case, this means
 * computing each possible letter space read string and doing a four layer
 * scan.
 */
static int
full_sw(int lena, int lenb, int threshscore, int *iret, int *jret,
	int *kret, bool revcmpl,
	struct anchor * anchors, uint anchors_cnt)
{
  int i, j, k, l, max_i, max_j, max_k;
  int score, ms, go, ge, tmp, resetval;
  //int sw_band, ne_band;
  int8_t tmp2;
  struct anchor rectangle;

  /* shut up gcc */
  max_i = max_j = max_k = j = 0;

  score = 0;
  go = gap_open;
  ge = gap_ext;

  for (j = 0; j < lena + 1; j++) {
    init_cell(j);
  }

  //for (j = 0; j < lenb + 1; j++) {
  //init_cell(j * (lena + 1));
  //}

  /*
   * Figure out our band.
   *   We can actually skip computation of a significant number of
   *   cells, which could never be part of an alignment corresponding
   *   to our threshhold score.
   */
  //sw_band = ((lenb * match - threshscore + match - 1) / match) + 1;
  //ne_band = lena - (lenb - sw_band);

  if (anchors != NULL && anchor_width >= 0) {
    join_anchors(anchors, anchors_cnt, &rectangle);
    widen_anchor(&rectangle, (uint)anchor_width);
  } else {
    struct anchor tmp_anchors[2];

    tmp_anchors[0].x = 0;
    tmp_anchors[0].y = (int16_t)((lenb * match - threshscore) / match);
    tmp_anchors[0].length = 1;
    tmp_anchors[0].width = 1;
    tmp_anchors[0].more_than_once = 0;

    tmp_anchors[1].x = (int16_t)(lena - 1);
    tmp_anchors[1].y = (int16_t)(lenb - 1 - tmp_anchors[0].y);
    tmp_anchors[1].length = 1;
    tmp_anchors[1].width = 1;
    tmp_anchors[1].more_than_once = 0;

    join_anchors(tmp_anchors, 2, &rectangle);
  }

  for (i = 0; i < lenb; i++) {
    /*
     * computing row i of virtual matrix, stored in row i+1
     */
    int x_min, x_max;

    get_x_range(&rectangle, lena, lenb, i, &x_min, &x_max);
    //if (x_min > 0) {
    init_cell((i + 1) * (lena + 1) + (x_min - 1) + 1);
    //}

    swcells += x_max - x_min + 1;

    for (j = x_min; j <= x_max; j++) {
      /*
       * computing column j of virtual matrix, stored in column j+1
       */
      struct swcell *cell_nw, *cell_n, *cell_w, *cell_cur;

      cell_nw  = &swmatrix[i * (lena + 1) + j];
      cell_n   = cell_nw + 1;
      cell_w   = cell_nw + (lena + 1);
      cell_cur = cell_w + 1;

      /* banding */
      //if (i >= sw_band + j) {
      //memset(cell_cur, 0, sizeof(*cell_cur));
      //continue;
      //}
      //if (j >= ne_band + i) {
      //memset(cell_cur, 0, sizeof(*cell_cur));
      //break;
      //}

      for (k = 0; k < 4; k++) {
	if (k != 0)
	  resetval = xover_penalty;
	else
	  resetval = 0;

	/*
	 * northwest
	 */
	ms = (db[j] == qr[k][i]) ? match : mismatch;

	if (!revcmpl) {
	  tmp  = cell_nw->from[k].score_nw + ms;
	  tmp2 = FROM_x(k, FROM_NORTHWEST_NORTHWEST);

	  if (cell_nw->from[k].score_n + ms > tmp) {
	    tmp  = cell_nw->from[k].score_n + ms;
	    tmp2 = FROM_x(k, FROM_NORTHWEST_NORTH);
	  }

	  if (cell_nw->from[k].score_w + ms > tmp) {
	    tmp  = cell_nw->from[k].score_w + ms;
	    tmp2 = FROM_x(k, FROM_NORTHWEST_WEST);
	  }
	} else {
	  tmp  = cell_nw->from[k].score_w + ms;
	  tmp2 = FROM_x(k, FROM_NORTHWEST_WEST);

	  if (cell_nw->from[k].score_n + ms > tmp) {
	    tmp  = cell_nw->from[k].score_n + ms;
	    tmp2 = FROM_x(k, FROM_NORTHWEST_NORTH);
	  }

	  if (cell_nw->from[k].score_nw + ms > tmp) {
	    tmp  = cell_nw->from[k].score_nw + ms;
	    tmp2 = FROM_x(k, FROM_NORTHWEST_NORTHWEST);
	  }
	}

	/* check neighbours */
	for (l = 0; l < 4; l++) {
	  if (l == k)
	    continue;

	  if (!revcmpl) {
	    /* northwest */
	    if (cell_nw->from[l].score_nw + ms + xover_penalty > tmp) {
	      tmp  = cell_nw->from[l].score_nw + ms + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTHWEST_NORTHWEST);
	    }

	    /* north */
	    if (cell_nw->from[l].score_n + ms + xover_penalty > tmp) {
	      tmp  = cell_nw->from[l].score_n + ms + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTHWEST_NORTH);
	    }

	    /* west */
	    if (cell_nw->from[l].score_w + ms + xover_penalty > tmp) {
	      tmp  = cell_nw->from[l].score_w + ms + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTHWEST_WEST);
	    }
	  } else {
	    /* west */
	    if (cell_nw->from[l].score_w + ms + xover_penalty > tmp) {
	      tmp  = cell_nw->from[l].score_w + ms + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTHWEST_WEST);
	    }

	    /* north */
	    if (cell_nw->from[l].score_n + ms + xover_penalty > tmp) {
	      tmp  = cell_nw->from[l].score_n + ms + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTHWEST_NORTH);
	    }

	    /* northwest */
	    if (cell_nw->from[l].score_nw + ms + xover_penalty > tmp) {
	      tmp  = cell_nw->from[l].score_nw + ms + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTHWEST_NORTHWEST);
	    }
	  }
	}

	if (tmp <= resetval) {
	  tmp = resetval;
	  tmp2 = 0;
	}

	cell_cur->from[k].score_nw = tmp;
	cell_cur->from[k].back_nw  = tmp2;


	/*
	 * north
	 */
	if (!revcmpl) {
	  tmp  = cell_n->from[k].score_nw - go - ge;
	  tmp2 = FROM_x(k, FROM_NORTH_NORTHWEST);

	  if (cell_n->from[k].score_n - ge > tmp) {
	    tmp  = cell_n->from[k].score_n - ge;
	    tmp2 = FROM_x(k, FROM_NORTH_NORTH);
	  }
	} else {
	  tmp  = cell_n->from[k].score_n - ge;
	  tmp2 = FROM_x(k, FROM_NORTH_NORTH);

	  if (cell_n->from[k].score_nw - go - ge > tmp) {
	    tmp  = cell_n->from[k].score_nw - go - ge;
	    tmp2 = FROM_x(k, FROM_NORTH_NORTHWEST);
	  }
	}

	/* check neighbours */
	for (l = 0; l < 4; l++) {
	  if (l == k)
	    continue;

	  if (!revcmpl) {
	    /* northwest */
	    if (cell_n->from[l].score_nw - go - ge + xover_penalty > tmp) {
	      tmp  = cell_n->from[l].score_nw - go - ge + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTH_NORTHWEST);
	    }

	    /* north */
	    if (cell_n->from[l].score_n - ge + xover_penalty > tmp) {
	      tmp  = cell_n->from[l].score_n - ge + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTH_NORTH);
	    }
	  } else {
	    /* north */
	    if (cell_n->from[l].score_n - ge + xover_penalty > tmp) {
	      tmp  = cell_n->from[l].score_n - ge + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTH_NORTH);
	    }

	    /* northwest */
	    if (cell_n->from[l].score_nw - go - ge + xover_penalty > tmp) {
	      tmp  = cell_n->from[l].score_nw - go - ge + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTH_NORTHWEST);
	    }
	  }
	}

	if (tmp <= resetval) {
	  tmp = resetval;
	  tmp2 = 0;
	}
					
	cell_cur->from[k].score_n = tmp;
	cell_cur->from[k].back_n  = tmp2;

				
	/*
	 * west
	 */
	if (!revcmpl) {
	  tmp  = cell_w->from[k].score_nw - go - ge;
	  tmp2 = FROM_x(k, FROM_WEST_NORTHWEST);

	  if (cell_w->from[k].score_w - ge > tmp) {
	    tmp  = cell_w->from[k].score_w - ge;
	    tmp2 = FROM_x(k, FROM_WEST_WEST);
	  }
	} else {
	  tmp  = cell_w->from[k].score_w - ge;
	  tmp2 = FROM_x(k, FROM_WEST_WEST);

	  if (cell_w->from[k].score_nw - go - ge > tmp) {
	    tmp  = cell_w->from[k].score_nw - go - ge;
	    tmp2 = FROM_x(k, FROM_WEST_NORTHWEST);
	  }
	}

	/*
	 * NB: It doesn't make sense to cross over on a
	 *     genomic gap, so we won't.
	 */

	if (tmp <= resetval) {
	  tmp = resetval;
	  tmp2 = 0;
	}

	cell_cur->from[k].score_w = tmp;
	cell_cur->from[k].back_w  = tmp2;


	/*
	 * max score
	 */
	if (!revcmpl) {
	  if (cell_cur->from[k].score_nw > score) {
	    score = cell_cur->from[k].score_nw;
	    max_i = i, max_j = j, max_k = k;
	  }
	  if (cell_cur->from[k].score_n > score) {
	    score = cell_cur->from[k].score_n;
	    max_i = i, max_j = j, max_k = k;
	  }
	  if (cell_cur->from[k].score_w > score) {
	    score = cell_cur->from[k].score_w;
	    max_i = i, max_j = j, max_k = k;
	  }
	} else {
	  if (cell_cur->from[k].score_w > score) {
	    score = cell_cur->from[k].score_w;
	    max_i = i, max_j = j, max_k = k;
	  }
	  if (cell_cur->from[k].score_n > score) {
	    score = cell_cur->from[k].score_n;
	    max_i = i, max_j = j, max_k = k;
	  }
	  if (cell_cur->from[k].score_nw > score) {
	    score = cell_cur->from[k].score_nw;
	    max_i = i, max_j = j, max_k = k;
	  }
	}

#ifdef DEBUG_SW
	fprintf(stderr, "i:%d j:%d k:%d score_nw:%d [%u,%s] score_n:%d [%u,%s] score_w:%d [%u,%s] \n", i+1, j+1, k,
		cell_cur->from[k].score_nw, cell_cur->from[k].back_nw & 0x3,
		(cell_cur->from[k].back_nw >> 2 == 0 ? "!" :
		 (cell_cur->from[k].back_nw >> 2 == FROM_NORTHWEST_NORTH ? "n" :
		  (cell_cur->from[k].back_nw >> 2 == FROM_NORTHWEST_NORTHWEST ? "nw" : "w"))),

		cell_cur->from[k].score_n, cell_cur->from[k].back_n & 0x3,
		(cell_cur->from[k].back_n >> 2 == 0 ? "!" :
		 (cell_cur->from[k].back_n >> 2 == FROM_NORTH_NORTH ? "n" : "nw")),

		cell_cur->from[k].score_w, cell_cur->from[k].back_w & 0x3,
		(cell_cur->from[k].back_w >> 2 == 0 ? "!" :
		 (cell_cur->from[k].back_w >> 2 == FROM_WEST_NORTHWEST ? "nw" : "w")));
#endif

      }
    }

    if (i+1 < lenb) {
      int next_x_min, next_x_max;

      get_x_range(&rectangle, lena, lenb, i+1, &next_x_min, &next_x_max);
      for (j = x_max + 1; j <= next_x_max; j++) {
	init_cell((i + 1) * (lena + 1) + (j + 1));
      }
    }
  }

  *iret = max_i;
  *jret = max_j;
  *kret = max_k;

  return (score);
}

/*
 * Fill in the backtrace in order to do a pretty printout.
 *
 * Returns the beginning matrix cell (i, j) in 'sfr->read_start' and
 * 'sfr->genome_start'.
 *
 * The return value is the first valid offset in the backtrace buffer.
 */
static int
do_backtrace(int lena, int i, int j, int k, struct sw_full_results *sfr)
{
  struct swcell *cell;
  int off, from, fromscore;

  off = (dblen + qrlen) - 1;

  cell = &swmatrix[(i + 1) * (lena + 1) + j + 1];

  from = cell->from[k].back_nw;
  fromscore = cell->from[k].score_nw;

  if (cell->from[k].score_w > fromscore) {
    from = cell->from[k].back_w;
    fromscore = cell->from[k].score_w;
  }
  if (cell->from[k].score_n > fromscore)
    from = cell->from[k].back_n;

  assert(from != 0);

  /* fill out the backtrace */
  while (i >= 0 && j >= 0) {
    assert(off >= 0);

    cell = NULL;

    /* common operations first */
    switch (from) {
    case FROM_A_NORTH_NORTH:
    case FROM_A_NORTH_NORTHWEST:
    case FROM_B_NORTH_NORTH:
    case FROM_B_NORTH_NORTHWEST:
    case FROM_C_NORTH_NORTH:
    case FROM_C_NORTH_NORTHWEST:
    case FROM_D_NORTH_NORTH:
    case FROM_D_NORTH_NORTHWEST:
      sfr->deletions++;
      sfr->read_start = i--;
      break;

    case FROM_A_WEST_WEST:
    case FROM_A_WEST_NORTHWEST:
    case FROM_B_WEST_WEST:
    case FROM_B_WEST_NORTHWEST:
    case FROM_C_WEST_WEST:
    case FROM_C_WEST_NORTHWEST:
    case FROM_D_WEST_WEST:
    case FROM_D_WEST_NORTHWEST:
      sfr->insertions++;
      sfr->genome_start = j--;
      break;

    case FROM_A_NORTHWEST_NORTH:
    case FROM_A_NORTHWEST_NORTHWEST:
    case FROM_A_NORTHWEST_WEST:
    case FROM_B_NORTHWEST_NORTH:
    case FROM_B_NORTHWEST_NORTHWEST:
    case FROM_B_NORTHWEST_WEST:
    case FROM_C_NORTHWEST_NORTH:
    case FROM_C_NORTHWEST_NORTHWEST:
    case FROM_C_NORTHWEST_WEST:
    case FROM_D_NORTHWEST_NORTH:
    case FROM_D_NORTHWEST_NORTHWEST:
    case FROM_D_NORTHWEST_WEST:
      if (db[j] == qr[k][i])
	sfr->matches++;
      else
	sfr->mismatches++;
      sfr->read_start = i--;
      sfr->genome_start = j--;
      break;

    default:
      fprintf(stderr, "INTERNAL ERROR: from = %d\n", from);
      assert(0);
    }

    /* handle match/mismatch and north */
    switch (from) {
    case FROM_A_NORTH_NORTH:
    case FROM_A_NORTH_NORTHWEST:
    case FROM_B_NORTH_NORTH:
    case FROM_B_NORTH_NORTHWEST:
    case FROM_C_NORTH_NORTH:
    case FROM_C_NORTH_NORTHWEST:
    case FROM_D_NORTH_NORTH:
    case FROM_D_NORTH_NORTHWEST:
      switch(k) {
      case 0:
	backtrace[off]= BACK_A_DELETION;
	break;
      case 1:
	backtrace[off]= BACK_B_DELETION;
	break;
      case 2:
	backtrace[off]= BACK_C_DELETION;
	break;
      case 3:
	backtrace[off]= BACK_D_DELETION;
	break;
      default:
	fprintf(stderr, "INTERNAL ERROR: k = %d\n", k);
	assert(0);
      }
      break;

    case FROM_A_WEST_WEST:
    case FROM_A_WEST_NORTHWEST:
    case FROM_B_WEST_WEST:
    case FROM_B_WEST_NORTHWEST:
    case FROM_C_WEST_WEST:
    case FROM_C_WEST_NORTHWEST:
    case FROM_D_WEST_WEST:
    case FROM_D_WEST_NORTHWEST:
      /* doesn't make sense to cross over on a genomic gap */
      backtrace[off] = BACK_INSERTION;
      break;

    case FROM_A_NORTHWEST_NORTH:
    case FROM_A_NORTHWEST_NORTHWEST:
    case FROM_A_NORTHWEST_WEST:
    case FROM_B_NORTHWEST_NORTH:
    case FROM_B_NORTHWEST_NORTHWEST:
    case FROM_B_NORTHWEST_WEST:
    case FROM_C_NORTHWEST_NORTH:
    case FROM_C_NORTHWEST_NORTHWEST:
    case FROM_C_NORTHWEST_WEST:
    case FROM_D_NORTHWEST_NORTH:
    case FROM_D_NORTHWEST_NORTHWEST:
    case FROM_D_NORTHWEST_WEST:
      switch(k) {
      case 0:
	backtrace[off] = BACK_A_MATCH_MISMATCH;
	break;
      case 1:
	backtrace[off] = BACK_B_MATCH_MISMATCH;
	break;
      case 2:
	backtrace[off] = BACK_C_MATCH_MISMATCH;
	break;
      case 3:
	backtrace[off] = BACK_D_MATCH_MISMATCH;
	break;
      default:
	fprintf(stderr, "INTERNAL ERROR: k = %d\n", k);
	assert(0);
      }
      break;

    default:
      fprintf(stderr, "INTERNAL ERROR: from = %d\n", from);
      assert(0);
    }

    /* set k */
    switch (from) {
    case FROM_A_NORTH_NORTH:
    case FROM_A_NORTH_NORTHWEST:
    case FROM_A_WEST_WEST:
    case FROM_A_WEST_NORTHWEST:
    case FROM_A_NORTHWEST_NORTH:
    case FROM_A_NORTHWEST_NORTHWEST:
    case FROM_A_NORTHWEST_WEST:
      if (k != 0) {
	backtrace[off] |= BT_CROSSOVER;
	sfr->crossovers++;
	k = 0;
      }
      break;

    case FROM_B_NORTH_NORTH:
    case FROM_B_NORTH_NORTHWEST:
    case FROM_B_WEST_WEST:
    case FROM_B_WEST_NORTHWEST:
    case FROM_B_NORTHWEST_NORTH:
    case FROM_B_NORTHWEST_NORTHWEST:
    case FROM_B_NORTHWEST_WEST:
      if (k != 1) {
	backtrace[off] |= BT_CROSSOVER;
	sfr->crossovers++;
	k = 1;
      }
      break;

    case FROM_C_NORTH_NORTH:
    case FROM_C_NORTH_NORTHWEST:
    case FROM_C_WEST_WEST:
    case FROM_C_WEST_NORTHWEST:
    case FROM_C_NORTHWEST_NORTH:
    case FROM_C_NORTHWEST_NORTHWEST:
    case FROM_C_NORTHWEST_WEST:
      if (k != 2) {
	backtrace[off] |= BT_CROSSOVER;
	sfr->crossovers++;
	k = 2;
      }
      break;

    case FROM_D_NORTH_NORTH:
    case FROM_D_NORTH_NORTHWEST:
    case FROM_D_WEST_WEST:
    case FROM_D_WEST_NORTHWEST:
    case FROM_D_NORTHWEST_NORTH:
    case FROM_D_NORTHWEST_NORTHWEST:
    case FROM_D_NORTHWEST_WEST:
      if (k != 3) {
	backtrace[off] |= BT_CROSSOVER;
	sfr->crossovers++;
	k = 3;
      }
      break;

    default:
      fprintf(stderr, "INTERNAL ERROR: from = %d\n", from);
      assert(0);
    }


    /*
     * Continue backtrace (nb: i,j and k have already been changed).
     */
    cell = &swmatrix[(i + 1) * (lena + 1) + j + 1];

    switch (from) {
    case FROM_A_NORTH_NORTH:
    case FROM_B_NORTH_NORTH:
    case FROM_C_NORTH_NORTH:
    case FROM_D_NORTH_NORTH:
      from = cell->from[k].back_n;
      break;

    case FROM_A_NORTH_NORTHWEST:
    case FROM_B_NORTH_NORTHWEST:
    case FROM_C_NORTH_NORTHWEST:
    case FROM_D_NORTH_NORTHWEST:
      from = cell->from[k].back_nw;
      break;

    case FROM_A_WEST_WEST:
    case FROM_B_WEST_WEST:
    case FROM_C_WEST_WEST:
    case FROM_D_WEST_WEST:
      from = cell->from[k].back_w;
      break;

    case FROM_A_WEST_NORTHWEST:
    case FROM_B_WEST_NORTHWEST:
    case FROM_C_WEST_NORTHWEST:
    case FROM_D_WEST_NORTHWEST:
      from = cell->from[k].back_nw;
      break;

    case FROM_A_NORTHWEST_NORTH:
    case FROM_B_NORTHWEST_NORTH:
    case FROM_C_NORTHWEST_NORTH:
    case FROM_D_NORTHWEST_NORTH:
      from = cell->from[k].back_n;
      break;

    case FROM_A_NORTHWEST_NORTHWEST:
    case FROM_B_NORTHWEST_NORTHWEST:
    case FROM_C_NORTHWEST_NORTHWEST:
    case FROM_D_NORTHWEST_NORTHWEST:
      from = cell->from[k].back_nw;
      break;

    case FROM_A_NORTHWEST_WEST:
    case FROM_B_NORTHWEST_WEST:
    case FROM_C_NORTHWEST_WEST:
    case FROM_D_NORTHWEST_WEST:
      from = cell->from[k].back_w;
      break;

    default:
      fprintf(stderr, "INTERNAL ERROR: from = %d\n", from);
      assert(0);
    }

    off--;

    if (from == 0)
      break;		  
  }

  off++;

  if (k != 0) {
    backtrace[off] |= BT_CROSSOVER;
    sfr->crossovers++;
  }

  return (off);
}

/*
 * Pretty print our alignment of 'db' and 'qr' in 'dbalign' and 'qralign'.
 *
 * i, j represent the beginning cell in the matrix.
 * k is the first valid offset in the backtrace buffer.
 */
static void
pretty_print(int i, int j, int k)
{
  char *d, *q;
  int l;

  d = dbalign;
  q = qralign;

  for (l = k; l < (dblen + qrlen); l++) {
    switch (BT_TYPE(backtrace[l])) {
    case BACK_A_DELETION:
      *d++ = '-';
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[0][i++], false));
      else	
	*q++ = base_translate(qr[0][i++], false);
      break;

    case BACK_B_DELETION:
      *d++ = '-';
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[1][i++], false));
      else	
	*q++ = base_translate(qr[1][i++], false);
      break;

    case BACK_C_DELETION:
      *d++ = '-';
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[2][i++], false));
      else	
	*q++ = base_translate(qr[2][i++], false);
      break;

    case BACK_D_DELETION:
      *d++ = '-';
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[3][i++], false));
      else	
	*q++ = base_translate(qr[3][i++], false);
      break;

    case BACK_INSERTION:
      *d++ = base_translate(db[j++], false);
      *q++ = '-';
      break;

    case BACK_A_MATCH_MISMATCH:
      *d++ = base_translate(db[j++], false);
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[0][i++], false));
      else	
	*q++ = base_translate(qr[0][i++], false);
      break;

    case BACK_B_MATCH_MISMATCH:
      *d++ = base_translate(db[j++], false);
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[1][i++], false));
      else	
	*q++ = base_translate(qr[1][i++], false);
      break;

    case BACK_C_MATCH_MISMATCH:
      *d++ = base_translate(db[j++], false);
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[2][i++], false));
      else	
	*q++ = base_translate(qr[2][i++], false);
      break;

    case BACK_D_MATCH_MISMATCH:
      *d++ = base_translate(db[j++], false);
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[3][i++], false));
      else	
	*q++ = base_translate(qr[3][i++], false);
      break;
	
    default:
      fprintf(stderr, "INTERNAL ERROR: backtrace[l] = 0x%x\n", backtrace[l]);
      assert(0);
    }
  }

  *d = *q = '\0';
}

int
sw_full_cs_setup(int _dblen, int _qrlen, int _gap_open, int _gap_ext,
		 int _match, int _mismatch, int _xover_penalty, bool reset_stats,
		 int _anchor_width)
{
  int i;

  dblen = _dblen;
  db = (int8_t *)malloc(dblen * sizeof(db[0]));
  if (db == NULL)
    return (1);

  qrlen = _qrlen;
  for (i = 0; i < 4; i++) {
    qr[i] = (int8_t *)malloc(qrlen * sizeof(qr[0]));
    if (qr[i] == NULL)
      return (1);
  }

  swmatrix = (struct swcell *)malloc((dblen + 1) * (qrlen + 1) *
				     sizeof(swmatrix[0]));
  if (swmatrix == NULL)
    return (1);

  backtrace = (uint8_t *)malloc((dblen + qrlen) * sizeof(backtrace[0]));
  if (backtrace == NULL)
    return (1);

  dbalign = (char *)malloc((dblen + qrlen + 1) * sizeof(dbalign[0]));
  if (dbalign == NULL)
    return (1);

  qralign = (char *)malloc((dblen + qrlen + 1) * sizeof(dbalign[0]));
  if (qralign == NULL)
    return (1);

  gap_open = -(_gap_open);
  gap_ext = -(_gap_ext);
  match = _match;
  mismatch = _mismatch;
  xover_penalty = _xover_penalty;

  if (reset_stats)
    swticks = swcells = swinvocs = 0;

  anchor_width = _anchor_width;

  initialised = 1;

  return (0);
}

void
sw_full_cs_stats(uint64_t *invoc, uint64_t *cells, uint64_t *ticks,
		 double *cellspersec)
{
	
  if (invoc != NULL)
    *invoc = swinvocs;
  if (cells != NULL)
    *cells = swcells;
  if (ticks != NULL)
    *ticks = swticks;
  if (cellspersec != NULL) {
    *cellspersec = (double)swcells / ((double)swticks / cpuhz());
    if (isnan(*cellspersec))
      *cellspersec = 0;
  }
}

void
sw_full_cs(uint32_t *genome_ls, int goff, int glen, uint32_t *read, int rlen,
	   int initbp, int threshscore, struct sw_full_results *sfr, bool revcmpl, bool is_rna,
	   struct anchor * anchors, uint anchors_cnt)
{
  struct sw_full_results scratch;
  uint64_t before;
  int i, j, k;

  before = rdtsc();

  if (!initialised)
    abort();

  swinvocs++;

  assert(glen > 0 && glen <= dblen);
  assert(rlen > 0 && rlen <= qrlen);

  if (sfr == NULL)
    sfr = &scratch;

  memset(sfr, 0, sizeof(*sfr));
  memset(backtrace, 0, (dblen + qrlen) * sizeof(backtrace[0]));

  dbalign[0] = qralign[0] = '\0';

  for (i = 0; i < glen; i++)
    db[i] = (int8_t)EXTRACT(genome_ls, goff + i);

  /*
   * Generate each possible letter space sequence from the colour space
   * read. qr[0] corresponds to initbp, which is given initial preference.
   */
  assert(initbp >= 0 && initbp <= 3);
  for (i = 0; i < 4; i++) {
    int letter = (i + initbp) % 4;

    for (j = 0; j < rlen; j++) {
      int base = EXTRACT(read, j);

      if (base == BASE_N || base == BASE_X) {
	qr[i][j] = BASE_N;
	letter = (i + initbp) % 4;
      } else {
	qr[i][j] = (int8_t)cstols(letter, base, is_rna);
	letter = qr[i][j];
      }
    }
  }

#ifdef DEBUG_SW
  fprintf(stderr, "db: ");
  for (j = 0; j < glen; j++)
    fprintf(stderr, "%c", base_translate(db[j], false));
  fprintf(stderr, "\n");
  for (i = 0; i < 4; i++) {
    fprintf(stderr, "qr[%u]: ", i);
    for (j = 0; j < rlen; j++)
      fprintf(stderr, "%c", base_translate(qr[i][j], false));
    fprintf(stderr, "\n");
  }
#endif

  sfr->score = full_sw(glen, rlen, threshscore, &i, &j, &k, revcmpl, anchors, anchors_cnt);
  k = do_backtrace(glen, i, j, k, sfr);
  pretty_print(sfr->read_start, sfr->genome_start, k);
  sfr->gmapped = j - sfr->genome_start + 1;
  sfr->genome_start += goff;
  sfr->rmapped = i - sfr->read_start + 1;
  sfr->dbalign = xstrdup(dbalign);
  sfr->qralign = xstrdup(qralign);

  //swcells += (glen * rlen);
  swticks += (rdtsc() - before);
}
