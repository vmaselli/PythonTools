/*	$Id: sw-full-cs.h 383 2009-09-30 22:47:56Z matei $	*/

#ifndef _SW_FULL_CS_H
#define _SW_FULL_CS_H

#include "anchors.h"

int	sw_full_cs_setup(int, int, int, int, int, int, int, bool, int);
void	sw_full_cs_stats(uint64_t *, uint64_t *, uint64_t *, double *);
void	sw_full_cs(uint32_t *, int, int, uint32_t *, int, int, int,
		   struct sw_full_results *, bool, bool, struct anchor *, uint);

#endif
