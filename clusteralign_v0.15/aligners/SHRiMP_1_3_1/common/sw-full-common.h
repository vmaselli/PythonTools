/*	$Id: sw-full-common.h 245 2008-06-06 18:24:28Z rumble $	*/

struct sw_full_results {
	/* Common fields */
	int read_start;				/* read index of map */
	int rmapped;				/* read mapped length */
	int genome_start;			/* genome index of map (abs) */
	int gmapped;				/* genome mapped len */
	int matches;				/* # of matches */
	int mismatches;				/* # of substitutions */
	int insertions;				/* # of insertions */
	int deletions;				/* # of deletions */
	int score;				/* final SW score */

	char *dbalign;				/* genome align string */
	char *qralign;				/* read align string */

	/* Colour space fields */
	int crossovers;				/* # of mat. xovers */
};

static inline bool
sw_full_results_equal(struct sw_full_results *sfr1, struct sw_full_results *sfr2)
{
	char *dbalign1, *dbalign2;
	char *qralign1, *qralign2;
	bool equal;

	dbalign1 = sfr1->dbalign;
	dbalign2 = sfr2->dbalign;
	qralign1 = sfr1->qralign;
	qralign2 = sfr2->qralign;

	sfr1->dbalign = sfr2->dbalign = NULL;
	sfr1->qralign = sfr2->qralign = NULL;

	equal = memcmp(sfr1, sfr2, sizeof(*sfr1)) == 0;

	sfr1->dbalign = dbalign1;
	sfr2->dbalign = dbalign2;
	sfr1->qralign = qralign1;
	sfr2->qralign = qralign2;

	return (equal);
}
