/*	$Id: fasta.h 348 2009-06-16 23:26:27Z rumble $	*/

#ifndef _FASTA_H_
#define _FASTA_H_

#define LETTER_SPACE	1
#define COLOUR_SPACE	2

/*
 * We're presently using 4 bits, so we max out at 16 different bases. This just
 * works out if we treat N and X as the same.
 *
 * from: http://genome.ucsc.edu/FAQ/FAQdownloads#download5
 */

#define BASE_A		0		/* Adenine */
#define BASE_C		1		/* Cytosine */
#define BASE_G		2		/* Guanine */
#define BASE_T		3		/* Thymine */
#define BASE_U		4		/* Uracil */
#define BASE_M		5		/* A or C */
#define BASE_R		6		/* A or G (Purine) */
#define BASE_W		7		/* A or T */
#define BASE_S		8		/* C or G */
#define BASE_Y		9		/* C or T (Pyrimidine) */
#define BASE_K		10		/* G or T */
#define BASE_V		11		/* A or C or G (not T) */
#define BASE_H		12		/* A or C or T (not G) */
#define BASE_D		13		/* A or G or T (not C) */
#define BASE_B		14		/* C or G or T (not A) */
#define BASE_X		15		/* G or A or T or C (any base) */
#define BASE_N		15		/* G or A or T or C (any base) */

#define BASE_LS_MIN	BASE_A
#define BASE_LS_MAX	BASE_B

#define BASE_0		0
#define BASE_1		1
#define BASE_2		2
#define BASE_3		3

#define BASE_CS_MIN	BASE_0
#define BASE_CS_MAX	BASE_3

typedef struct _fasta_t {
	gzFile fp;
	char  *file;
	int    space;
	char   buffer[8*1024*1024];
	char   translate[256];
	bool   leftover;
} * fasta_t;

typedef struct _fasta_stats_t {
	uint64_t	total_ticks;
} * fasta_stats_t;

fasta_t	  fasta_open(const char *, int);
void	  fasta_close(fasta_t);
bool	  fasta_get_next(fasta_t, char **, char **, bool *);
int	  fasta_get_initial_base(fasta_t, char *);
uint32_t *fasta_bitfield_to_colourspace(fasta_t, uint32_t *, uint32_t, bool);
uint32_t *fasta_sequence_to_bitfield(fasta_t, char *);
fasta_stats_t fasta_stats(void);
char      base_translate(int, bool);

#endif
