/*	$Id: fasta.c 348 2009-06-16 23:26:27Z rumble $	*/

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
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
#include "../common/util.h"

static uint64_t total_ticks;

fasta_t
fasta_open(const char *file, int space)
{
	fasta_t fasta = NULL;
	struct stat sb;
	gzFile fp;
	uint64_t before = rdtsc();

	assert(space == COLOUR_SPACE || space == LETTER_SPACE);

	if (stat(file, &sb))
		goto out;

	if (!S_ISREG(sb.st_mode))
		goto out;

	fp = gzopen(file, "r");
	if (fp == NULL)
		goto out;

	fasta = (fasta_t)xmalloc(sizeof(*fasta));
	memset(fasta, 0, sizeof(*fasta));

	fasta->fp = fp;
	fasta->file = xstrdup(file);
	fasta->space = space;
	memset(fasta->translate, -1, sizeof(fasta->translate));

	if (space == COLOUR_SPACE) {
		fasta->translate[(int)'0'] = BASE_0;
		fasta->translate[(int)'1'] = BASE_1;
		fasta->translate[(int)'2'] = BASE_2;
		fasta->translate[(int)'3'] = BASE_3;
		fasta->translate[(int)'4'] = BASE_N;
		fasta->translate[(int)'N'] = BASE_N;
		fasta->translate[(int)'n'] = BASE_N;
		fasta->translate[(int)'.'] = BASE_N;
		fasta->translate[(int)'X'] = BASE_X;
		fasta->translate[(int)'x'] = BASE_X;
	} else {	
		fasta->translate[(int)'A'] = BASE_A;
		fasta->translate[(int)'a'] = BASE_A;
		fasta->translate[(int)'C'] = BASE_C;
		fasta->translate[(int)'c'] = BASE_C;
		fasta->translate[(int)'G'] = BASE_G;
		fasta->translate[(int)'g'] = BASE_G;
		fasta->translate[(int)'T'] = BASE_T;
		fasta->translate[(int)'t'] = BASE_T;
		fasta->translate[(int)'U'] = BASE_U;
		fasta->translate[(int)'u'] = BASE_U;
		fasta->translate[(int)'M'] = BASE_M;
		fasta->translate[(int)'m'] = BASE_M;
		fasta->translate[(int)'R'] = BASE_R;
		fasta->translate[(int)'r'] = BASE_R;
		fasta->translate[(int)'W'] = BASE_W;
		fasta->translate[(int)'w'] = BASE_W;
		fasta->translate[(int)'S'] = BASE_S;
		fasta->translate[(int)'s'] = BASE_S;
		fasta->translate[(int)'Y'] = BASE_Y;
		fasta->translate[(int)'y'] = BASE_Y;
		fasta->translate[(int)'K'] = BASE_K;
		fasta->translate[(int)'k'] = BASE_K;
		fasta->translate[(int)'V'] = BASE_V;
		fasta->translate[(int)'v'] = BASE_V;
		fasta->translate[(int)'H'] = BASE_H;
		fasta->translate[(int)'h'] = BASE_H;
		fasta->translate[(int)'D'] = BASE_D;
		fasta->translate[(int)'d'] = BASE_D;
		fasta->translate[(int)'B'] = BASE_B;
		fasta->translate[(int)'b'] = BASE_B;
		fasta->translate[(int)'N'] = BASE_N;
		fasta->translate[(int)'n'] = BASE_N;
		fasta->translate[(int)'.'] = BASE_N;
		fasta->translate[(int)'X'] = BASE_X;
		fasta->translate[(int)'x'] = BASE_X;
	}

 out:
	total_ticks += (rdtsc() - before);
	return (fasta);
}

void
fasta_close(fasta_t fasta)
{
	uint64_t before = rdtsc();

	gzclose(fasta->fp);
	free(fasta->file);
	free(fasta);

	total_ticks += (rdtsc() - before);
}

fasta_stats_t
fasta_stats()
{
	fasta_stats_t fs;

	fs = (fasta_stats_t)xmalloc(sizeof(*fs));
	fs->total_ticks = total_ticks;

	return (fs);
}

static char *
extract_name(char *buffer)
{
	char *extracted;
	char *ret;
	int len;

	assert(buffer[0] == '>');

	/* Funny business for valgrind. See bottom of fasta_get_next. */
	extracted = strtrim(&buffer[1]);
	len = strlen(extracted);
	ret = (char *)xmalloc(len + 17);
	memcpy(ret, extracted, len);
	memset(ret + len, 0, 17);

	return (ret);
}

bool
fasta_get_next(fasta_t fasta, char **name, char **sequence, bool *is_rna)
{
	static int categories[256] = { 42 };
	enum { CAT_SPACE, CAT_NEWLINE, CAT_NUL, CAT_ELSE };

	int i;
	bool gotname = false;
	uint32_t sequence_length = 0;
	uint32_t max_sequence_length = 0;
	uint64_t before = rdtsc();
	char *readinseq;

	/* isspace(3) is really slow, so build a lookup instead */
	if (categories[0] == 42) {
		for (i = 0; i < 256; i++) {
			if (i == '\n')
				categories[i] = CAT_NEWLINE;
			else if (i == '\0')
				categories[i] = CAT_NUL;
			else if (isspace(i))
				categories[i] = CAT_SPACE;
			else if (i == '-')	/* simply ignore dashes in haplome alignments */
				categories[i] = CAT_SPACE;
			else
				categories[i] = CAT_ELSE;
		}
	}

	*name = *sequence = readinseq = NULL;

	if (fasta->leftover) {
		*name = extract_name(fasta->buffer);
		fasta->leftover = false;
		gotname = true;
	}

	while (fast_gzgets(fasta->fp, fasta->buffer, sizeof(fasta->buffer)) != NULL) {
		if (fasta->buffer[0] == '#')
			continue;

		if (fasta->buffer[0] == '>') {
			if (gotname) {
				fasta->leftover = true;
				break;
			}

			*name = extract_name(fasta->buffer);
			gotname = true;
			continue;
		}

		if (max_sequence_length == 0 ||
		    max_sequence_length - sequence_length <= sizeof(fasta->buffer)) {
			if (max_sequence_length == 0)
				max_sequence_length = 16*1024*1024;
			else if (max_sequence_length >= 128*1024*1024)
				max_sequence_length += 128*1024*1024;
			else
				max_sequence_length *= 2;
			readinseq = (char *)xrealloc(readinseq, max_sequence_length);
		}

		for (i = 0; fasta->buffer[i] != '\0'; i++) {
			int act = categories[(int)fasta->buffer[i]];

			assert(i < (int)sizeof(fasta->buffer));

			switch (act) {
			case CAT_NUL:
			case CAT_NEWLINE:
				fasta->buffer[i] = '\0';
				goto out;
			case CAT_SPACE:
				continue;
			case CAT_ELSE:
			default:
				/* fall down below */
				break;
			}

			readinseq[sequence_length++] = fasta->buffer[i];

			continue;
 out:
			break;
		}

		assert(sequence_length <= max_sequence_length);
	}

	if (!gotname) {
		if (readinseq != NULL)
			free(readinseq);
		total_ticks += (rdtsc() - before);
		return (false);
	}

	/* if we read in a tag name but there's no sequence, don't just return false */ 
	if (gotname && readinseq == NULL) {
		readinseq = xstrdup("");
		sequence_length = 1;
	}

	if (readinseq == NULL)
		return (false);

	/*
	 * Ensure nul-termination and allocate extra space to appease valgrind.
	 * I think strlen may do > char-sized loads and valgrind thinks it's
	 * conditionally jumping on uninitialised memory. That explanation
	 * doesn't seem right to me, but doing this makes it happy again.
	 */
	*sequence = (char *)xmalloc(sequence_length + 17);
	memcpy(*sequence, readinseq, sequence_length);
	memset(*sequence + sequence_length, 0, 17);

	/* check if the sequence is rna (contains uracil and not thymine) */
	if (is_rna != NULL) {
		bool got_uracil = false;
		bool got_thymine = false;
		unsigned int j;

		for (j = 0; j < sequence_length; j++) {
			unsigned char chr = readinseq[j];
			got_thymine |= (chr == 'T' || chr == 't');
			got_uracil  |= (chr == 'U' || chr == 'u');
		}

		if (got_uracil && got_thymine)
			fprintf(stderr, "WARNING: sequence has both uracil and "
			    "thymine!?!\n");

		*is_rna = (got_uracil && !got_thymine);
	}

	free(readinseq);
	assert(*name != NULL);
	assert(*sequence != NULL);
	total_ticks += (rdtsc() - before);
	return (true);
}

int
fasta_get_initial_base(fasta_t fasta, char *sequence)
{

	assert(fasta->space == COLOUR_SPACE);

	switch (*sequence) {
	case 'A':
	case 'a':
		return (BASE_A);
	case 'C':
	case 'c':
		return (BASE_C);
	case 'G':
	case 'g':
		return (BASE_G);
	case 'T':
	case 't':
		return (BASE_T);
	}

	return (-1);
}

uint32_t *
fasta_bitfield_to_colourspace(fasta_t fasta, uint32_t *source, uint32_t length, bool is_rna)
{
	int a, lastbp = BASE_T;
	uint32_t *dst;
	uint32_t i;
	uint64_t before = rdtsc();

	assert(fasta->space == LETTER_SPACE);

	dst = (uint32_t *)xmalloc(BPTO32BW(length) * sizeof(uint32_t));
	memset(dst, 0, BPTO32BW(length) * sizeof(uint32_t));

	for (i = 0; i < length; i++) {
		a = EXTRACT(source, i);
		bitfield_insert(dst, i, lstocs(lastbp, a, is_rna));
		lastbp = a;
	}

	total_ticks += (rdtsc() - before);
	return (dst);
}

uint32_t *
fasta_sequence_to_bitfield(fasta_t fasta, char *sequence)
{
	uint32_t i, length, idx;
	uint32_t *bitfield;
	uint64_t before = rdtsc();
	int a;
	char c;

	length = strlen(sequence);
	bitfield = (uint32_t *)xmalloc(BPTO32BW(length) * sizeof(uint32_t));
	memset(bitfield, 0, BPTO32BW(length) * sizeof(uint32_t));

	i = 0;
	c = (char)sequence[0];
	if (fasta->space == COLOUR_SPACE) {
		if (c != 'A' && c != 'a' &&
		    c != 'C' && c != 'c' && 
		    c != 'G' && c != 'g' &&
		    c != 'T' && c != 't') {
			free(bitfield);
			total_ticks += (rdtsc() - before);
			return (NULL);
		}
		
		i = 1;
	}

	for (idx = 0; i < length; i++) {
		a = fasta->translate[(int)sequence[i]];
		if (a == -1) {
			fprintf(stderr, "error: invalid character ");
			if (isprint(a))
				fprintf(stderr, "(%c) ", a);
			else if (a != -1)
				fprintf(stderr, "(0x%x) ", a);
			fprintf(stderr, "in input file [%s]\n", fasta->file);
			fprintf(stderr, "       (Did you mix up letter "
			    "space and colour space programs?)\n");
			exit(1);
		}

		if (fasta->space == COLOUR_SPACE) {
			assert((a >= BASE_CS_MIN && a <= BASE_CS_MAX) ||
			    (a == BASE_N || a == BASE_X));
		} else {
			assert((a >= BASE_LS_MIN && a <= BASE_LS_MAX) ||
			    (a == BASE_N || a == BASE_X));
		}

		bitfield_append(bitfield, idx++, a);
	}

	if (fasta->space == COLOUR_SPACE)
		assert(idx == length - 1);
	else
		assert(idx == length);

	total_ticks += (rdtsc() - before);
	return (bitfield);
}

/*
 * Give BASE_x, return the appropriate character.
 *
 * NB: Since we're limited to 4-bits, BASE_X returns 'N'.
 */
char
base_translate(int base, bool use_colours)
{
	/*
	 * NB: colour-space only valid for 0-3 and BASE_N/BASE_X
	 *     BASE_N is reported as a skipped cycle: '.' in CS.
	 */
	char cstrans[] = { '0', '1', '2', '3', '!', '@', '#', '$',
			   '%', '^', '&', '*', '?', '~', ';', '.' };
	char lstrans[] = { 'A', 'C', 'G', 'T', 'U', 'M', 'R', 'W',
			   'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N' };

	if (use_colours) {
		assert((base >= BASE_CS_MIN && base <= BASE_CS_MAX) ||
		    (base == BASE_N || base == BASE_X));
		return (cstrans[base]);
	} else {
		assert((base >= BASE_LS_MIN && base <= BASE_LS_MAX) ||
		    (base == BASE_N || base == BASE_X));
		return (lstrans[base]);
	}
}
