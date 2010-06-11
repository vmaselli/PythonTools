/*	$Id: output.c 383 2009-09-30 22:47:56Z matei $	*/

#include <assert.h>
#include <ctype.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "../common/fasta.h"
#include "../common/sw-full-common.h"
#include "../common/input.h"
#include "../common/output.h"
#include "../common/util.h"

char *
readtostr(const uint32_t *read, u_int len, bool use_colours, int initbp)
{
	char *buf;
	u_int i, j;

	buf = (char *)xmalloc(len + 2);

	i = 0;
	if (use_colours)
		buf[i++] = base_translate(initbp, false);

	for (j = 0; j < len; j++) {
		if (use_colours)
			buf[i + j] = base_translate(EXTRACT(read, j), true);
		else
			buf[i + j] = base_translate(EXTRACT(read, j), false); 
	}

	buf[i + j] = '\0';

	return (buf);
}

/*
 * Create a Phil Edit String:
 *	The proposal for the edit strings look like this:
 *	<number> = size of a matching substring
 *	<letter> = mismatch, value is the tag letter
 *	(<letters>) = gap in the reference, value shows the letters in the tag
 * 	- = one-base gap in the tag
 * 	x = crossover (inserted between the appropriate two bases)
 *
 *	A perfect match for 25-bp tags is: 25
 *	A SNP at the 16th base of the tag is: 15A9
 *	A four-base insertion in the reference: 3(TGCT)20
 *	A four-base deletion in the reference: 5----20
 *	Two sequencing errors: 4x15x6	(25 total matches)
 *	etc.
 */
static char *
editstr(const char *dbalign, const char *qralign)
{
	strbuf_t sb;
	char *str;
	int i, len, consec;
	bool refgap = false;

	len = strlen(dbalign);
	assert(len == (int)strlen(qralign));

	sb = strbuf_create();

	consec = 0;
	for (i = 0; i <= len; i++) {
		if (i != len && dbalign[i] == qralign[i] && dbalign[i] != '-') {
			consec++;
			continue;
		}

		if (refgap && (consec != 0 || dbalign[i] != '-')) {
			strbuf_append(sb, ")");
			refgap = false;
		}

		if (consec != 0) {
			strbuf_append(sb, "%d", consec);
			consec = 0;
		}

		if (i == len)
			break;

		if (dbalign[i] == '-') {
			if (islower((int)qralign[i]))
				strbuf_append(sb, "x");
			if (!refgap)
				strbuf_append(sb, "(");
			strbuf_append(sb, "%c", toupper((int)qralign[i]));
			refgap = true;
			continue;
		}

		if (qralign[i] == '-') {
			strbuf_append(sb, "-");
		} else {
			assert(i > 0 || dbalign[i] == toupper((int)qralign[i]));
			if (dbalign[i] == toupper((int)qralign[i])) {
				strbuf_append(sb, "x");
				consec++;
			} else if (islower((int)qralign[i])) {
				strbuf_append(sb, "x");
				strbuf_append(sb, "%c", toupper((int)qralign[i]));
			} else {
				strbuf_append(sb, "%c", qralign[i]);
			}
		}
	}

	str = strbuf_string(sb, NULL);
	strbuf_destroy(sb);
	return (str);
}

/* This is so ugly it hurts. */
char *
output_pretty(const char *readname, const char *contigname,
    const struct sw_full_results *sfr, uint32_t *genome, uint32_t genome_len,
    bool use_colours, uint32_t *read, u_int readlen, int initbp, bool revcmpl)
{
	char *str, *gpre, *gpost, *lspre, *lspost, *mpre;
	char const * nospace = "";
	const char *dbalign, *qralign;
	strbuf_t sb;
	u_int j, len;
	uint32_t genome_start, genome_end;		/* start w/in *genome */
	uint32_t idx_genome_start, idx_genome_end;	/* relative idx user sees */
	uint32_t read_start, read_end;

	sb = strbuf_create();

	dbalign = sfr->dbalign;
	qralign = sfr->qralign;

	/* shut up, icc */
	(void)readname;
	(void)genome;
	(void)contigname;

	genome_start = sfr->genome_start;
	genome_end   = sfr->genome_start + sfr->gmapped - 1;

	assert(genome_len > genome_start);
	assert(genome_len > genome_end);
	assert(genome_end > genome_start);

	if (revcmpl) {
		idx_genome_start = genome_len - genome_end - 1;
		idx_genome_end = genome_len - genome_start - 1;
	} else {
		idx_genome_start = genome_start;
		idx_genome_end = genome_end;
	}

	read_start = sfr->read_start;
	read_end   = sfr->read_start + sfr->rmapped - 1;

	len = strlen(dbalign);
	assert(len == strlen(qralign));

	gpre = gpost = lspre = lspost = mpre = (char *)nospace;
	if (read_start > 0) {
		gpre  = (char *)xmalloc(read_start + 1);
		lspre = (char *)xmalloc(read_start + 1);
		mpre  = (char *)xmalloc(read_start + 1);
		for (j = 0; j < read_start; j++) {
			if (genome_start + j > read_start)
				gpre[j] = base_translate(EXTRACT(genome,
				    genome_start - read_start + j), false);
			else
				gpre[j] = '-';
			lspre[j] = '-';
			mpre[j] = ' ';
		}
		gpre[j] = lspre[j] = mpre[j] = '\0';
	}
	if (read_end < (readlen - 1)) {
		gpost  = (char *)xmalloc(readlen - read_end);
		lspost = (char *)xmalloc(readlen - read_end);
		for (j = 0; j < (readlen - read_end - 1); j++) {
			if (genome_end + 1 + j < genome_len)
				gpost[j] = base_translate(EXTRACT(genome,
				    genome_end + 1 + j), false);
			else
				gpost[j] = '-';
			lspost[j] = '-';
		}
		gpost[j] = lspost[j] = '\0';
	}

	/* NB: internally 0 is first position, output uses 1. adjust. */

	strbuf_append(sb, "G: %10" PRId64 "    %s%s%s    %-10" PRId64 "\n",
	    (revcmpl) ? (int64_t)idx_genome_end + 1 : (int64_t)idx_genome_start + 1,
	    gpre, dbalign, gpost,
	    (revcmpl) ? (int64_t)idx_genome_start + 1 : (int64_t)idx_genome_end + 1);

	strbuf_append(sb, "%16s %s", "", mpre);
	for (j = 0; j < len; j++) {
		if (dbalign[j] == qralign[j] && dbalign[j] != '-') {
			strbuf_append(sb, "|");
		} else {
			assert(j > 0 || dbalign[j] == toupper((int)qralign[j]));
			if (dbalign[j] == toupper((int)qralign[j]))
				strbuf_append(sb, "X");
			else if (islower((int)qralign[j]))
				strbuf_append(sb, "x");
			else
				strbuf_append(sb, " ");
		}
	}
	strbuf_append(sb, "\n");

	if (use_colours) {
		strbuf_append(sb, "T: %10s    %s%s%s\n", "", lspre, qralign, lspost);
	} else {
		strbuf_append(sb, "R: %10u    %s%s%s    %-10u\n", read_start + 1,
		    lspre, qralign, lspost, read_end + 1);
	}

	if (use_colours) {
		char *rstr, *rstrb;

		strbuf_append(sb, "R: %10u   ", read_start + 1);
		rstr = rstrb = readtostr(read, readlen, use_colours, initbp);
		assert(strlen(rstr) > 1);
		strbuf_append(sb, "%c", *rstr++);
		for (j = 0; j < read_start; j++) {
			assert(*rstr != '\0');
			strbuf_append(sb, "%c", *rstr++);
		}
		for (j = 0; *rstr != '\0';) {
			if (qralign[j] == '-')
				strbuf_append(sb, "-");
			else
				strbuf_append(sb, "%c", *rstr++);
			if (qralign[j] != '\0')
				j++;
		}
		strbuf_append(sb, "    %-10u\n", read_end + 1);
		free(rstrb);
	}

	if (gpre != nospace)
		free(gpre);
	if (gpost != nospace)
		free(gpost);
	if (lspre != nospace)
		free(lspre);
	if (lspost != nospace)
		free(lspost);
	if (mpre != nospace)
		free(mpre);

	str = strbuf_string(sb, NULL);
	strbuf_destroy(sb);
	return (str);
}

char *
output_format_line(bool inc_read)
{
	strbuf_t sb = strbuf_create();
	char *str;

	strbuf_append(sb, "#FORMAT: readname contigname strand contigstart contigend readstart readend "
	    "readlength score editstring");
	if (inc_read)
		strbuf_append(sb, " readsequence");

	str = strbuf_string(sb, NULL);
	strbuf_destroy(sb);
	return (str);
}

char *
output_normal(const char *readname, const char *contigname,
    const struct sw_full_results *sfr, uint32_t genome_len, bool use_colours,
    uint32_t *read, u_int readlen, int initbp, bool revcmpl, bool inc_read)
{
	struct sw_full_results sfr_tmp;
	uint32_t genome_start, genome_end;
	uint32_t idx_genome_start, idx_genome_end;
	const char *dbalign, *qralign;
	char *edit, *readseq, *str;
	strbuf_t sb;
	bool ret;

	sb = strbuf_create();

	dbalign = sfr->dbalign;
	qralign = sfr->qralign;

	assert(readname != NULL);
	assert(contigname != NULL);
	assert(read != NULL);

	genome_start = sfr->genome_start;
	genome_end = sfr->genome_start + sfr->gmapped - 1;

	assert(genome_len > genome_start);
	assert(genome_len > genome_end);
	assert(genome_end > genome_start);

	if (inc_read)
		readseq = readtostr(read, readlen, use_colours, initbp);
	else
		readseq = xstrdup("");

	if (revcmpl) {
		idx_genome_start = genome_len - genome_end - 1;
		idx_genome_end = genome_len - genome_start - 1;
	} else {
		idx_genome_start = genome_start;
		idx_genome_end = genome_end;
	}

	strbuf_append(sb, ">%s\t%s\t%c", readname, contigname,
	    (revcmpl) ? '-' : '+');

	edit = editstr(dbalign, qralign);
	ret = editstr_to_sfr(edit, &sfr_tmp);
	assert(ret &&
	    sfr->matches == sfr_tmp.matches &&
	    sfr->mismatches == sfr_tmp.mismatches &&
	    sfr->crossovers == sfr_tmp.crossovers &&
	    sfr->insertions == sfr_tmp.insertions &&
	    sfr->deletions == sfr_tmp.deletions);

	strbuf_append(sb, "\t%u\t%u\t%d\t%d\t%d\t%d\t%s\t%s", idx_genome_start + 1,
	    idx_genome_end + 1, sfr->read_start + 1,
	    sfr->read_start + sfr->rmapped - 1 + 1, readlen, sfr->score, edit, readseq);

	free(readseq);
	free(edit);

	str = strbuf_string(sb, NULL);
	strbuf_destroy(sb);
	return (str);
}
