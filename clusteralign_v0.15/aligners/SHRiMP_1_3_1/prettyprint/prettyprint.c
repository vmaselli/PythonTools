/*	$Id: prettyprint.c 383 2009-09-30 22:47:56Z matei $	*/

#include <assert.h>
#include <dirent.h>
#include <errno.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
 
#include <sys/types.h>
#include <sys/stat.h>

#include "../common/dag_glue.h"
#include "../common/fasta.h"
#include "../common/dynhash.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/sw-full-ls.h"
#include "../common/input.h"
#include "../common/output.h"
#include "../common/util.h"
#include "../common/version.h"

#include "../rmapper/rmapper.h"		/* for External parameters below */

/* External parameters */
static int match_value    = DEF_MATCH_VALUE;
static int mismatch_value = DEF_MISMATCH_VALUE;
static int a_gap_open	  = DEF_A_GAP_OPEN;
static int b_gap_open	  = DEF_B_GAP_OPEN;
static int a_gap_extend   = DEF_A_GAP_EXTEND;
static int b_gap_extend   = DEF_B_GAP_EXTEND;
static int xover_penalty  = DEF_XOVER_PENALTY;

static dynhash_t read_list;		/* cache of reads we need */
static dynhash_t contig_list;		/* cache of reads per contig list */

static uint32_t nread_files;
static uint64_t nread_bases;
static uint32_t ncontig_files;
static uint64_t ncontig_bases;
static uint32_t longest_read_len;	/* longest read we have */

/* probcalc/rmapper output cache */
struct fpo {
	struct input input;		/* input fields */
	struct fpo  *next_ordered;	/* global list in input order */
	struct fpo  *next_contig;	/* per-contig/revcmpl list */
	char        *output_normal;	/* alignment output */
	char        *output_pretty;	/* alignment output */
};

struct contig_ll {
	struct fpo  *head;		/* per-contig list */
	struct fpo  *head_revcmpl;
};

/* for read and contig lists */
struct sequence {
	char	 *name;
	uint32_t *sequence;
	uint32_t  sequence_len;
	int	  initbp;		/* initial base pair (colourspace) */
	bool	  revcmpl;
	bool	  is_rna;
};

static struct fpo *alignments_ordered;	/* alignments in order of input */

static uint64_t    nalignments;
static uint64_t    nalignments_revcmpl;

static bool Rflag = false;		/* don't output read sequence */
static bool Tflag = false;		/* reverse tie-breaks on neg strand */

static bool seen_probs = false;		/* set if ever seen normodds, etc */

#define MAX_READ_LEN	5000		/* ridiculously high upper bound */

/*
 * Run through our linked list, doing S-W on the appropriate contig
 * and read and pretty print to stdout.
 */
static void
compute_alignment(struct fpo *fpo, struct sequence *contig)
{
	struct sw_full_results sfr;
	struct sequence *read;
	uint32_t genome_start, genome_len;
	bool revcmpl;

	if (!dynhash_find(read_list, fpo->input.read, NULL,
	    (void **)(void *)&read)) {
		fprintf(stderr, "error: read [%s] is missing\n",
		    fpo->input.read);
		exit(1);
	}

	revcmpl = INPUT_IS_REVCMPL(&fpo->input);
	assert(contig->revcmpl == revcmpl);

	genome_start = fpo->input.genome_start;
	genome_len = fpo->input.genome_end - genome_start + 1;

	/* offset is given relative to the positive strand */
	if (revcmpl) {
		genome_start =
		    contig->sequence_len - fpo->input.genome_end - 1;
	}

	if (shrimp_mode == MODE_COLOUR_SPACE) {
		sw_full_cs(contig->sequence, genome_start, genome_len,
		    read->sequence, read->sequence_len, read->initbp,
		    fpo->input.score, &sfr, revcmpl && Tflag,
		    contig->is_rna, NULL, 0);
	} else {
		sw_full_ls(contig->sequence, genome_start, genome_len,
		    read->sequence, read->sequence_len,
		    fpo->input.score, fpo->input.score, &sfr, revcmpl && Tflag,
		    NULL, 0);
	}

	if (sfr.score != fpo->input.score) {
		fprintf(stderr, "warning: score differs from input "
		    "file (read=\"%s\", genome=\"%s\")\n",
		    fpo->input.read, fpo->input.genome);
		fprintf(stderr, "         Are you using different S-W "
		    "parameters than before?\n");
	}

	fpo->output_normal = output_normal(read->name, contig->name, &sfr,
	    contig->sequence_len, (shrimp_mode == MODE_COLOUR_SPACE),
	    read->sequence, read->sequence_len, read->initbp, revcmpl,
	    Rflag);

	fpo->output_pretty = output_pretty(read->name, contig->name, &sfr,
	    contig->sequence, contig->sequence_len,
	    (shrimp_mode == MODE_COLOUR_SPACE), read->sequence,
	    read->sequence_len, read->initbp, revcmpl);

	assert(fpo->output_normal != NULL);
	assert(fpo->output_pretty != NULL);

	free(sfr.dbalign);
	free(sfr.qralign);
}

static void
print_alignments()
{
	struct fpo *fpo;
	char *fmt;

	fmt = output_format_line(Rflag);
	printf("%s%s\n", fmt, (seen_probs) ? " normodds pgenome pchance" : "");
	free(fmt);

	for (fpo = alignments_ordered; fpo != NULL; fpo = fpo->next_ordered) {
		if (fpo->output_normal == NULL || fpo->output_pretty == NULL) {
			fprintf(stderr, "warning: could not align read [%s] to contig"
			    "[%s] - missing contig file!\n", fpo->input.read,
			    fpo->input.genome);
			continue;
		}

		printf("%s", fpo->output_normal);
		if (INPUT_HAS_NORMODDS(&fpo->input))
			printf("\t%e", fpo->input.normodds);
		if (INPUT_HAS_PGENOME(&fpo->input))
			printf("\t%e", fpo->input.pgenome);
		if (INPUT_HAS_PCHANCE(&fpo->input))
			printf("\t%e", fpo->input.pchance);
		puts("\n");
		puts(fpo->output_pretty);
	}
}

static void
load_output_file(char *file)
{
	struct input inp;
	gzFile fp;
	struct fpo *fpo, *lastfpo = NULL;

	fp = gzopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "error: failed to open probcalc/rmapper output "
		    "file [%s]: %s\n", file, strerror(errno));
		exit(1);
	}

	while (input_parseline(fp, &inp)) {
		bool found;
		struct contig_ll *cll;

		if (INPUT_HAS_NORMODDS(&inp) || INPUT_HAS_PGENOME(&inp) ||
		    INPUT_HAS_PCHANCE(&inp))
			seen_probs = true;

		fpo = (struct fpo *)xmalloc(sizeof(*fpo));
		memset(fpo, 0, sizeof(*fpo));
		memcpy(&fpo->input, &inp, sizeof(fpo->input));

		/*
		 * Cache the contig name. We'll use it later when loading
		 * contigs.
		 */
		found = dynhash_find(contig_list, fpo->input.genome, NULL,
		    (void **)(void *)&cll);
		if (found) {
			assert(cll != NULL);
			if (INPUT_IS_REVCMPL(&inp)) {
				fpo->next_contig = cll->head_revcmpl;
				cll->head_revcmpl = fpo;
			} else {
				fpo->next_contig = cll->head;
				cll->head = fpo;
			}
			dynhash_remove(contig_list, fpo->input.genome, NULL, NULL);
		} else {
			cll = (struct contig_ll *)xmalloc(sizeof(*cll));
			memset(cll, 0, sizeof(*cll));
			if (INPUT_IS_REVCMPL(&inp))
				cll->head_revcmpl = fpo;
			else
				cll->head = fpo;
		}

		if (!dynhash_add(contig_list, fpo->input.genome, cll)) {
			fprintf(stderr, "error: failed to add read to contig "
			    "list - probably out of memory\n");
			exit(1);
		}

		if (lastfpo == NULL)
			alignments_ordered = fpo;
		else
			lastfpo->next_ordered = fpo;
		lastfpo = fpo;

		if (INPUT_IS_REVCMPL(&inp))
			nalignments_revcmpl++;
		else
			nalignments++;

		/*
		 * Cache the read name. We'll use it later when loading reads
		 * so that we only take in what we need.
		 */ 
		if (!dynhash_find(read_list, fpo->input.read, NULL, NULL)) {
			if (!dynhash_add(read_list, fpo->input.read, NULL)) {
				fprintf(stderr, "error: failed to add read to list - "
				    "probably out of memory\n");
				exit(1);
			}
		}

		longest_read_len = MAX(longest_read_len,
		   fpo->input.read_length);
	}

	gzclose(fp);
}

static void
load_genome_file(char *fpath, struct stat *sb, void *arg)
{
	fasta_t fasta;
	struct sequence *s;
	char *name, *seq;
	bool is_rna = false;

	/* shut up, icc */
	(void)sb;
	(void)arg;

	fasta = fasta_open(fpath, LETTER_SPACE);
	if (fasta == NULL) {
		fprintf(stderr, "error: failed to parse contig fasta file "
		    "[%s]\n", fpath);
		exit(1);
	}

	while (fasta_get_next(fasta, &name, &seq, &is_rna)) {
		bool found;
		struct contig_ll *cll;

		s = (struct sequence *)xmalloc(sizeof(*s));
		s->sequence = fasta_sequence_to_bitfield(fasta, seq);
		s->sequence_len = strlen(seq);
		s->name = name;
		s->initbp = -1;
		s->revcmpl = false;
		s->is_rna = is_rna;
		free(seq);

		/* add to our dynhash */
		found = dynhash_find(contig_list, s->name, NULL,
		    (void **)(void *)&cll);
		if (found) {
			struct fpo *fpo;

			assert(cll != NULL);

			/*
			 * Compute alignments for all inputs that required this contig.
			 */
			for (fpo = cll->head; fpo != NULL; fpo = fpo->next_contig)
				compute_alignment(fpo, s);

			reverse_complement(s->sequence, NULL, s->sequence_len, s->is_rna);
			s->revcmpl = true;

			for (fpo = cll->head_revcmpl; fpo != NULL; fpo = fpo->next_contig)
				compute_alignment(fpo, s);
		}

		ncontig_bases += s->sequence_len;
		free(s->sequence);
		free(s->name);
		free(s);
	}

	fasta_close(fasta);

	ncontig_files++;
}

static void
load_reads_file(char *fpath, struct stat *sb, void *arg)
{
	fasta_t fasta;
	struct sequence *s;
	char *name, *seq;
	int space;

	/* shut up, icc */
	(void)sb;
	(void)arg;

	if (shrimp_mode == MODE_COLOUR_SPACE)
		space = COLOUR_SPACE;
	else
		space = LETTER_SPACE;

	fasta = fasta_open(fpath, space);
	if (fasta == NULL) {
		fprintf(stderr, "error: failed to parse reads fasta file "
		    "[%s]\n", fpath);
		exit(1);
	}

	while (fasta_get_next(fasta, &name, &seq, NULL)) {
		void *key, *val;
		bool found;

		s = (struct sequence *)xmalloc(sizeof(*s));
		s->sequence = fasta_sequence_to_bitfield(fasta, seq);
		s->sequence_len = strlen(seq);
		s->name = name;
		if (shrimp_mode == MODE_COLOUR_SPACE) {
			s->initbp = fasta_get_initial_base(fasta, seq);
			s->sequence_len--;
			assert(s->sequence_len > 1);
		}
		s->revcmpl = false;
		free(seq);

		longest_read_len = MAX(longest_read_len, s->sequence_len);

		/*
		 * See if this read is in our dynhash (i.e. that we have an
		 * alignment that requires it.
		 */
		found = dynhash_find(read_list, s->name, &key, &val);
		if (found) {
			if (val != NULL) {
				fprintf(stderr, "error: read [%s] occurs multiple "
				    "times in the read input files\n", s->name);
				exit(1);
			}

			dynhash_remove(read_list, key, NULL, NULL);
			free(s->name);
			s->name = (char *)key;

			if (!dynhash_add(read_list, s->name, s)) {
				fprintf(stderr, "error: failed to add read to list - "
				    "probably out of memory\n");
				exit(1);
			}

			nread_bases += s->sequence_len;
		} else {
			free(s->sequence);
			free(s->name);
			free(s);
		}
	}

	fasta_close(fasta);

	nread_files++;
}

static uint32_t
keyhasher(void *k)
{

	return (hash_string((char *)k));
}

static int
keycomparer(void *a, void *b)
{

	return (strcmp((char *)a, (char *)b));
}

static void
usage(char *progname)
{
	char *slash;

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s [parameters] [options] shrimp_output_file "
	    "genome_dir|genome_file reads_dir|read_file\n", progname);

	fprintf(stderr, "Parameters:\n");

	fprintf(stderr,
	    "    -m    S-W Match Value                         (default: %d)\n",
	    DEF_MATCH_VALUE);

	fprintf(stderr,
	    "    -i    S-W Mismatch Value                      (default: %d)\n",
	    DEF_MISMATCH_VALUE);

	fprintf(stderr,
	    "    -g    S-W Gap Open Penalty (Reference)        (default: %d)\n",
	    a_gap_open);

	fprintf(stderr,
	    "    -q    S-W Gap Open Penalty (Query)            (default: %d)\n",
	    b_gap_open);

	fprintf(stderr,
	    "    -e    S-W Gap Extend Penalty (Reference)      (default: %d)\n",
	    a_gap_extend);

	fprintf(stderr,
	    "    -f    S-W Gap Extend Penalty (Query)          (default: %d)\n",
	    b_gap_extend);

	if (shrimp_mode == MODE_COLOUR_SPACE) {
		fprintf(stderr,
		    "    -x    S-W Crossover Penalty                   ("
		    "default: %d)\n", DEF_XOVER_PENALTY);
	}

	fprintf(stderr, "\nOptions:\n");

	fprintf(stderr,
	    "    -R    Print Reads in Output (if in input)     (default: "
	    "disabled)\n");

	exit(1);
}

int
main(int argc, char **argv)
{
	char *fpout, *readsdir, *genomedir, *progname;
	char const * optstr;
	bool a_gap_open_set, b_gap_open_set;
	bool a_gap_extend_set, b_gap_extend_set;
	int ch, ret;

	a_gap_open_set = b_gap_open_set = false;
	a_gap_extend_set = b_gap_extend_set = false;

	set_mode_from_argv(argv);

	if (shrimp_mode == MODE_HELICOS_SPACE) {
		match_value	= DEF_MATCH_VALUE_DAG;
		mismatch_value	= DEF_MISMATCH_VALUE_DAG;
		a_gap_open	= DEF_A_GAP_OPEN_DAG;
		b_gap_open	= DEF_B_GAP_OPEN_DAG;
		a_gap_extend	= DEF_A_GAP_EXTEND_DAG;
		b_gap_extend	= DEF_B_GAP_EXTEND_DAG;
	}

	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");
	fprintf(stderr, "prettyprint: %s SPACE.\nSHRiMP %s [%s]\n",
	    get_mode_string(), SHRIMP_VERSION_STRING, get_compiler());
	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");

	progname = argv[0];

	if (shrimp_mode == MODE_COLOUR_SPACE)
		optstr = "m:i:g:e:x:RT";
	else
		optstr = "m:i:g:e:RT";

	while ((ch = getopt(argc, argv, optstr)) != -1) {
		switch (ch) {
		case 'm':
			match_value = atoi(optarg);
			break;
		case 'i':
			mismatch_value = atoi(optarg);
			break;
		case 'g':
			a_gap_open_set = true;
			a_gap_open = atoi(optarg);
			break;
		case 'q':
			b_gap_open_set = true;
			b_gap_open = atoi(optarg);
			break;
		case 'e':
			a_gap_extend_set = true;
			a_gap_extend = atoi(optarg);
			break;
		case 'f':
			b_gap_extend_set = true;
			b_gap_extend = atoi(optarg);
			break;
		case 'x':
			assert(shrimp_mode == MODE_COLOUR_SPACE);
			xover_penalty = atoi(optarg);
			break;
		case 'R':
			Rflag = true;
			break;
		case 'T':
			Tflag = true;
			break;
		default:
			usage(progname);
		}
	}
	argc -= optind;
	argv += optind;

	if (argc != 3)
		usage(progname);

	if ((a_gap_open_set && !b_gap_open_set)
	    || (a_gap_extend_set && !b_gap_extend_set))
		fputc('\n', stderr);
	if (a_gap_open_set && !b_gap_open_set) {
		fprintf(stderr, "Notice: Gap open penalty set for reference but not query; assuming symmetry.\n");
		b_gap_open = a_gap_open;
	}
	if (a_gap_extend_set && !b_gap_extend_set) {
		fprintf(stderr, "Notice: Gap extend penalty set for reference but not query; assuming symmetry.\n");
		b_gap_extend = a_gap_extend;
	}
	if ((a_gap_open_set && !b_gap_open_set)
	    || (a_gap_extend_set && !b_gap_extend_set))
		fputc('\n', stderr);

	fprintf(stderr, "Settings:\n");
	fprintf(stderr, "    S-W Match Value:                  %d\n", match_value);
	fprintf(stderr, "    S-W Mismatch Value:               %d\n", mismatch_value);
	fprintf(stderr, "    S-W Gap Open Penalty (Ref):       %d\n", a_gap_open);
	fprintf(stderr, "    S-W Gap Open Penalty (Qry):       %d\n", b_gap_open);
	fprintf(stderr, "    S-W Gap Extend Penalty (Ref):     %d\n", a_gap_extend);
	fprintf(stderr, "    S-W Gap Extend Penalty (Qry):     %d\n", b_gap_extend);
	
	if (shrimp_mode == MODE_COLOUR_SPACE) {
		fprintf(stderr, "    S-W Crossover Penalty:            %d\n",
		    xover_penalty);
	}
	fputc('\n', stderr);

	read_list   = dynhash_create(keyhasher, keycomparer);
	contig_list = dynhash_create(keyhasher, keycomparer);
	if (read_list == NULL || contig_list == NULL) {
		fprintf(stderr, "error: failed to allocate read and contig "
		    "lists\n");
		exit(1);
	}

	fpout     = argv[0];
	genomedir = argv[1];
	readsdir  = argv[2];

	fprintf(stderr, "Loading shrimp output file...\n");
	load_output_file(fpout);
	fprintf(stderr, "Loaded %" PRIu64 " alignments from shrimp output (%u unique reads)\n",
	    nalignments + nalignments_revcmpl, dynhash_count(read_list));

	fprintf(stderr, "Loading read file(s)...\n");
	file_iterator(readsdir, load_reads_file, NULL);
	fprintf(stderr, "Loaded %u %s reads in %u file(s) (%" PRIu64 " total bases)\n",
	    dynhash_count(read_list),
	    (shrimp_mode == MODE_COLOUR_SPACE) ? "colourspace" : "letterspace", nread_files,
	    nread_bases);

	if (shrimp_mode == MODE_COLOUR_SPACE) {
/* XXX - a vs. b gap */
		ret = sw_full_cs_setup(longest_read_len * 10, longest_read_len,
		    a_gap_open, a_gap_extend, match_value, mismatch_value,
		    xover_penalty, false, -1);
	} else {
		ret = sw_full_ls_setup(longest_read_len * 10, longest_read_len,
		    a_gap_open, a_gap_extend, b_gap_open, b_gap_extend, match_value,
		    mismatch_value, false, -1);
	}
	if (ret) {
		fprintf(stderr, "failed to initialise scalar Smith-Waterman "
		    "(%s)\n", strerror(errno));
		exit(1);
	}

	fprintf(stderr, "Loading genome contig file(s)...\n");
	file_iterator(genomedir, load_genome_file, NULL);
	fprintf(stderr, "Loaded %u genomic contigs in %u file(s) (%" PRIu64 " total "
	    "bases)\n", dynhash_count(contig_list), ncontig_files, ncontig_bases);

	print_alignments();

	return (0);
}
