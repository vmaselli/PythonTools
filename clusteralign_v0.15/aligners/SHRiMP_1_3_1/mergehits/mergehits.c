/*	$Id: mergehits.c 348 2009-06-16 23:26:27Z rumble $	*/

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

#include "../common/fasta.h"
#include "../common/dynhash.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/sw-full-ls.h"
#include "../common/input.h"
#include "../common/output.h"
#include "../common/util.h"
#include "../common/version.h"

/* alignment list member */
struct alist_item {
	struct input       input;	/* input fields */
	struct alist_item *next;	/* next alignment in list */
};

/* cached data associated with a read */
struct sequence {
	char	 *name;
	uint32_t *sequence;
	uint32_t  sequence_len;
	int	  initbp;		      /* initial base pair (cs only) */
	struct    alist_item *alignments;     /* start of alignment list */
	struct    alist_item *last_alignment; /* end of alignment list */
	struct    sequence *next;	      /* next sequence in list */
};

static struct sequence *read_list;	/* list of read names and alignments */
static dynhash_t read_index;		/* name index for read_list */
static uint64_t nalignments;		/* number of alignments */

static bool Rflag = false;		/* don't output read sequence */

#define MAX_READ_LEN	5000		/* ridiculously high upper bound */

static void
load_output_file(char *file)
{
	struct input inp;
	gzFile fp;
	struct alist_item *alignment;
	struct sequence *seq, *last_read;

	fp = gzopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "error: failed to open shrimp output file "
		    "[%s]: %s\n", file, strerror(errno));
		exit(1);
	}

	assert(read_list == NULL);
	last_read = NULL;
	while (input_parseline(fp, &inp)) {
		alignment = (struct alist_item *)xmalloc(sizeof(*alignment));
		memcpy(&alignment->input, &inp, sizeof(alignment->input));
		alignment->next = NULL;

		if (dynhash_find(read_index, inp.read, NULL,
		    (void **)(void *)&seq)) {
			/*
			 * we've already seen an alignment for the same
			 * read, so add this one to that read's list
			 */
			assert(seq != NULL);
			assert(seq->last_alignment != NULL);
			assert(seq->last_alignment->next == NULL);
			seq->last_alignment->next = alignment;
			seq->last_alignment = alignment;
		} else {
			/*
		         * this alignment is for a read we haven't seen
			 * yet, so allocate a new list
			 */
			seq = (struct sequence *)xmalloc(sizeof(*seq));
			seq->name = xstrdup(inp.read);
			seq->sequence = NULL;
			seq->sequence_len = 0;
			seq->initbp = 0;
			seq->alignments = alignment;
			seq->last_alignment = alignment;
			seq->next = NULL;
			if (last_read == NULL) {
				read_list = seq;
			} else {
				last_read->next = seq;
			}
			last_read = seq;

			if (dynhash_add(read_index, seq->name, seq) == false) {
				fprintf(stderr, "error: failed to add "
					"read to index - "
					"probably out of memory\n");
				exit(1);
			}
			nalignments++;
		}
	}

	gzclose(fp);
}

static void
store_sequence(char *name, uint32_t *read, uint32_t read_len, int initbp)
{
	struct sequence *seq;

	if (!dynhash_find(read_index, name, NULL, (void **)(void *)&seq)) {
		/* this read has no alignments so skip it */
		return;
	}
	assert(seq != NULL);
	seq->sequence = read;
	seq->sequence_len = read_len;
	seq->initbp = initbp;
}

static void
load_reads_file(char *fpath, struct stat *sb, void *arg)
{
	fasta_t fasta;
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
		uint32_t *bf = fasta_sequence_to_bitfield(fasta, seq);
		uint32_t seqlen = strlen(seq);
		int initbp = -1;

		assert(seqlen > 0);
		if (shrimp_mode == MODE_COLOUR_SPACE) {
			seqlen -= 1;
			initbp = fasta_get_initial_base(fasta, seq);
		}
		assert(seqlen > 0);

		store_sequence(name, bf, seqlen, initbp);
		free(name);
		free(seq);
	}

	fasta_close(fasta);
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
print_alignments()
{
	struct sequence   *seq;
	struct alist_item *alignment;
	char              *read_str;
	int32_t		   genome_align_start;
	uint32_t	   genome_align_size;
	uint32_t	   read_align_size;

	for (seq = read_list; seq != NULL; seq = seq->next) {
		printf(">%s", seq->name);
		for (alignment = seq->alignments; alignment != NULL;
		     alignment = alignment->next) {
			/*
			 * NB: internally 0 is first position, output uses 1.
			 *     adjust.
			 */
			genome_align_start =
			    INPUT_IS_REVCMPL(&alignment->input) ? 
			    -(alignment->input.genome_start + 1) :
			alignment->input.genome_start + 1;
			genome_align_size = alignment->input.genome_end -
			    alignment->input.genome_start + 1;
			read_align_size = alignment->input.read_end -
			    alignment->input.read_start + 1;
			printf(",%s.%d.%u.%hu.%u.%d", alignment->input.genome,
			    genome_align_start, genome_align_size,
			    alignment->input.read_start + 1, read_align_size,
			    alignment->input.score);
		}
		printf("\n");
		if (seq->sequence_len > 0) {
			read_str = readtostr(seq->sequence, seq->sequence_len,
			    (shrimp_mode == MODE_COLOUR_SPACE), seq->initbp);
			puts(read_str);
		}
	}
}

static void
usage(char *progname)
{
	char *slash;

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s [-R reads_dir] shrimp_output_file\n",
	    progname);

	exit(1);
}

int
main(int argc, char **argv)
{
	char *fpout, *readsdir, *progname;
	struct stat sb;
	int ch;

	set_mode_from_argv(argv);

	/* shut up, gcc */
	readsdir = NULL;

	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");
	fprintf(stderr, "mergehits: %s.\nSHRiMP %s [%s]\n",
	    get_mode_string(), SHRIMP_VERSION_STRING,
	    get_compiler());
	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");

	progname = argv[0];

	while ((ch = getopt(argc, argv, "R:")) != -1) {
		switch (ch) {
		case 'R':
			Rflag = true;
			readsdir = optarg;
			break;
		default:
			usage(progname);
		}
	}
	argc -= optind;
	argv += optind;

	if (argc != 1)
		usage(progname);

	read_index   = dynhash_create(keyhasher, keycomparer);
	if (read_index == NULL) {
		fprintf(stderr, "error: failed to allocate read index\n");
		exit(1);
	}

	fpout = argv[0];

	xstat(fpout, &sb);
	if (!S_ISREG(sb.st_mode)) {
		fprintf(stderr, "error: cannot open shrimp output file [%s]: "
		    "not a regular file or symlink to a regular file.\n",fpout);
		exit(1);
	}

	fprintf(stderr, "Loading shrimp output file...\n");
	load_output_file(fpout);
	fprintf(stderr, "Loaded %" PRIu64 " alignments from shrimp output\n",
	    nalignments);

	if (Rflag) {
		fprintf(stderr, "Loading read file(s)...\n");
		file_iterator(readsdir, load_reads_file, NULL);
	}

	print_alignments();

	return (0);
}
