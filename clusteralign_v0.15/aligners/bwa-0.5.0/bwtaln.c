#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"

#ifdef HAVE_PTHREAD
#define THREAD_BLOCK_SIZE 1024
#include <pthread.h>
static pthread_mutex_t g_seq_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

gap_opt_t *gap_init_opt()
{
	gap_opt_t *o;
	o = (gap_opt_t*)calloc(1, sizeof(gap_opt_t));
	/* IMPORTANT: s_mm*10 should be about the average base error
	   rate. Voilating this requirement will break pairing! */
	o->s_mm = 3; o->s_gapo = 11; o->s_gape = 4;
	o->max_diff = -1; o->max_gapo = 1; o->max_gape = 6;
	o->indel_end_skip = 5; o->max_del_occ = 10; o->max_entries = 2000000;
	o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
	o->seed_len = 0x7fffffff; o->max_seed_diff = 2;
	o->fnr = 0.04;
	o->n_threads = 1;
	o->max_top2 = 30;
	return o;
}

int bwa_cal_maxdiff(int l, double err, double thres)
{
	double elambda = exp(-l * err);
	double sum, y = 1.0;
	int k, x = 1;
	for (k = 1, sum = elambda; k < 1000; ++k) {
		y *= l * err;
		x *= k;
		sum += elambda * y / x;
		if (1.0 - sum < thres) return k;
	}
	return 2;
}

// width must be filled as zero
static int bwt_cal_width(const bwt_t *rbwt, int len, const ubyte_t *str, bwt_width_t *width)
{
	bwtint_t k, l, ok, ol;
	int i, bid;
	bid = 0;
	k = 0; l = rbwt->seq_len;
	for (i = 0; i < len; ++i) {
		ubyte_t c = str[i];
		if (c < 4) {
			bwt_2occ(rbwt, k - 1, l, c, &ok, &ol);
			k = rbwt->L2[c] + ok + 1;
			l = rbwt->L2[c] + ol;
		}
		if (k > l || c > 3) { // then restart
			k = 0;
			l = rbwt->seq_len;
			++bid;
		}
		width[i].w = l - k + 1;
		width[i].bid = bid;
	}
	width[len].w = 0;
	width[len].bid = ++bid;
	return bid;
}

static void bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt[2], int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt)
{
	int i, max_l = 0, max_len;
	gap_stack_t *stack;
	bwt_width_t *w[2], *seed_w[2];
	const ubyte_t *seq[2];
	gap_opt_t local_opt = *opt;

	// initiate priority stack
	for (i = max_len = 0; i != n_seqs; ++i)
		if (seqs[i].len > max_len) max_len = seqs[i].len;
	if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
	if (local_opt.max_diff < local_opt.max_gapo) local_opt.max_gapo = local_opt.max_diff;
	stack = gap_init_stack(opt->max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);

	seed_w[0] = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
	seed_w[1] = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
	w[0] = w[1] = 0;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;
#ifdef HAVE_PTHREAD
		if (opt->n_threads > 1) {
			pthread_mutex_lock(&g_seq_lock);
			if (p->tid < 0) { // unassigned
				int j;
				for (j = i; j < n_seqs && j < i + THREAD_BLOCK_SIZE; ++j)
					seqs[j].tid = tid;
			} else if (p->tid != tid) {
				pthread_mutex_unlock(&g_seq_lock);
				continue;
			}
			pthread_mutex_unlock(&g_seq_lock);
		}
#endif
		p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
		seq[0] = p->seq; seq[1] = p->rseq;
		if (max_l < p->len) {
			max_l = p->len;
			w[0] = (bwt_width_t*)calloc(max_l + 1, sizeof(bwt_width_t));
			w[1] = (bwt_width_t*)calloc(max_l + 1, sizeof(bwt_width_t));
		}
		bwt_cal_width(bwt[0], p->len, seq[0], w[0]);
		bwt_cal_width(bwt[1], p->len, seq[1], w[1]);
		if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
		if (opt->seed_len >= p->len) local_opt.seed_len = 0x7fffffff;
		if (p->len > opt->seed_len) {
			bwt_cal_width(bwt[0], opt->seed_len, seq[0] + (p->len - opt->seed_len), seed_w[0]);
			bwt_cal_width(bwt[1], opt->seed_len, seq[1] + (p->len - opt->seed_len), seed_w[1]);
		}
		// core function
		p->aln = bwt_match_gap(bwt, p->len, seq, w, p->len <= opt->seed_len? 0 : seed_w, &local_opt, &p->n_aln, stack);
		// store the alignment
		free(p->name); free(p->seq); free(p->rseq); free(p->qual);
		p->name = 0; p->seq = p->rseq = p->qual = 0;
	}
	free(seed_w[0]); free(seed_w[1]);
	free(w[0]); free(w[1]);
	gap_destroy_stack(stack);
}

#ifdef HAVE_PTHREAD
typedef struct {
	int tid;
	bwt_t *bwt[2];
	int n_seqs;
	bwa_seq_t *seqs;
	const gap_opt_t *opt;
} thread_aux_t;

static void *worker(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	bwa_cal_sa_reg_gap(d->tid, d->bwt, d->n_seqs, d->seqs, d->opt);
	return 0;
}
#endif

void bwa_aln_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt)
{
	int i, n_seqs, tot_seqs = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t;
	bwt_t *bwt[2];

	// initialization
	ks = bwa_seq_open(fn_fa);

	{ // load BWT
		char *str = (char*)calloc(strlen(prefix) + 10, 1);
		strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
		strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);
		free(str);
	}

	// core loop
	fwrite(opt, sizeof(gap_opt_t), 1, stdout);
	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode & BWA_MODE_COMPREAD)) != 0) {
		tot_seqs += n_seqs;
		t = clock();

		fprintf(stderr, "[bwa_aln_core] calculate SA coordinate... ");

#ifdef HAVE_PTHREAD
		if (opt->n_threads <= 1) { // no multi-threading at all
			bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
		} else {
			pthread_t *tid;
			pthread_attr_t attr;
			thread_aux_t *data;
			int j;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
			tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
			for (j = 0; j < opt->n_threads; ++j) {
				data[j].tid = j; data[j].bwt[0] = bwt[0]; data[j].bwt[1] = bwt[1];
				data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = opt;
				pthread_create(&tid[j], &attr, worker, data + j);
			}
			for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
			free(data); free(tid);
		}
#else
		bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
#endif

		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		t = clock();
		fprintf(stderr, "[bwa_aln_core] write to the disk... ");
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
			fwrite(&p->n_aln, 4, 1, stdout);
			if (p->n_aln) fwrite(p->aln, sizeof(bwt_aln1_t), p->n_aln, stdout);
		}
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		bwa_free_read_seq(n_seqs, seqs);
		fprintf(stderr, "[bwa_aln_core] %d sequences have been processed.\n", tot_seqs);
	}

	// destroy
	bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
	bwa_seq_close(ks);
}

int bwa_aln(int argc, char *argv[])
{
	int c, opte = -1;
	gap_opt_t *opt;

	opt = gap_init_opt();
	while ((c = getopt(argc, argv, "n:o:e:i:d:l:k:cLR:m:t:NM:O:E:")) >= 0) {
		switch (c) {
		case 'n':
			if (strstr(optarg, ".")) opt->fnr = atof(optarg), opt->max_diff = -1;
			else opt->max_diff = atoi(optarg), opt->fnr = -1.0;
			break;
		case 'o': opt->max_gapo = atoi(optarg); break;
		case 'e': opte = atoi(optarg); break;
		case 'M': opt->s_mm = atoi(optarg); break;
		case 'O': opt->s_gapo = atoi(optarg); break;
		case 'E': opt->s_gape = atoi(optarg); break;
		case 'd': opt->max_del_occ = atoi(optarg); break;
		case 'i': opt->indel_end_skip = atoi(optarg); break;
		case 'l': opt->seed_len = atoi(optarg); break;
		case 'k': opt->max_seed_diff = atoi(optarg); break;
		case 'm': opt->max_entries = atoi(optarg); break;
		case 't': opt->n_threads = atoi(optarg); break;
		case 'L': opt->mode |= BWA_MODE_LOGGAP; break;
		case 'R': opt->max_top2 = atoi(optarg); break;
		case 'c': opt->mode &= ~BWA_MODE_COMPREAD; break;
		case 'N': opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0x7fffffff; break;
		default: return 1;
		}
	}
	if (opte > 0) {
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa aln [options] <prefix> <in.fq>\n\n");
		fprintf(stderr, "Options: -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n",
				BWA_AVG_ERR, opt->fnr);
		fprintf(stderr, "         -o INT    maximum number or fraction of gap opens [%d]\n", opt->max_gapo);
		fprintf(stderr, "         -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
		fprintf(stderr, "         -i INT    do not put an indel within INT bp towards the ends [%d]\n", opt->indel_end_skip);
		fprintf(stderr, "         -d INT    maximum occurrences for extending a long deletion [%d]\n", opt->max_del_occ);
		fprintf(stderr, "         -l INT    seed length [%d]\n", opt->seed_len);
		fprintf(stderr, "         -k INT    maximum differences in the seed [%d]\n", opt->max_seed_diff);
		fprintf(stderr, "         -m INT    maximum entries in the queue [%d]\n", opt->max_entries);
		fprintf(stderr, "         -t INT    number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "         -M INT    mismatch penalty [%d]\n", opt->s_mm);
		fprintf(stderr, "         -O INT    gap open penalty [%d]\n", opt->s_gapo);
		fprintf(stderr, "         -E INT    gap extension penalty [%d]\n", opt->s_gape);
		fprintf(stderr, "         -R INT    stop searching when there are >INT equally best hits [%d]\n", opt->max_top2);
		fprintf(stderr, "         -c        input sequences are in the color space\n");
		fprintf(stderr, "         -L        log-scaled gap penalty for long deletions\n");
		fprintf(stderr, "         -N        non-iterative mode: search for all n-difference hits (slooow)\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (opt->fnr > 0.0) {
		int i, k;
		for (i = 17, k = 0; i <= 250; ++i) {
			int l = bwa_cal_maxdiff(i, BWA_AVG_ERR, opt->fnr);
			if (l != k) fprintf(stderr, "[bwa_aln] %dbp reads: max_diff = %d\n", i, l);
			k = l;
		}
	}
	bwa_aln_core(argv[optind], argv[optind+1], opt);
	free(opt);
	return 0;
}
