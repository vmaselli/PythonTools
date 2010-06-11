#include <zlib.h>
#include "bwtaln.h"
#include "utils.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

extern unsigned char nst_nt4_table[256];

struct __bwa_seqio_t {
	kseq_t *ks;
};

bwa_seqio_t *bwa_seq_open(const char *fn)
{
	gzFile fp;
	bwa_seqio_t *bs;
	bs = (bwa_seqio_t*)calloc(1, sizeof(bwa_seqio_t));
	fp = xzopen(fn, "r");
	bs->ks = kseq_init(fp);
	return bs;
}

void bwa_seq_close(bwa_seqio_t *bs)
{
	if (bs == 0) return;
	gzclose(bs->ks->f->f);
	kseq_destroy(bs->ks);
	free(bs);
}

void seq_reverse(int len, ubyte_t *seq, int is_comp)
{
	int i;
	if (is_comp) {
		for (i = 0; i < len>>1; ++i) {
			char tmp = seq[len-1-i];
			if (tmp < 4) tmp = 3 - tmp;
			seq[len-1-i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
			seq[i] = tmp;
		}
		if (len&1) seq[i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
	} else {
		for (i = 0; i < len>>1; ++i) {
			char tmp = seq[len-1-i];
			seq[len-1-i] = seq[i]; seq[i] = tmp;
		}
	}
}

bwa_seq_t *bwa_read_seq(bwa_seqio_t *bs, int n_needed, int *n, int is_comp)
{
	bwa_seq_t *seqs, *p;
	kseq_t *seq = bs->ks;
	int n_seqs, l, i;

	n_seqs = 0;
	seqs = (bwa_seq_t*)calloc(n_needed, sizeof(bwa_seq_t));
	while ((l = kseq_read(seq)) >= 0) {
		p = &seqs[n_seqs++];
		p->tid = -1; // no assigned to a thread
		p->qual = 0;
		p->len = l;
		p->seq = (ubyte_t*)calloc(p->len, 1);
		for (i = 0; i != p->len; ++i)
			p->seq[i] = nst_nt4_table[(int)seq->seq.s[i]];
		p->rseq = (ubyte_t*)calloc(p->len, 1);
		memcpy(p->rseq, p->seq, p->len);
		seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
		seq_reverse(p->len, p->rseq, is_comp);
		p->name = strdup((const char*)seq->name.s);
		{ // trim /[12]$
			int t = strlen(p->name);
			if (t > 2 && p->name[t-2] == '/' && (p->name[t-1] == '1' || p->name[t-1] == '2')) p->name[t-2] = '\0';
		}
		if (seq->qual.l) // copy quality
			p->qual = (ubyte_t*)strdup((char*)seq->qual.s);
		if (n_seqs == n_needed) break;
	}
	*n = n_seqs;
	if (n_seqs == 0) {
		free(seqs);
		return 0;
	}
	return seqs;
}

void bwa_free_read_seq(int n_seqs, bwa_seq_t *seqs)
{
	int i;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;
		free(p->name);
		free(p->seq); free(p->rseq); free(p->qual); free(p->aln); free(p->md);
		free(p->cigar);
	}
	free(seqs);
}
