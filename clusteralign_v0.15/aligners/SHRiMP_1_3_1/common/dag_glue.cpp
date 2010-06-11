/*	$Id: dag_glue.cpp 383 2009-09-30 22:47:56Z matei $	*/

#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <string>
#include <climits>
#include <fstream>
#include <zlib.h>

#include "../common/dag_align.h"
#include "../common/dag_glue.h"
//extern "C" {
#include "../common/util.h"
//}

using namespace std;

static uint64_t aligner_ticks;
static uint64_t aligner_invocations;
static bool	is_set_up = false;		/* be extra careful this time*/
static uint64_t kmers_invocations;
static uint64_t kmers_ticks;
static uint64_t kmers_total_kmers;

/*
 * This file contains a C wrapper interface for the C++ DAG alignment
 * implementation.
 */

void
dag_setup(int read_match, int read_mismatch, int read_gap, int dag_match,
    int dag_snp, int dag_half_match, int dag_neither_match, int dag_match_deletion,
    int dag_mismatch_deletion, int dag_error_insertion)
{
	is_set_up = true;

	Column::setscore(read_match, read_gap, read_mismatch, dag_match,
	    dag_snp, dag_half_match, dag_neither_match, dag_match_deletion,
	    dag_mismatch_deletion, dag_error_insertion);
}

dag_cookie_t
dag_build_kmer_graph(char *read1, char *read2, int epsilon)
{
	uint64_t before = rdtsc();

	assert(is_set_up);

	Graph g1(read1);
	Graph g2(read2);
	CrossProduct cp(g1, g2);
	Graph *kg = new Graph(cp, epsilon);
	dag_cookie_t dct = (dag_cookie_t)xmalloc(sizeof(*dct));
	dct->dag = (void *)kg;
	dct->loops_added = false;

	kmers_ticks += rdtsc() - before;

	return (dct);
}

void
dag_free_kmer_graph(dag_cookie_t dct)
{
	Graph *g = (Graph *)dct->dag;
	uint64_t before = rdtsc();

	assert(is_set_up);

	delete g;
	free(dct);

	kmers_ticks += rdtsc() - before;
}

char **
dag_get_kmers(dag_cookie_t dct, int len)
{
	int i;
	char **list;
	kmersInfo *ki;
	set <string>::iterator setit;
	Graph *kg = (Graph *)dct->dag;
	uint64_t before = rdtsc();

	assert(is_set_up);
	kmers_invocations++;

	if ((len - 1) < 0)
		return (NULL);

	ki = kg->getkmers(len);

	if ((len - 1) >= (int)(ki->kmers).size()) {
		delete ki;
		return (NULL);
	}

	i = 1;
	list = (char **)xmalloc(sizeof(char *));
	list[0] = NULL;

	for (setit = (ki->kmers).at(len - 1).begin(); setit != (ki->kmers).at(len - 1).end(); setit++) {
		list = (char **)xrealloc(list, (i + 1) * sizeof(char *));
		list[i - 1] = xstrdup(setit->c_str());
		list[i++] = NULL;
	}

	kmers_total_kmers += i - 1;

	delete ki;

	kmers_ticks += rdtsc() - before;

	return (list);
}

void
dag_free_kmers(char **kmers)
{
	int i;

	assert(is_set_up);

	for (i = 0; kmers[i] != NULL; i++)
		free(kmers[i]);
	free(kmers);
}

struct dag_alignment *
dag_build_alignment(char *genome, dag_cookie_t dct)
{
	Graph *kg;
	Alignment *al;
	SmallCrossProduct *scp;
	SmallAlignNode *bestend;
	Node *bestendap,  *bestendbp;
	struct dag_alignment *dap;
	uint64_t before;

	assert(is_set_up);

	before = rdtsc();
	aligner_invocations++;

	dap = (struct dag_alignment *)xmalloc(sizeof(*dap));
	kg = (Graph *)dct->dag;
	Graph g(genome);
	g.AddSelfLoops();
	if (!dct->loops_added) {
		kg->AddSelfLoops();
		dct->loops_added = true;
	}
	scp = new SmallCrossProduct(g, *kg);
	dap->score = scp->DijkstraForward(&bestend, &bestendap, &bestendbp);
	al = scp->GetBestPath(bestend, bestendap, bestendbp);
	dap->start_index = al->seqstart;
	dap->end_index = al->seqend;
	dap->read1 = (char *)xstrdup(al->read1.c_str());
	dap->read2 = (char *)xstrdup(al->read2.c_str());
	dap->sequence = (char *)xstrdup(al->seq.c_str());
	
	delete scp;
	delete al;

	aligner_ticks += rdtsc() - before;

	return (dap);
}

void
dag_free_alignment(struct dag_alignment *dap)
{
	uint64_t before;

	assert(is_set_up);

	before = rdtsc();
	free(dap->read1);
	free(dap->read2);
	free(dap->sequence);
	free(dap);
	aligner_ticks += rdtsc() - before;
}

struct dag_statistics *
dag_get_statistics()
{
	struct dag_statistics *dsp;

	assert(is_set_up);

	dsp = (struct dag_statistics *)xmalloc(sizeof(*dsp));
	dsp->aligner_invocations = aligner_invocations;
	dsp->aligner_seconds = (double)aligner_ticks / cpuhz();
	dsp->kmers_invocations = kmers_invocations;
	dsp->kmers_total_kmers = kmers_total_kmers;
	dsp->kmers_seconds = (double)kmers_ticks / cpuhz();

	return (dsp);
}

void
dag_free_statistics(struct dag_statistics *dsp)
{
	assert(is_set_up);
	free(dsp);
}
