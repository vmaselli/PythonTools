/*	$Id: dag_glue.h 383 2009-09-30 22:47:56Z matei $	*/

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _dag_cookie_t {
	void *dag;
	bool  loops_added;
} * dag_cookie_t;

struct dag_alignment {
	int score;
	int start_index;
	int end_index;
	
	char *sequence;
	char *read1;
	char *read2;
};

struct dag_statistics {
	double   aligner_seconds;
	uint64_t aligner_invocations;
	uint64_t kmers_invocations;
	uint64_t kmers_total_kmers;
	double   kmers_seconds;
};

void		       dag_setup(int, int, int, int, int, int, int, int, int, int);
dag_cookie_t	       dag_build_kmer_graph(char *, char *, int);
void                   dag_free_kmer_graph(dag_cookie_t);
char                 **dag_get_kmers(dag_cookie_t, int);
void                   dag_free_kmers(char **);
extern struct dag_alignment  *dag_build_alignment(char *, dag_cookie_t);
void		       dag_free_alignment(struct dag_alignment *);
struct dag_statistics *dag_get_statistics(void);
void		       dag_free_statistics(struct dag_statistics *);

#ifdef __cplusplus
}
#endif
