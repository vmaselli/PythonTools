/*	$Id: dynhash.h 176 2007-12-31 17:09:10Z rumble $	*/

#ifndef _DYNHASH_H_
#define _DYNHASH_H_

#define DYNHASH_INIT_LENGTH	1024
#define DYNHASH_BULGE_FACTOR	2

struct dynhash_entry {
	void		     *key;
	void		     *val;
	struct dynhash_entry *next;
};

struct dynhash {
	struct dynhash_entry  **table;
	uint32_t		table_count;
	uint32_t		table_length;

	uint32_t	       (*hashfn)(void *);
	int		       (*keycmp)(void *, void *);

	bool			iterating;
};

typedef struct dynhash * dynhash_t;

dynhash_t dynhash_create(uint32_t (*)(void *), int (*)(void *, void*));
void	  dynhash_destroy(dynhash_t);
bool	  dynhash_find(dynhash_t, void *, void **, void **);
bool	  dynhash_add(dynhash_t, void *, void *);
bool	  dynhash_remove(dynhash_t, void *, void **, void **);
uint32_t  dynhash_count(dynhash_t);
void	  dynhash_iterate(dynhash_t, void (*)(void *, void *, void *), void *);

#endif /* !_DYNHASH_H_ */
