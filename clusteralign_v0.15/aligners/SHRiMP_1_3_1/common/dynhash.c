/*	$Id: dynhash.c 245 2008-06-06 18:24:28Z rumble $	*/

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "../common/dynhash.h"
#include "../common/util.h"

/*
 * Safety assertion: we don't want to allow any changes to internal state
 * while iterating.
 */
#define ASSERT_NOT_ITERATING(dh)					\
	do {								\
		if ((dh)->iterating) {					\
			fprintf(stderr, "error: %s while iterating\n",	\
			    __func__);					\
			abort();					\
		}							\
	} while (0)

/*
 * Expand and rehash our dynhash. If we fail to allocate for whatever reason,
 * just return. In most cases this means the application will die somewhere
 * else, possibly when dynhash_add cannot allocate a dhe and returns false...
 */
static void
dynhash_expand(dynhash_t dh)
{
	struct dynhash_entry **new_table, *dhe;
	uint32_t new_length, i, idx;

	assert(dh != NULL);
	assert(dh->table != NULL);

	ASSERT_NOT_ITERATING(dh);

	new_length = dh->table_length * DYNHASH_BULGE_FACTOR;

	/* don't bother if we've overflowed */
	if (new_length <= dh->table_length)
		return;

	new_table = (struct dynhash_entry **)
	    malloc(new_length * sizeof(struct dynhash_entry *));
	if (new_table == NULL)
		return;

	memset(new_table, 0, new_length * sizeof(struct dynhash_entry *));

	for (i = 0; i < dh->table_length; i++) {
		dhe = dh->table[i];
		while (dhe != NULL) {
			struct dynhash_entry *dhe_tmp;

			dhe_tmp = dhe->next;

			idx = dh->hashfn(dhe->key) % new_length;
			dhe->next = new_table[idx];
			new_table[idx] = dhe;

			dhe = dhe_tmp;
		}
	}

	free(dh->table);
	dh->table = new_table;
	dh->table_length = new_length;
}

/*
 * Create a new dynhash structure.
 */
dynhash_t
dynhash_create(uint32_t (*hashfn)(void *), int (*keycmp)(void *, void *))
{
	dynhash_t dh;

	if (hashfn == NULL || keycmp == NULL)
		return (NULL);

	dh = (dynhash_t)malloc(sizeof(*dh));
	if (dh == NULL)
		return (NULL);

	memset(dh, 0, sizeof(*dh));

	dh->table = (struct dynhash_entry **)malloc(DYNHASH_INIT_LENGTH *
	    sizeof(struct dynhash_entry *));
	if (dh->table == NULL) {
		free(dh);
		return (NULL);
	}

	memset(dh->table, 0, DYNHASH_INIT_LENGTH *
	    sizeof(struct dynhash_entry *));
	
	dh->table_length = DYNHASH_INIT_LENGTH;
	dh->table_count = 0;
	dh->keycmp = keycmp;
	dh->hashfn = hashfn;
	dh->iterating = false;

	return (dh);
}

/*
 * We permit this to be called on a non-empty structure, so the caller had
 * better be sure that they won't leave any dangling pointers.
 */
void
dynhash_destroy(dynhash_t dh)
{
	struct dynhash_entry *dhe;
	uint32_t i;

	assert(dh != NULL);
	assert(dh->table != NULL);

	ASSERT_NOT_ITERATING(dh);

	for (i = 0; i < dh->table_length; i++) {
		dhe = dh->table[i];
		while (dhe != NULL) {
			struct dynhash_entry *dhe_tmp;

			dhe_tmp = dhe;
			dhe = dhe->next;
			free(dhe_tmp);

			assert(dh->table_count != 0);
			dh->table_count--;
		}
	}

	free(dh->table);
	free(dh);
}

/*
 * Return the key/value associated with 'key', if it exists.
 */
bool
dynhash_find(dynhash_t dh, void *key, void **rkey, void **rvalue)
{
	struct dynhash_entry *dhe;

	assert(key != NULL);
	assert(dh != NULL);
	assert(dh->table != NULL);

	dhe = dh->table[dh->hashfn(key) % dh->table_length];
	while (dhe != NULL) {
		if (dh->keycmp(dhe->key, key) == 0) {
			if (rkey != NULL)
				*rkey = dhe->key;
			if (rvalue != NULL)
				*rvalue = dhe->val;
			break;
		}
		dhe = dhe->next;
	}

	return (dhe != NULL);
}

/*
 * Add a key/value pair to the dynhash structure.
 */
bool
dynhash_add(dynhash_t dh, void *key, void *value)
{
	struct dynhash_entry *dhe;
	uint32_t idx;

	assert(key != NULL);
	assert(dh != NULL);
	assert(dh->table != NULL);

	ASSERT_NOT_ITERATING(dh);

	if (dynhash_find(dh, key, NULL, NULL))
		return (false);

	if (dh->table_count == dh->table_length)
		dynhash_expand(dh);

	dhe = (struct dynhash_entry *)malloc(sizeof(*dhe));
	if (dhe == NULL)
		return (false);

	dhe->key = key;
	dhe->val = value;

	idx = dh->hashfn(key) % dh->table_length;

	dhe->next = dh->table[idx];
	dh->table[idx] = dhe;

	dh->table_count++;
	assert(dh->table_count != 0);

	return (true);
}

/*
 * Remove a key from the dynhash structure, returning the original key and
 * value, which may be pointers to comparable, but physically distinct
 * keys/values.
 */
bool
dynhash_remove(dynhash_t dh, void *key, void **rkey, void **rvalue)
{
	struct dynhash_entry *dhe, *dhe_prev;
	uint32_t idx;

	assert(key != NULL);
	assert(dh != NULL);
	assert(dh->table != NULL);

	ASSERT_NOT_ITERATING(dh);

	if (dynhash_find(dh, key, rkey, rvalue) == false)
		return (false);

	assert(dh->table_count > 0);

	idx = dh->hashfn(key) % dh->table_length;

	dhe_prev = NULL;
	dhe = dh->table[idx];
	while (dhe != NULL) {
		if (dh->keycmp(dhe->key, key) == 0) {
			if (rkey != NULL)
				*rkey = dhe->key;
			if (rvalue != NULL)
				*rvalue = dhe->val;
			if (dhe_prev != NULL)
				dhe_prev->next = dhe->next;
			else
				dh->table[idx] = dhe->next;
			free(dhe);
			break;
		}
		dhe_prev = dhe;
		dhe = dhe->next;
	}

	if (dhe == NULL)
		return (false);

	assert(dh->table_count != 0);
	dh->table_count--;

	return (true);
}

/*
 * Return the number of elements in our structure.
 */
uint32_t
dynhash_count(dynhash_t dh)
{
	
	assert(dh != NULL);

	return (dh->table_count);
}

/*
 * Process all key/value pairs in our structure with the provided
 * callback function.
 *
 * Do _NOT_ pass this function any callbacks which will augment 'dh',
 * otherwise very bad things may happen. For example: do not dynhash_iterate
 * over all entries and dynhash_remove them to clean up. That's asking for
 * trouble.
 */
void
dynhash_iterate(dynhash_t dh, void (*iter_func)(void *, void *, void *),
    void *arg)
{
	struct dynhash_entry *dhe;
	uint32_t i;

	assert(dh != NULL);
	assert(dh->table != NULL);
	assert(iter_func != NULL);

	dh->iterating = true;

	for (i = 0; i < dh->table_length; i++) {
		for (dhe = dh->table[i]; dhe != NULL; dhe = dhe->next)
			iter_func(arg, dhe->key, dhe->val);
	}

	dh->iterating = false;
}
