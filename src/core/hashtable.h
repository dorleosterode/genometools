/*
  Copyright (c) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <inttypes.h>
#include <stdlib.h>

#include "core/error.h"
#include "core/fptr.h"

typedef struct Hashtable Hashtable;

enum iterator_op
{
  CONTINUE_ITERATION,
  STOP_ITERATION,
  DELETED_ELEM,
  MODIFIED_KEY,
  REDO_ITERATION,
};

typedef enum iterator_op (*Elemvisitfunc)(void *elem, void *data,
                                             GT_Error *err);

typedef void (*FreeFuncWData)(void *elem, void *table_data);

typedef uint32_t htsize_t;
typedef htsize_t (*HashFunc)(const void *elem);

struct HashElemInfo
{
  HashFunc keyhash;
  union
  {
    FreeFunc free_elem;
    FreeFuncWData free_elem_with_data;
  } free_op;                           /**< either of these can be
                                        * used for the individual
                                        * destructors or set to NULL */
  size_t elem_size;
  Compare cmp;
  void *table_data;             /**< per table data, passed to
                                 * free_elem_with_data */
  FreeFunc table_data_free;     /**< called on hashtable_delete with
                                 * table_data as argument if != NULL */
};

typedef struct HashElemInfo HashElemInfo;

extern Hashtable *
hashtable_new(HashElemInfo);
extern Hashtable *
hashtable_new_with_start_size(HashElemInfo htype, unsigned short size_log);
void*      hashtable_get(Hashtable*, const void *elem);
/**
 * @return 1 if add succeeded, 0 if elem is already in table.
 */
int        hashtable_add(Hashtable*, const void *elem);
int        hashtable_remove(Hashtable*, const void *elem);
/**
 * @brief iterate over the hashtable in key order given by compare
 * function <cmp>
 * @return 0 => no error, -1 => error occured
 */
extern int
hashtable_foreach_ordered(Hashtable *ht, Elemvisitfunc iter, void *data,
                          Compare cmp, GT_Error *err);
/**
 * @brief iterate over the hashtable in implementation-defined order
 * @return 0 => no error, -1 => error occured
 */
extern int
hashtable_foreach(Hashtable *ht, Elemvisitfunc iter, void *data,
                  GT_Error *err);
/* iterate over the hashtable in default order. Requires that the hashtable
   was constructed with an ordering compare function. */
extern int
hashtable_foreach_in_default_order(Hashtable*, Elemvisitfunc, void *data,
                                   GT_Error *err);
size_t
hashtable_fill(Hashtable *);
void       hashtable_reset(Hashtable*);
int        hashtable_unit_test(GT_Error*);
void       hashtable_delete(Hashtable*);

/*
 * helper functions for users constructing their own HashElemInfos
 */
static inline uint32_t
uint32_key_mul_hash(uint32_t key);

/**
 * @brief Hash pointer by address value
 * @param elem treated as a void ** so that *(void **)elem is used as
 * key to hash
 */
extern uint32_t
ht_ptr_elem_hash(const void *elem);

/**
 * @brief Hash unsigned long by value
 * @param elem *(unsigned long *)elem is used as key to hash
 */
extern uint32_t
ht_ul_elem_hash(const void *elem);

/**
 * @brief Hash string by mixing hash value of characters
 *
 * @param elem *(char **)elem is used as key to hash
 */
extern uint32_t
ht_cstr_elem_hash(const void *elem);
/**
 * @brief hash binary data
 */
extern uint32_t
uint32_data_hash(const void *data, size_t length);
extern int
ht_ptr_elem_cmp(const void *elemA, const void *elemB);
static inline int
ht_ul_cmp(unsigned long a, unsigned long b);
extern int
ht_ul_elem_cmp(const void *elemA, const void *elemB);
extern int
ht_cstr_elem_cmp(const void *elemA, const void *elemB);

#include "core/hashtable-siop.h"

#endif
