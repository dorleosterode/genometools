/*
  Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

#include "core/bitpackarray.h"
#include "core/ensure.h"
#include "core/log.h"
#include "core/yarandom.h"

enum {
/*   MAX_RND_NUMS = 10, */
  MAX_RND_NUMS = 100000,
};

int gt_bitpackarray_unit_test(GT_Error *err)
{
  struct BitPackArray *bitStore = NULL;
  int had_err = 0;
  {
    uint32_t *randSrc = NULL; /*< create random ints here for input as bit
                               *  store */
    uint32_t *randCmp = NULL, /*< used for random ints read back */
       mask;
    size_t i, numRnd;
    unsigned bits;
    numRnd = random() % MAX_RND_NUMS + 1;
    bits = random() % (sizeof (randSrc[0]) * CHAR_BIT + 1);
    if (bits == 32)
      mask = ~(uint32_t)0;
    else
      mask = ~((~(uint32_t)0)<<bits);

    log_log("numRnd=%lu\n", (long unsigned)numRnd);
    randSrc = ma_malloc(sizeof (uint32_t)*numRnd);
    bitStore = bitpackarray_new(bits, numRnd);
    randCmp = ma_malloc(sizeof (uint32_t)*numRnd);
    for (i = 0; i < numRnd; ++i)
    {
      uint32_t v = randSrc[i] = random();
      bitpackarray_store_uint32(bitStore, i, v);
    }
    for (i = 0; i < numRnd; ++i)
    {
      uint32_t v = randSrc[i];
      uint32_t r = bitpackarray_get_uint32(bitStore, i);
      ensure(had_err, (v & mask) == r);
      if (had_err)
      {
        log_log("bsStoreUInt32/bitpackarray_get_uint32: "
                "Expected %"PRIu32", got %"PRIu32", i = %lu, bits=%u\n",
                v & mask, r, (unsigned long)i, bits);
        ma_free(randSrc);
        ma_free(randCmp);
        bitpackarray_delete(bitStore);
        return had_err;
      }
    }
  ma_free(randSrc);
  ma_free(randCmp);
  bitpackarray_delete(bitStore);
  }
  log_log("bitpackarray_store_uint32/bitpackarray_get_uint32: passed\n");
  {
    uint64_t *randSrc = NULL; /*< create random ints here for input as bit
                        *  store */
    uint64_t *randCmp = NULL, /*< used for random ints read back */
       mask;
    size_t i, numRnd;
    unsigned bits;
    numRnd = random() % MAX_RND_NUMS + 1;
    bits = random() % (sizeof (randSrc[0]) * CHAR_BIT + 1);
    if (bits == (sizeof (randSrc[0]) * CHAR_BIT))
      mask = ~(uint64_t)0;
    else
      mask = ~((~(uint64_t)0)<<bits);
    ensure(had_err, (randSrc = ma_malloc(sizeof (uint64_t)*numRnd))
           && (bitStore = bitpackarray_new(bits, numRnd))
           && (randCmp = ma_malloc(sizeof (uint64_t)*numRnd)));
    if (had_err)
    {
      perror("Storage allocations failed");
      if (randSrc)
        ma_free(randSrc);
      if (randCmp)
        ma_free(randCmp);
      if (bitStore)
        bitpackarray_delete(bitStore);
      return had_err;
    }
    for (i = 0; i < numRnd; ++i)
    {
      uint64_t v = randSrc[i] = ((uint64_t)random() << 32 | random());
      bitpackarray_store_uint64(bitStore, i, v);
    }
    for (i = 0; i < numRnd; ++i)
    {
      uint64_t v = randSrc[i];
      uint64_t r = bitpackarray_get_uint64(bitStore, i);
      ensure(had_err, (v & mask) == r);
      if (had_err)
      {
        log_log("bsStoreUInt64/bitpackarray_get_uint64: "
                "Expected %llu, got %llu, i = %lu, bits=%u\n",
                (unsigned long long)(v & mask),
                (unsigned long long)r, (unsigned long)i, bits);
        ma_free(randSrc);
        ma_free(randCmp);
        bitpackarray_delete(bitStore);
        return had_err;
      }
    }
    ma_free(randSrc);
    ma_free(randCmp);
    bitpackarray_delete(bitStore);
  }
  log_log("bitpackarray_store_uint64/bitpackarray_get_uint64: passed\n");
  return had_err;
}
