/*
** autogenerated content - DO NOT EDIT
*/
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

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include <time.h>
#include <sys/time.h>

#include "core/bitpackstring.h"
#include "core/error.h"
#include "core/ensure.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/yarandom.h"

enum {
/*   MAX_RND_NUMS = 10, */
  MAX_RND_NUMS = 100000,
};

static inline int
icmp(uint32_t a, uint32_t b)
{
  if (a > b)
    return 1;
  else if (a < b)
    return -1;
  else /* if (a == b) */
    return 0;
}

/**
 * \brief bit count reference
 * @param v count the number of bits set in v
 */
static inline int
genBitCount(uint32_t v)
{
  unsigned c; /* c accumulates the total bits set in v */
  for (c = 0; v; c++)
    v &= v - 1; /* clear the least significant bit set */
  return c;
}

#define freeResourcesAndReturn(retval) \
  do {                                 \
    ma_free(numBitsList);              \
    ma_free(randSrc);                  \
    ma_free(randCmp);                  \
    ma_free(bitStore);                 \
    ma_free(bitStoreCopy);             \
    return retval;                     \
  } while (0)

int
bitPackStringInt32_unit_test(GT_Error *err)
{
  BitString bitStore = NULL;
  BitString bitStoreCopy = NULL;
  uint32_t *randSrc = NULL; /*< create random ints here for input as bit
                                *  store */
  uint32_t *randCmp = NULL; /*< used for random ints read back */
  unsigned *numBitsList = NULL;
  size_t i, numRnd;
  BitOffset offsetStart, offset;
  int had_err = 0;
  offset = offsetStart = random()%(sizeof (uint32_t) * CHAR_BIT);
  numRnd = random() % (MAX_RND_NUMS + 1);
  log_log("offset=%lu, numRnd=%lu\n",
          (long unsigned)offsetStart, (long unsigned)numRnd);
  {
    BitOffset numBits = sizeof (uint32_t) * CHAR_BIT * numRnd + offsetStart;
    randSrc = ma_malloc(sizeof (uint32_t)*numRnd);
    bitStore = ma_malloc(bitElemsAllocSize(numBits) * sizeof (BitElem));
    bitStoreCopy = ma_calloc(bitElemsAllocSize(numBits), sizeof (BitElem));
    randCmp = ma_malloc(sizeof (uint32_t)*numRnd);
  }
  /* first test unsigned types */
  log_log("bsStoreUInt32/bsGetUInt32: ");
  for (i = 0; i < numRnd; ++i)
  {
#if 32 > 32 && LONG_BIT < 32
    uint32_t v = randSrc[i] = (uint32_t)random() << 32 | random();
#else /* 32 > 32 && LONG_BIT < 32 */
    uint32_t v = randSrc[i] = random();
#endif /* 32 > 32 && LONG_BIT < 32 */
    int bits = requiredUInt32Bits(v);
    bsStoreUInt32(bitStore, offset, bits, v);
    offset += bits;
  }
  offset = offsetStart;
  for (i = 0; i < numRnd; ++i)
  {
    uint32_t v = randSrc[i];
    int bits = requiredUInt32Bits(v);
    uint32_t r = bsGetUInt32(bitStore, offset, bits);
    ensure(had_err, r == v);
    if (had_err)
    {
      log_log("Expected %"PRIu32", got %"PRIu32", i = %lu\n",
              v, r, (unsigned long)i);
      freeResourcesAndReturn(had_err);
    }
    offset += bits;
  }
  log_log("passed\n");
  if (numRnd > 0)
  {
    uint32_t v = randSrc[0], r = 0;
    unsigned numBits = requiredUInt32Bits(v);
    BitOffset i = offsetStart + numBits;
    uint32_t mask = ~(uint32_t)0;
    if (numBits < 32)
      mask = ~(mask << numBits);
    log_log("bsSetBit, bsClearBit, bsToggleBit, bsGetBit: ");
    while (v)
    {
      int lowBit = v & 1;
      v >>= 1;
      ensure(had_err, lowBit == (r = bsGetBit(bitStore, --i)));
      if (had_err)
      {
        log_log("Expected %d, got %d, i = %llu\n",
                lowBit, (int)r, (unsigned long long)i);
        freeResourcesAndReturn(had_err);
      }
    }
    i = offsetStart + numBits;
    bsClear(bitStoreCopy, offsetStart, numBits, random()&1);
    v = randSrc[0];
    while (i)
    {
      int lowBit = v & 1;
      v >>= 1;
      if (lowBit)
        bsSetBit(bitStoreCopy, --i);
      else
        bsClearBit(bitStoreCopy, --i);
    }
    v = randSrc[0];
    r = bsGetUInt32(bitStoreCopy, offsetStart, numBits);
    ensure(had_err, r == v);
    if (had_err)
    {
      log_log("Expected %"PRIu32", got %"PRIu32"\n", v, r);
      freeResourcesAndReturn(had_err);
    }
    for (i = 0; i < numBits; ++i)
      bsToggleBit(bitStoreCopy, offsetStart + i);
    r = bsGetUInt32(bitStoreCopy, offsetStart, numBits);
    ensure(had_err, r == (v = (~v & mask)));
    if (had_err)
    {
      log_log("Expected %"PRIu32", got %"PRIu32"\n", v, r);
      freeResourcesAndReturn(had_err);
    }
    log_log("passed\n");
  }
  if (numRnd > 1)
  {
    log_log("bsCompare: ");
    {
      uint32_t v0 = randSrc[0];
      int bits0 = requiredUInt32Bits(v0);
      uint32_t r0;
      offset = offsetStart;
      r0 = bsGetUInt32(bitStore, offset, bits0);
      for (i = 1; i < numRnd; ++i)
      {
        uint32_t v1 = randSrc[i];
        int bits1 = requiredUInt32Bits(v1);
        uint32_t r1 = bsGetUInt32(bitStore, offset + bits0, bits1);
        int result = -2;   /*< -2 is not a return value of bsCompare, thus
                            *   if it is displayed, there was an earlier
                            *   error. */
        ensure(had_err, r0 == v0 && r1 == v1);
        ensure(had_err, icmp(v0, v1) ==
               (result = bsCompare(bitStore, offset, bits0,
                                   bitStore, offset + bits0, bits1)));
        if (had_err)
        {
          log_log("Expected v0 %s v1, got v0 %s v1,\n for v0=%"
                  PRIu32" and v1=%"PRIu32",\n"
                  "i = %lu, bits0=%u, bits1=%u\n",
                  (v0 > v1?">":(v0 < v1?"<":"==")),
                  (result > 0?">":(result < 0?"<":"==")), v0, v1,
                  (unsigned long)i, bits0, bits1);
          freeResourcesAndReturn(had_err);
        }
        offset += bits0;
        bits0 = bits1;
        v0 = v1;
        r0 = r1;
      }
    }
    log_log("passed\n");
  }
  log_log("bsStoreUniformUInt32Array/bsGetUInt32: ");
  {
    unsigned numBits = random()%32 + 1;
    uint32_t mask = ~(uint32_t)0;
    if (numBits < 32)
      mask = ~(mask << numBits);
    offset = offsetStart;
    bsStoreUniformUInt32Array(bitStore, offset, numBits, numRnd, randSrc);
    for (i = 0; i < numRnd; ++i)
    {
      uint32_t v = randSrc[i] & mask;
      uint32_t r = bsGetUInt32(bitStore, offset, numBits);
      ensure(had_err, r == v);
      if (had_err)
      {
        log_log("Expected %"PRIu32", got %"PRIu32",\n"
                "i = %lu, bits=%u\n", v, r, (unsigned long)i, numBits);
        freeResourcesAndReturn(had_err);
      }
      offset += numBits;
    }
    log_log("passed\n");
    log_log("bsStoreUniformUInt32Array/bsGetUniformUInt32Array: ");
    bsGetUniformUInt32Array(bitStore, offset = offsetStart,
                               numBits, numRnd, randCmp);
    for (i = 0; i < numRnd; ++i)
    {
      uint32_t v = randSrc[i] & mask;
      uint32_t r = randCmp[i];
      ensure(had_err, r == v);
      if (had_err)
      {
        log_log( "Expected %"PRIu32", got %"PRIu32",\n"
                " i = %lu, bits=%u\n",
                v, r, (unsigned long)i, numBits);
        freeResourcesAndReturn(had_err);
      }
    }
    if (numRnd > 1)
    {
      uint32_t v = randSrc[0] & mask;
      uint32_t r;
      bsGetUniformUInt32Array(bitStore, offsetStart,
                                 numBits, 1, &r);
      if (r != v)
      {
        log_log("Expected %"PRIu32", got %"PRIu32","
                " one value extraction\n",
                v, r);
        freeResourcesAndReturn(had_err);
      }
    }
    log_log(" passed\n");
  }
  /* int types */
  log_log("bsStoreInt32/bsGetInt32: ");
  for (i = 0; i < numRnd; ++i)
  {
    int32_t v = (int32_t)randSrc[i];
    unsigned bits = requiredInt32Bits(v);
    bsStoreInt32(bitStore, offset, bits, v);
    offset += bits;
  }
  offset = offsetStart;
  for (i = 0; i < numRnd; ++i)
  {
    int32_t v = randSrc[i];
    unsigned bits = requiredInt32Bits(v);
    int32_t r = bsGetInt32(bitStore, offset, bits);
    ensure(had_err, r == v);
    if (had_err)
    {
      log_log("Expected %"PRId32", got %"PRId32",\n"
                  "i = %lu, bits=%u\n",
                  v, r, (unsigned long)i, bits);
      freeResourcesAndReturn(had_err);
    }
    offset += bits;
  }
  log_log("passed\n");
  log_log("bsStoreUniformInt32Array/bsGetInt32: ");
  {
    unsigned numBits = random()%32 + 1;
    int32_t mask = ~(int32_t)0;
    if (numBits < 32)
      mask = ~(mask << numBits);
    offset = offsetStart;
    bsStoreUniformInt32Array(bitStore, offset, numBits, numRnd,
                                (int32_t *)randSrc);
    for (i = 0; i < numRnd; ++i)
    {
      int32_t m = (int32_t)1 << (numBits - 1);
      int32_t v = (int32_t)((randSrc[i] & mask) ^ m) - m;
      int32_t r = bsGetInt32(bitStore, offset, numBits);
      ensure(had_err, r == v);
      if (had_err)
      {
        log_log("Expected %"PRId32", got %"PRId32",\n"
                    "i = %lu, numBits=%u\n",
                    v, r, (unsigned long)i, numBits);
        freeResourcesAndReturn(had_err);
      }
      offset += numBits;
    }
    log_log("passed\n");
    log_log("bsStoreUniformInt32Array/bsGetUniformInt32Array: ");
    bsGetUniformInt32Array(bitStore, offset = offsetStart,
                              numBits, numRnd, (int32_t *)randCmp);
    for (i = 0; i < numRnd; ++i)
    {
      int32_t m = (int32_t)1 << (numBits - 1);
      int32_t v = (int32_t)((randSrc[i] & mask) ^ m) - m;
      int32_t r = randCmp[i];
      ensure(had_err, r == v);
      if (had_err)
      {
        log_log("Expected %"PRId32", got %"PRId32", i = %lu\n",
                v, r, (unsigned long)i);
        freeResourcesAndReturn(had_err);
      }
    }
    if (numRnd > 0)
    {
      int32_t m = (int32_t)1 << (numBits - 1);
      int32_t v = (int32_t)((randSrc[0] & mask) ^ m) - m;
      int32_t r = 0;
      bsGetUniformInt32Array(bitStore, offsetStart,
                                numBits, 1, &r);
      ensure(had_err, r == v);
      if (had_err)
      {
        log_log("Expected %"PRId32", got %"PRId32
                ", one value extraction\n",
                v, r);
        freeResourcesAndReturn(had_err);
      }
    }
    log_log("passed\n");
  }

  log_log("bsStoreNonUniformUInt32Array/bsGetUInt32: ");
  {
    BitOffset bitsTotal = 0;
    numBitsList = ma_malloc(sizeof(unsigned) * numRnd);
    for (i = 0; i < numRnd; ++i)
      bitsTotal += (numBitsList[i] = random()%32 + 1);
    offset = offsetStart;
    bsStoreNonUniformUInt32Array(bitStore, offset, numRnd, bitsTotal,
                                     numBitsList, randSrc);
    for (i = 0; i < numRnd; ++i)
    {
      unsigned numBits = numBitsList[i];
      uint32_t mask = (numBits < 32)?
        ~((~(uint32_t)0) << numBits):~(uint32_t)0;
      uint32_t v = randSrc[i] & mask;
      uint32_t r = bsGetUInt32(bitStore, offset, numBits);
      ensure(had_err, r == v);
      if (had_err)
      {
        log_log("Expected %"PRIu32", got %"PRIu32",\n"
                "i = %lu, bits=%u\n",
                v, r, (unsigned long)i, numBits);
        freeResourcesAndReturn(had_err);
      }
      offset += numBits;
    }
    log_log("passed\n");
    log_log("bsStoreNonUniformUInt32Array/"
            "bsGetNonUniformUInt32Array: ");
    bsGetNonUniformUInt32Array(bitStore, offset = offsetStart,
                                   numRnd, bitsTotal, numBitsList, randCmp);
    for (i = 0; i < numRnd; ++i)
    {
      unsigned numBits = numBitsList[i];
      uint32_t mask = (numBits < 32)?
        ~((~(uint32_t)0) << numBits):~(uint32_t)0;
      uint32_t v = randSrc[i] & mask,
        r = randCmp[i];
      ensure(had_err, r == v);
      if (had_err)
      {
        log_log( "Expected %"PRIu32", got %"PRIu32",\n"
                " i = %lu, bits=%u\n",
                v, r, (unsigned long)i, numBits);
        freeResourcesAndReturn(had_err);
      }
    }
    if (numRnd > 1)
    {
      unsigned numBits = numBitsList[0];
      uint32_t mask = (numBits < 32)?
        ~((~(uint32_t)0) << numBits):~(uint32_t)0;
      uint32_t v = randSrc[0] & mask;
      uint32_t r;
      bsGetNonUniformUInt32Array(bitStore, offsetStart, 1, numBits,
                                     numBitsList, &r);
      if (r != v)
      {
        log_log("Expected %"PRIu32", got %"PRIu32", "
                " one value extraction\n",
                v, r);
        freeResourcesAndReturn(had_err);
      }
    }
    log_log(" passed\n");
    ma_free(numBitsList);
    numBitsList = NULL;
  }
  log_log("bsNonStoreUniformInt32Array/bsGetInt32: ");
  {
    BitOffset bitsTotal = 0;
    numBitsList = ma_malloc(sizeof(unsigned) * numRnd);
    for (i = 0; i < numRnd; ++i)
      bitsTotal += (numBitsList[i] = random()%32 + 1);
    offset = offsetStart;
    bsStoreNonUniformInt32Array(bitStore, offset, numRnd, bitsTotal,
                                     numBitsList, (int32_t *)randSrc);
    for (i = 0; i < numRnd; ++i)
    {
      unsigned numBits = numBitsList[i];
      int32_t mask = (numBits < 32)
        ? ~((~(int32_t)0) << numBits) : ~(int32_t)0;
      int32_t m = (int32_t)1 << (numBits - 1);
      int32_t v = (int32_t)((randSrc[i] & mask) ^ m) - m;
      int32_t r = bsGetInt32(bitStore, offset, numBits);
      ensure(had_err, r == v);
      if (had_err)
      {
        log_log("Expected %"PRId32", got %"PRId32",\n"
                    "i = %lu, numBits=%u\n",
                    v, r, (unsigned long)i, numBits);
        freeResourcesAndReturn(had_err);
      }
      offset += numBits;
    }
    log_log("passed\n");
    log_log("bsStoreNonUniformInt32Array/"
            "bsGetNonUniformInt32Array: ");
    bsGetNonUniformInt32Array(bitStore, offset = offsetStart, numRnd,
                                   bitsTotal, numBitsList,
                                   (int32_t *)randCmp);
    for (i = 0; i < numRnd; ++i)
    {
      unsigned numBits = numBitsList[i];
      int32_t mask = (numBits < 32)
        ? ~((~(int32_t)0) << numBits) : ~(int32_t)0;
      int32_t m = (int32_t)1 << (numBits - 1);
      int32_t v = (int32_t)((randSrc[i] & mask) ^ m) - m;
      int32_t r = randCmp[i];
      ensure(had_err, r == v);
      if (had_err)
      {
        log_log("Expected %"PRId32", got %"PRId32", i = %lu\n",
                v, r, (unsigned long)i);
        freeResourcesAndReturn(had_err);
      }
    }
    if (numRnd > 0)
    {
      unsigned numBits = numBitsList[0];
      int32_t mask = (numBits < 32)
        ? ~((~(int32_t)0) << numBits) : ~(int32_t)0;
      int32_t m = (int32_t)1 << (numBits - 1);
      int32_t v = (int32_t)((randSrc[0] & mask) ^ m) - m;
      int32_t r = 0;
      bsGetNonUniformInt32Array(bitStore, offsetStart,
                                     1, numBits, numBitsList, &r);
      ensure(had_err, r == v);
      if (had_err)
      {
        log_log("Expected %"PRId32", got %"PRId32
                ", one value extraction\n",
                v, r);
        freeResourcesAndReturn(had_err);
      }
    }
    log_log("passed\n");
    ma_free(numBitsList);
    numBitsList = NULL;
  }

  if (numRnd > 0)
  {
    log_log("bsCopy: ");
    {
      /* first decide how many of the values to use and at which to start */
      size_t numValueCopies, copyStart;
      BitOffset numCopyBits = 0, destOffset;
      unsigned numBits = random()%32 + 1;
      uint32_t mask = ~(uint32_t)0;
      if (numBits < 32)
        mask = ~(mask << numBits);
      if (random()&1)
      {
        numValueCopies = random()%(numRnd + 1);
        copyStart = random()%(numRnd - numValueCopies + 1);
      }
      else
      {
        copyStart = random() % numRnd;
        numValueCopies = random()%(numRnd - copyStart) + 1;
      }
      assert(copyStart + numValueCopies <= numRnd);
      offset = offsetStart + (BitOffset)copyStart * numBits;
      bsStoreUniformUInt32Array(bitStore, offset, numBits, numValueCopies,
                                    randSrc);
      destOffset = random()%(offsetStart + 32
                             * (BitOffset)(numRnd - numValueCopies) + 1);
      numCopyBits = (BitOffset)numBits * numValueCopies;
      /* the following bsCopy should be equivalent to:
       * bsStoreUniformUInt32Array(bitStoreCopy, destOffset,
       *                              numBits, numValueCopies, randSrc); */
      bsCopy(bitStore, offset, bitStoreCopy, destOffset, numCopyBits);
      ensure(had_err,
             bsCompare(bitStore, offset, numCopyBits,
                       bitStoreCopy, destOffset, numCopyBits) == 0);
      if (had_err)
      {
        log_log("Expected equality on bitstrings\n"
                    "offset = %llu, destOffset = %llu,"
                    " numCopyBits=%llu\n",
                    (unsigned long long)offset,
                    (unsigned long long)destOffset,
                    (unsigned long long)numCopyBits);
        /* FIXME: implement bitstring output function */
        freeResourcesAndReturn(had_err);
      }
      log_log("passed\n");
    }
  }
  if (numRnd > 0)
  {
    log_log("bsClear: ");
    {
      /* first decide how many of the values to use and at which to start */
      size_t numResetValues, resetStart;
      BitOffset numResetBits = 0;
      unsigned numBits = random()%32 + 1;
      int bitVal = random()&1;
      int32_t cmpVal = bitVal?-1:0;
      uint32_t mask = ~(uint32_t)0;
      if (numBits < 32)
        mask = ~(mask << numBits);
      if (random()&1)
      {
        numResetValues = random()%(numRnd + 1);
        resetStart = random()%(numRnd - numResetValues + 1);
      }
      else
      {
        resetStart = random() % numRnd;
        numResetValues = random()%(numRnd - resetStart) + 1;
      }
      assert(resetStart + numResetValues <= numRnd);
      offset = offsetStart;
      bsStoreUniformInt32Array(bitStore, offset, numBits, numRnd,
                                    (int32_t *)randSrc);
      numResetBits = (BitOffset)numBits * numResetValues;
      bsClear(bitStore, offset + (BitOffset)resetStart * numBits,
              numResetBits, bitVal);
      {
        int32_t m = (int32_t)1 << (numBits - 1);
        for (i = 0; i < resetStart; ++i)
        {
          int32_t v = (int32_t)((randSrc[i] & mask) ^ m) - m;
          int32_t r = bsGetInt32(bitStore, offset, numBits);
          ensure(had_err, r == v);
          if (had_err)
          {
            log_log( "Expected %"PRId32", got %"PRId32",\n"
                     "i = %lu, numBits=%u\n",
                     v, r, (unsigned long)i, numBits);
            freeResourcesAndReturn(had_err);
          }
          offset += numBits;
        }
        for (; i < resetStart + numResetValues; ++i)
        {
          int32_t r = bsGetInt32(bitStore, offset, numBits);
          ensure(had_err, r == cmpVal);
          if (had_err)
          {
            log_log("Expected %"PRId32", got %"PRId32",\n"
                    "i = %lu, numBits=%u\n",
                    cmpVal, r, (unsigned long)i, numBits);
            freeResourcesAndReturn(had_err);
          }
          offset += numBits;
        }
        for (; i < numRnd; ++i)
        {
          int32_t v = (int32_t)((randSrc[i] & mask) ^ m) - m;
          int32_t r = bsGetInt32(bitStore, offset, numBits);
          ensure(had_err, r == v);
          if (had_err)
          {
            log_log("Expected %"PRId32", got %"PRId32",\n"
                    "i = %lu, numBits=%u\n",
                    v, r, (unsigned long)i, numBits);
            freeResourcesAndReturn(had_err);
          }
          offset += numBits;
        }
      }
    }
    log_log("passed\n");
  }
  if (numRnd > 0)
  {
    log_log("bs1BitsCount: ");
    {
      /* first decide how many of the values to use and at which to start */
      size_t numCountValues, countStart;
      BitOffset numCountBits = 0, bitCountRef = 0, bitCountCmp;
      unsigned numBits = random()%32 + 1;
      uint32_t mask = ~(uint32_t)0;
      if (numBits < 32)
        mask = ~(mask << numBits);
      if (random()&1)
      {
        numCountValues = random()%(numRnd + 1);
        countStart = random()%(numRnd - numCountValues + 1);
      }
      else
      {
        countStart = random() % numRnd;
        numCountValues = random()%(numRnd - countStart) + 1;
      }
      assert(countStart + numCountValues <= numRnd);
      offset = offsetStart;
      bsStoreUniformUInt32Array(bitStore, offset, numBits, numRnd, randSrc);
      numCountBits = (BitOffset)numBits * numCountValues;
      bitCountCmp = bs1BitsCount(bitStore,
                                 offset + (BitOffset)countStart * numBits,
                                 numCountBits);
      for (i = countStart; i < countStart + numCountValues; ++i)
      {
        uint32_t v = (uint32_t)randSrc[i] & mask;
        bitCountRef += genBitCount(v);
      }
      ensure(had_err, bitCountRef == bitCountCmp);
      if (had_err)
      {
        log_log("Expected %llu, got %llu,\n"
                "numBits=%u\n", (unsigned long long)bitCountRef,
                (unsigned long long)bitCountCmp, numBits);
        freeResourcesAndReturn(had_err);
      }
      offset += numBits;
    }
    log_log("passed\n");
  }
  freeResourcesAndReturn(had_err);
}
