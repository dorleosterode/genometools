/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/intbits.h"
#include "core/types_api.h"
#include "core/codetype.h"
#include "core/divmodmul.h"
#include "core/encseq.h"
#include "kmercodes.h"
#include "firstcodes-buf.h"
#include "firstcodes-insert.h"

#define GT_FIRSTCODES_INSERTSUFFIXES(BUF,CODE,SEQNUM,RELPOS)\
        {\
          if ((BUF)->currentmincode <= (CODE) &&\
              (CODE) <= (BUF)->currentmaxcode &&\
              GT_MARKSUBSTRING_CHECKMARK((BUF)->markprefix,CODE) &&\
              GT_MARKSUBSTRING_CHECKMARK((BUF)->marksuffix,CODE))\
          {\
            if ((BUF)->nextfree == (BUF)->allocated)\
            {\
              (BUF)->flush_function((BUF)->fciptr);\
            }\
            gt_assert ((BUF)->nextfree < (BUF)->allocated);\
            (BUF)->spaceGtUlongPair[(BUF)->nextfree].a = CODE;\
            (BUF)->spaceGtUlongPair[(BUF)->nextfree++].b\
              = gt_seqnumrelpos_encode((BUF)->snrp,SEQNUM,RELPOS);\
          }\
        }

static void gt_firstcodes_insert_kmerscan_range(
                                         const GtBitsequence *twobitencoding,
                                         unsigned int kmersize,
                                         unsigned int minmatchlength,
                                         unsigned long startpos,
                                         unsigned long endpos,
                                         unsigned long fseqnum,
                                         unsigned long rseqnum,
                                         unsigned long maxunitindex,
                                         GtCodeposbuffer *buf)
{
  unsigned long position, unitindex, frelpos;
  unsigned int shiftright;
  const unsigned int shiftleft = GT_MULT2(kmersize-1);
  const unsigned long maskright = GT_MASKRIGHT(kmersize);
  const unsigned long lastpossiblepos = endpos - startpos - minmatchlength;
  const unsigned long lastfrelpos = endpos - startpos - kmersize;
  GtTwobitencoding currentencoding;
  GtCodetype cc, marksubstringtmpcode, fcode, rccode;

  gt_assert(kmersize <= (unsigned int) GT_UNITSIN2BITENC);
  position = startpos;
  fcode = gt_kmercode_at_position(twobitencoding, position, kmersize);
  rccode = gt_kmercode_complement(gt_kmercode_reverse(fcode,kmersize),
                                  maskright);
  GT_FIRSTCODES_INSERTSUFFIXES(buf,fcode,fseqnum,0);
  if (lastfrelpos <= lastpossiblepos)
  {
    GT_FIRSTCODES_INSERTSUFFIXES(buf,rccode,rseqnum,lastfrelpos);
  }
  unitindex = GT_DIVBYUNITSIN2BITENC(startpos + kmersize);
  currentencoding = twobitencoding[unitindex];
  shiftright = (unsigned int)
               GT_MULT2(GT_UNITSIN2BITENC - 1 -
                        GT_MODBYUNITSIN2BITENC(startpos + kmersize));
  gt_assert(endpos >= (unsigned long) kmersize);
  endpos -= kmersize;
  frelpos = 1UL;
  while (position < endpos)
  {
    position++;
    cc = (GtCodetype) (currentencoding >> shiftright) & 3;
    fcode = ((fcode << 2) | cc) & maskright;
    rccode = (rccode >> 2) | ((cc ^ 3UL) << shiftleft);
    gt_assert(lastfrelpos >= frelpos);
    if (frelpos <= lastpossiblepos)
    {
      GT_FIRSTCODES_INSERTSUFFIXES(buf,fcode,fseqnum,frelpos);
    }
    if (lastfrelpos - frelpos <= lastpossiblepos)
    {
      GT_FIRSTCODES_INSERTSUFFIXES(buf,rccode,rseqnum,lastfrelpos-frelpos);
    }
    if (shiftright > 0)
    {
      shiftright -= 2;
    } else
    {
      gt_assert(unitindex < maxunitindex-1 || position == endpos);
      if (unitindex < maxunitindex-1)
      {
        currentencoding = twobitencoding[++unitindex];
        shiftright = (unsigned int) (GT_INTWORDSIZE-2);
      }
    }
    frelpos++;
  }
}

static void gt_firstcodes_insert_kmerscan_eqlen(
                                     const GtBitsequence *twobitencoding,
                                     unsigned long equallength,
                                     unsigned long totallength,
                                     unsigned long numofsequences,
                                     unsigned long maxunitindex,
                                     unsigned int kmersize,
                                     unsigned int minmatchlength,
                                     GtCodeposbuffer *buf)
{
  unsigned long startpos, fseqnum;

  if (equallength > (unsigned long) kmersize)
  {
    for (startpos = 0, fseqnum = 0; startpos < totallength;
         startpos += equallength+1, fseqnum++)
    {
      gt_firstcodes_insert_kmerscan_range(twobitencoding,
                                          kmersize,
                                          minmatchlength,
                                          startpos,
                                          startpos + equallength,
                                          fseqnum,
                                          numofsequences - 1 - fseqnum,
                                          maxunitindex,
                                          buf);
    }
  }
}

static void gt_firstcodes_insert_kmerscan(const GtEncseq *encseq,
                                          const GtBitsequence *twobitencoding,
                                          unsigned long totallength,
                                          unsigned long numofsequences,
                                          unsigned long maxunitindex,
                                          unsigned int kmersize,
                                          unsigned int minmatchlength,
                                          GtCodeposbuffer *buf)
{
  unsigned long laststart = 0, fseqnum = 0;

  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;

    sri = gt_specialrangeiterator_new(encseq,true);
    while (gt_specialrangeiterator_next(sri,&range)
           && range.start < totallength)
    {
      gt_assert(range.start >= laststart);
      if (range.start - laststart >= (unsigned long) kmersize)
      {
        gt_firstcodes_insert_kmerscan_range(twobitencoding,
                                            kmersize,
                                            minmatchlength,
                                            laststart,
                                            range.start,
                                            fseqnum,
                                            numofsequences - 1 - fseqnum,
                                            maxunitindex,
                                            buf);
      }
      laststart = range.end;
      fseqnum++;
    }
    gt_specialrangeiterator_delete(sri);
  }
  if (totallength - laststart >= (unsigned long) kmersize)
  {
    gt_firstcodes_insert_kmerscan_range(twobitencoding,
                                        kmersize,
                                        minmatchlength,
                                        laststart,
                                        totallength,
                                        fseqnum,
                                        numofsequences - 1 - fseqnum,
                                        maxunitindex,
                                        buf);
  }
}

void gt_firstcodes_insert_runkmerscan(const GtEncseq *encseq,
                                      unsigned int kmersize,
                                      unsigned int minmatchlength,
                                      GtCodeposbuffer *buf)
{
  const GtTwobitencoding *twobitencoding
    = gt_encseq_twobitencoding_export(encseq);
  unsigned long totallength, maxunitindex, numofsequences;

  if (gt_encseq_is_mirrored(encseq))
  {
    totallength = (gt_encseq_total_length(encseq)-1)/2;
  } else
  {
    totallength = gt_encseq_total_length(encseq);
  }
  numofsequences = gt_encseq_num_of_sequences(encseq);
  maxunitindex = gt_unitsoftwobitencoding(totallength) - 1;
  if (gt_encseq_accesstype_get(encseq) == GT_ACCESS_TYPE_EQUALLENGTH)
  {
    unsigned long equallength = gt_encseq_equallength(encseq);

    gt_firstcodes_insert_kmerscan_eqlen(twobitencoding,
                                        equallength,
                                        totallength,
                                        numofsequences,
                                        maxunitindex,
                                        kmersize,
                                        minmatchlength,
                                        buf);
  } else
  {
    gt_firstcodes_insert_kmerscan(encseq,
                                  twobitencoding,
                                  totallength,
                                  numofsequences,
                                  maxunitindex,
                                  kmersize,
                                  minmatchlength,
                                  buf);
  }
}