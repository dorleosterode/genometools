/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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
#include "core/symboldef.h"
#include "core/xansi.h"
#include "core/arraydef.h"
#include "seqpos-def.h"
#include "sfx-lcpsub.h"
#include "lcpoverflow.h"

void outlcpvalues(Lcpsubtab *lcpsubtab,
                  unsigned long bucketleft,
                  unsigned long bucketright,
                  Seqpos posoffset,
                  FILE *fplcptab,
                  FILE *fpllvtab)
{
  unsigned long idx;
  Seqpos lcpvalue;
  Largelcpvalue *largelcpvalueptr;

  lcpsubtab->largelcpvalues.nextfreeLargelcpvalue = 0;
  if (lcpsubtab->numoflargelcpvalues >=
      (Seqpos) lcpsubtab->largelcpvalues.allocatedLargelcpvalue)
  {
    lcpsubtab->largelcpvalues.spaceLargelcpvalue
      = gt_realloc(lcpsubtab->largelcpvalues.spaceLargelcpvalue,
                   sizeof (Largelcpvalue) * lcpsubtab->numoflargelcpvalues);
    lcpsubtab->largelcpvalues.allocatedLargelcpvalue
      = (unsigned long) lcpsubtab->numoflargelcpvalues;
  }
  for (idx=bucketleft; idx<=bucketright; idx++)
  {
    lcpvalue = lcpsubtab->spaceSeqpos[idx];
    if (lcpsubtab->maxbranchdepth < lcpvalue)
    {
      lcpsubtab->maxbranchdepth = lcpvalue;
    }
    if (lcpvalue < (Seqpos) LCPOVERFLOW)
    {
      lcpsubtab->smalllcpvalues[idx] = (Uchar) lcpvalue;
    } else
    {
      if (lcpsubtab->largelcpvalues.nextfreeLargelcpvalue >=
          lcpsubtab->largelcpvalues.allocatedLargelcpvalue)
      {
        fprintf(stderr,"at pos %lu: lcpvalue= %lu, "
                       "nextfreeLargelcpvalue = %lu >= %lu\n",
                 (unsigned long) (posoffset+idx),
                 (unsigned long) lcpvalue,
                 lcpsubtab->largelcpvalues.nextfreeLargelcpvalue,
                 lcpsubtab->largelcpvalues.allocatedLargelcpvalue);
        exit(EXIT_FAILURE);
      }
      gt_assert(lcpsubtab->largelcpvalues.nextfreeLargelcpvalue <
                lcpsubtab->largelcpvalues.allocatedLargelcpvalue);
      largelcpvalueptr = lcpsubtab->largelcpvalues.spaceLargelcpvalue +
                         lcpsubtab->largelcpvalues.nextfreeLargelcpvalue++;
      largelcpvalueptr->position = posoffset+idx;
      largelcpvalueptr->value = lcpvalue;
      lcpsubtab->smalllcpvalues[idx] = LCPOVERFLOW;
    }
  }
  lcpsubtab->countoutputlcpvalues += (bucketright - bucketleft + 1);
  gt_xfwrite(lcpsubtab->smalllcpvalues,
             sizeof (Uchar),(size_t) (bucketright - bucketleft + 1),fplcptab);
  if (lcpsubtab->largelcpvalues.nextfreeLargelcpvalue > 0)
  {
    printf("add %lu large lcps to obtain",
            (unsigned long) lcpsubtab->numoflargelcpvalues);
    lcpsubtab->totalnumoflargelcpvalues += lcpsubtab->numoflargelcpvalues;
    printf("%lu\n",
            (unsigned long) lcpsubtab->totalnumoflargelcpvalues);
    gt_xfwrite(lcpsubtab->largelcpvalues.spaceLargelcpvalue,
               sizeof (Largelcpvalue),
               (size_t)
               lcpsubtab->largelcpvalues.nextfreeLargelcpvalue,
               fpllvtab);
  }
}

#define NUMBEROFZEROS 1024

void outmany0lcpvalues(Lcpsubtab *lcpsubtab,Seqpos totallength,
                       FILE *outfplcptab)
{
  Seqpos i, countout, many;
  Uchar outvalues[NUMBEROFZEROS] = {0};

  many = totallength + 1 - lcpsubtab->countoutputlcpvalues;
  countout = many/NUMBEROFZEROS;
  for (i=0; i<countout; i++)
  {
    gt_xfwrite(outvalues,sizeof (Uchar),(size_t) NUMBEROFZEROS,outfplcptab);
  }
  gt_xfwrite(outvalues,sizeof (Uchar),(size_t) many % NUMBEROFZEROS,
             outfplcptab);
  lcpsubtab->countoutputlcpvalues += many;
}
