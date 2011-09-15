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

#include "core/assert_api.h"
#include "core/ma_api.h"
#include "core/codetype.h"
#include "bcktab.h"
#include "sfx-partssuf.h"

typedef struct
{
  unsigned long nextidx,
                widthofpart,
                suftaboffset,
                sumofwidth;
} GtSuftabpartcomponent;

struct GtSuftabparts
{
  GtSuftabpartcomponent *components;
  unsigned int numofparts;
  unsigned long largestsizemappedpartwise,
                largestsuftabwidth;
};

#ifdef SKDEBUG
static void showrecord(const GtSuftabpartcomponent *component)
{
  printf("# part: width=%lu offset=%lu sumwidth=%lu nextidx=%lu\n",
          component->widthofpart,component->suftaboffset,
                                 component->sumofwidth,
                                 component->nextidx);
}

static void showallrecords(const GtSuftabparts *suftabparts)
{
  unsigned int idx;

  gt_assert(suftabparts != NULL);
  for (idx = 0; idx < suftabparts->numofparts; idx++)
  {
    showrecord(suftabparts->components + idx);
  }
}
#endif

static void gt_suftabparts_removeemptyparts(GtSuftabparts *suftabparts,
                                            GtLogger *logger)
{
#ifdef SKDEBUG
  printf("# before removial\n");
  showallrecords(suftabparts);
#endif
  gt_assert(suftabparts != NULL);
  if (suftabparts->numofparts > 0)
  {
    unsigned int destpart, srcpart;
    for (destpart = 0, srcpart = 0; srcpart < suftabparts->numofparts;
         srcpart++)
    {
      if (suftabparts->components[srcpart].widthofpart > 0)
      {
        if (destpart < srcpart)
        {
          suftabparts->components[destpart] = suftabparts->components[srcpart];
        }
        destpart++;
      }
    }
    if (suftabparts->components[suftabparts->numofparts-1].widthofpart == 0)
    {
      gt_assert(suftabparts->numofparts > 1U);
      destpart = suftabparts->numofparts-2;
      while (true)
      {
        if (suftabparts->components[destpart].widthofpart > 0)
        {
          suftabparts->components[destpart].nextidx
            = suftabparts->components[suftabparts->numofparts-1].nextidx;
          suftabparts->numofparts = destpart + 1;
          break;
        }
        if (destpart > 0)
        {
          destpart--;
        } else
        {
          gt_assert(false);
        }
      }
    } else
    {
      if (destpart < srcpart)
      {
        suftabparts->numofparts -= (srcpart - destpart);
        gt_assert(suftabparts->numofparts == destpart);
      }
    }
    for (srcpart = 0; srcpart < suftabparts->numofparts; srcpart++)
    {
      gt_assert(suftabparts->components[srcpart].widthofpart > 0);
      gt_logger_log(logger,"widthofpart[%u]=%lu",
                    srcpart,
                    suftabparts->components[srcpart].widthofpart);
    }
  }
#ifdef SKDEBUG
  printf("#after removal\n");
  showallrecords(suftabparts);
#endif
}

GtSuftabparts *gt_suftabparts_new(unsigned int numofparts,
                               const GtBcktab *bcktab,
                               const GtSfxmappedrange *mappedmarkprefixbuckets,
                               const GtSfxmappedrangelist *sfxmrlist,
                               unsigned long numofsuffixestoinsert,
                               unsigned long fullspecials,
                               GtLogger *logger)
{
  GtSuftabparts *suftabparts;
  unsigned long sizemapped;

  suftabparts = gt_malloc(sizeof *suftabparts);
  suftabparts->largestsizemappedpartwise = 0;
  gt_assert(suftabparts != NULL);
  if (numofsuffixestoinsert == 0)
  {
    suftabparts->numofparts = 0;
  } else
  {
    if (numofsuffixestoinsert < (unsigned long) numofparts ||
        gt_bcktab_prefixlength(bcktab) == 1U)
    {
      suftabparts->numofparts = 1U;
    } else
    {
      suftabparts->numofparts = numofparts;
    }
  }
  if (suftabparts->numofparts == 0)
  {
    suftabparts->largestsuftabwidth = fullspecials/numofparts+1;
    suftabparts->largestsizemappedpartwise = gt_bcktab_size_lb_cs(bcktab);
    if (mappedmarkprefixbuckets != NULL)
    {
      suftabparts->largestsizemappedpartwise
        += gt_Sfxmappedrange_size_entire(mappedmarkprefixbuckets);
    }
    if (suftabparts->largestsizemappedpartwise !=
        gt_Sfxmappedrangelist_size_entire(sfxmrlist))
    {
      fprintf(stderr,"largestsizemapped = %lu != %lu = size_entire\n",
              suftabparts->largestsizemappedpartwise,
              gt_Sfxmappedrangelist_size_entire(sfxmrlist));
      exit(EXIT_FAILURE);
    }
    suftabparts->components = NULL;
  } else
  {
    unsigned int part, remainder;
    unsigned long suftaboffset = 0, sumofwidth = 0;
    const unsigned long widthofsuftabpart
      = numofsuffixestoinsert/suftabparts->numofparts;

    suftabparts->components
      = gt_malloc(sizeof (*suftabparts->components) * suftabparts->numofparts);
    remainder = (unsigned int) (numofsuffixestoinsert %
                                (unsigned long) suftabparts->numofparts);
    suftabparts->largestsuftabwidth = 0;
    for (part=0; part < suftabparts->numofparts; part++)
    {
      if (remainder > 0)
      {
        suftaboffset += widthofsuftabpart + 1;
        remainder--;
      } else
      {
        suftaboffset += widthofsuftabpart;
      }
      if (part == suftabparts->numofparts - 1)
      {
        suftabparts->components[part].nextidx
          = gt_bcktab_numofallcodes(bcktab);
      } else
      {
        suftabparts->components[part].nextidx
          = gt_bcktab_findfirstlarger(bcktab,suftaboffset);
      }
      if (part == 0)
      {
        suftabparts->components[part].widthofpart
          = gt_bcktab_get_leftborder(bcktab,
                                     suftabparts->components[part].nextidx);
        suftabparts->components[part].suftaboffset = 0;
      } else
      {
        suftabparts->components[part].widthofpart
          = gt_bcktab_get_leftborder(bcktab,
                                     suftabparts->components[part].nextidx) -
            gt_bcktab_get_leftborder(bcktab,
                                     suftabparts->components[part-1].nextidx);
        suftabparts->components[part].suftaboffset
          = gt_bcktab_get_leftborder(bcktab,
                                     suftabparts->components[part-1].nextidx);
      }
      if (suftabparts->largestsuftabwidth <
          suftabparts->components[part].widthofpart)
      {
        suftabparts->largestsuftabwidth
          = suftabparts->components[part].widthofpart;
      }
      sumofwidth += suftabparts->components[part].widthofpart;
      suftabparts->components[part].sumofwidth = sumofwidth;
      sizemapped
        = gt_bcktab_mapped_range_size(bcktab,
                                      gt_suftabparts_mincode(part,suftabparts),
                                      gt_suftabparts_maxcode(part,suftabparts));
      if (mappedmarkprefixbuckets != NULL)
      {
        sizemapped += gt_Sfxmappedrange_size_mapped(mappedmarkprefixbuckets,
                                      gt_suftabparts_mincode(part,suftabparts),
                                      gt_suftabparts_maxcode(part,suftabparts));
      }
      if (suftabparts->largestsizemappedpartwise < sizemapped)
      {
        suftabparts->largestsizemappedpartwise = sizemapped;
      }
    }
    gt_assert(sumofwidth == numofsuffixestoinsert);
  }
  gt_suftabparts_removeemptyparts(suftabparts,logger);
  return suftabparts;
}

void gt_suftabparts_delete(GtSuftabparts *suftabparts)
{
  if (suftabparts != NULL)
  {
    gt_free(suftabparts->components);
    gt_free(suftabparts);
  }
}

GtCodetype gt_suftabparts_mincode(unsigned int part,
                                  const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && part < suftabparts->numofparts);
  /*  XXX: use nextidx as index */
  return (part == 0) ? 0 : suftabparts->components[part-1].nextidx + 1;
}

GtCodetype gt_suftabparts_maxcode(unsigned int part,
                                  const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && part < suftabparts->numofparts);
  return (part == suftabparts->numofparts - 1)
           ? suftabparts->components[part].nextidx - 1
           : suftabparts->components[part].nextidx;
}

unsigned long gt_suftabparts_offset(unsigned int part,
                                    const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && part < suftabparts->numofparts);
  return suftabparts->components[part].suftaboffset;
}

unsigned long gt_suftabparts_sumofwdith(unsigned int part,
                                        const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && part < suftabparts->numofparts);
  return suftabparts->components[part].sumofwidth;
}

unsigned long gt_suftabparts_widthofpart(unsigned int part,
                                         const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && part < suftabparts->numofparts);
  return suftabparts->components[part].widthofpart;
}

unsigned long gt_suftabparts_largest_width(const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL);
  return suftabparts->largestsuftabwidth;
}

unsigned long gt_suftabparts_largestsizemappedpartwise(
                                          const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL);
  return suftabparts->largestsizemappedpartwise;
}

unsigned int gt_suftabparts_numofparts(const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL);
  return suftabparts->numofparts;
}
