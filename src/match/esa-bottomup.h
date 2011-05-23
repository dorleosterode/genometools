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

#ifndef ESA_BOTTOMUP_H
#define ESA_BOTTOMUP_H

#include "core/logger_api.h"
#include "core/error_api.h"
#include "match/esa-seqread.h"

typedef struct GtBUinfo GtBUinfo;
typedef struct GtBUstate GtBUstate;

int gt_esa_bottomup(Sequentialsuffixarrayreader *ssar,
                    GtBUinfo *(*allocateBUinfo)(GtBUstate *),
                    void(*freeBUinfo)(GtBUinfo *,GtBUstate *),
                    int (*processleafedge)(bool,unsigned long,
                                           unsigned long,
                                           GtBUinfo *,
                                           unsigned long,
                                           GtBUstate *,
                                           GtError *err),
                    int (*processbranchingedge)(bool firstsucc,
                                                unsigned long,
                                                unsigned long,
                                                GtBUinfo *,
                                                unsigned long,
                                                unsigned long,
                                                GtBUinfo *,
                                                GtBUstate *,
                                                GtError *),
                    GtBUstate *bustate,
                    GtError *err);

#endif
