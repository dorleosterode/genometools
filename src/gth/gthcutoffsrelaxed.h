/*
  Copyright (c) 2003-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef GTHCUTOFFSRELAXED_H
#define GTHCUTOFFSRELAXED_H

#include "gth/gthtravalign.h"
#include "gth/sa.h"

/*
  This struct is used in gthdeterminecutoffs() to call
  gthlibtraversealignment().
*/
typedef struct {
  Cutoffs *cutoffs;
  bool breakforloop;
} Relaxedcutoffsdata;

void gt_initRelaxedcutoffsTravfunctions(Traversealignmentfunctions*);
void gt_initRelaxedcutoffsdata(Relaxedcutoffsdata*, Cutoffs*);

#endif
