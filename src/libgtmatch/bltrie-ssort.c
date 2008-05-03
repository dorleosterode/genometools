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

#include "libgtcore/arraydef.h"
#include "divmodmul.h"
#include "seqpos-def.h"
#include "spacedef.h"
#include "encseq-def.h"
#include "bltrie-ssort.h"
#include "stamp.h"

#undef SKDEBUG
#ifdef SKDEBUG
#include "sfx-cmpsuf.pr"

#define NODENUM(PTR)\
        ((PTR) == NULL ? 99UL\
                       : (unsigned long) ((PTR) - trierep->spaceBlindtrienode))
#endif

typedef struct Blindtrienode
{
  Seqpos depth;
  struct Blindtrienode *rightsibling;
  union
  {
    struct Blindtrienode *firstchild;
    Seqpos startpos;
  } either;
  bool isleaf;
  Uchar firstchar;
} Blindtrienode;

typedef Blindtrienode * Nodeptr;

DECLAREARRAYSTRUCT(Nodeptr);

struct Blindtrierep
{
  const Encodedsequence *encseq;
  Encodedsequencescanstate *esr1, *esr2;
  Readmode readmode;
  Seqpos totallength;
  Nodeptr root;
  unsigned long allocatedBlindtrienode,
                nextfreeBlindtrienode;
  Blindtrienode *spaceBlindtrienode;
  ArrayNodeptr stack;
};

static Nodeptr newBlindtrienode(Blindtrierep *trierep)
{
  assert(trierep->nextfreeBlindtrienode < trierep->allocatedBlindtrienode);
  return trierep->spaceBlindtrienode + trierep->nextfreeBlindtrienode++;
}

static Blindtrienode *makenewleaf(Blindtrierep *trierep,
                                  Seqpos startpos,
                                  Uchar firstchar)
{
  Blindtrienode *newleaf;

  newleaf = newBlindtrienode(trierep);
  newleaf->either.startpos = startpos;
  newleaf->depth = 0;
  newleaf->isleaf = true;
  newleaf->firstchar = firstchar;
  newleaf->rightsibling = NULL;
  return newleaf;
}

static Nodeptr makeroot(Blindtrierep *trierep,Seqpos startpos)
{
  Blindtrienode *root, *newleaf;
  Uchar firstchar;

  root = newBlindtrienode(trierep);
  root->depth = 0;
  root->firstchar = 0; /* undefined */
  root->rightsibling = NULL;
  root->isleaf = false;
  if (startpos >= trierep->totallength)
  {
    firstchar = (Uchar) SEPARATOR;
  } else
  {
    firstchar = getencodedchar(trierep->encseq, /* Random access */
                               startpos,
                               trierep->readmode);
    if (firstchar == (Uchar) WILDCARD)
    {
      firstchar = (Uchar) SEPARATOR;
    }
  }
  newleaf = makenewleaf(trierep,startpos,firstchar);
  root->either.firstchild = newleaf;
  return root;
}

static Nodeptr extractleafnode(Nodeptr head)
{
  assert(!head->isleaf);
  do
  {
    head = head->either.firstchild;
  } while (!head->isleaf);
  return head;
}

static int comparecharacters(Uchar oldchar,Uchar newchar)
{
  if (oldchar > newchar)
  {
    return 1;
  }
  if (oldchar < newchar || ISSPECIAL(oldchar))
  {
    return -1;
  }
  return 0;
}

static Nodeptr findsucc(Nodeptr node,Uchar newchar)
{
  int retval;
  for (;;)
  {
    retval = comparecharacters(node->firstchar,newchar);
    if (retval == 0)
    {              /* found branch corresponding to newchar */
      return node;
    }
    if (retval == 1)
    {               /* found branch which is already greater than newchar */
      return NULL;
    }
    node = node->rightsibling;
    if (node == NULL) /* no other branches: mismatch */
    {
      return NULL;
    }
  }
}

static Nodeptr findcompanion(Blindtrierep *trierep,Seqpos startpos)
{
  Uchar newchar;
  Nodeptr head, succ;

  trierep->stack.nextfreeNodeptr = 0;
  head = trierep->root;
  while (!head->isleaf)
  {
    STOREINARRAY (&trierep->stack, Nodeptr, 128, head);
    if (startpos + head->depth >= trierep->totallength)
    {
      newchar = (Uchar) SEPARATOR;
    } else
    {
      newchar = getencodedchar(trierep->encseq, /* Random access */
                               startpos + head->depth,
                               trierep->readmode);
      if (newchar == (Uchar) WILDCARD)
      {
        newchar = (Uchar) SEPARATOR;
      }
    }
    if (ISSPECIAL(newchar))
    {
      return extractleafnode(head);
    }
    succ = findsucc(head->either.firstchild,newchar);
    if (succ == NULL)
    {
      return extractleafnode(head);
    }
    head = succ;
  }
  STOREINARRAY (&trierep->stack, Nodeptr, 128, head);
  return head;
}

static void insertsuffixintoblindtrie(Blindtrierep *trierep,
                                      Nodeptr oldnode,
                                      Uchar mm_oldsuffix,
                                      Seqpos lcp,
                                      Uchar mm_newsuffix,
                                      Seqpos currentstartpos)
{
  Nodeptr newleaf, newnode, *pp;
  int retval;

  assert(ISSPECIAL(mm_oldsuffix) || ISSPECIAL(mm_newsuffix) ||
         mm_oldsuffix != mm_newsuffix || oldnode->isleaf ||
         oldnode->depth == lcp);

  /* insert a new node before node oldnode if necessary */
  if (oldnode->depth != lcp)
  {
    newnode = newBlindtrienode(trierep);
    newnode->firstchar = mm_oldsuffix;
    newnode->depth = oldnode->depth; /* newnode inherits depth+children */
    newnode->isleaf = oldnode->isleaf;
    newnode->either = oldnode->either;
    newnode->rightsibling = NULL;
    oldnode->depth = lcp;
    oldnode->isleaf = false;
    oldnode->either.firstchild = newnode; /* oldnode has newnode as only child*/
  }
  assert(oldnode->depth == lcp);

  /* search the position of s[n] among *h offsprings */
  pp = &(oldnode->either.firstchild);
  while ((*pp) != NULL)
  {
    retval = comparecharacters((*pp)->firstchar,mm_newsuffix);
    if (retval >= 0)
    {
      break;
    }
    pp = &((*pp)->rightsibling);
  }
  /* insert new node containing suf */
  newleaf = newBlindtrienode(trierep);
  newleaf->depth = 0;
  newleaf->isleaf = true;
  newleaf->firstchar = mm_newsuffix;
  newleaf->rightsibling = *pp;
  *pp = newleaf;
  newleaf->either.startpos = currentstartpos;
}

static Seqpos getlcp(Uchar *mm_oldsuffix,
                     Uchar *mm_newsuffix,
                     const Encodedsequence *encseq,
                     Encodedsequencescanstate *esr1,
                     Encodedsequencescanstate *esr2,
                     Readmode readmode,
                     Seqpos totallength,
                     Seqpos leafpos,
                     Seqpos currentstartpos)
{
  Seqpos idx1, idx2;
  Uchar cc1, cc2;

  initEncodedsequencescanstate(esr1,encseq,readmode,leafpos);
  initEncodedsequencescanstate(esr2,encseq,readmode,currentstartpos);
  for (idx1 = leafpos, idx2=currentstartpos;
       /* Nothing */;
       idx1++, idx2++)
  {
    if (idx1 < totallength)
    {
      cc1 = sequentialgetencodedchar(encseq,esr1,idx1,readmode);
      if (cc1 == (Uchar) WILDCARD)
      {
        cc1 = (Uchar) SEPARATOR;
      }
    } else
    {
      cc1 = (Uchar) SEPARATOR;
    }
    if (idx2 < totallength)
    {
      cc2 = sequentialgetencodedchar(encseq,esr2,idx2,readmode);
      if (cc2 == (Uchar) WILDCARD)
      {
        cc2 = (Uchar) SEPARATOR;
      }
    } else
    {
      cc2 = (Uchar) SEPARATOR;
    }
    if (comparecharacters(cc1,cc2) != 0)
    {
      *mm_oldsuffix = cc1;
      *mm_newsuffix = cc2;
      break;
    }
  }
  return idx1 - leafpos;
}

#define SETCURRENT(VAL)\
        currentnode = VAL;\
        currentnodeisleaf = (VAL)->isleaf ? true : false

static unsigned long enumeratetrieleaves (Seqpos *suffixtable,
                                          Seqpos *lcpsubtab,
                                          Blindtrierep *trierep,
                                          Seqpos offset)
{
  bool readyforpop = false, currentnodeisleaf;
  Nodeptr currentnode, siblval, lcpnode = trierep->root;
  unsigned long nextfree = 0;

  trierep->stack.nextfreeNodeptr = 0;
  STOREINARRAY (&trierep->stack, Nodeptr, 128, trierep->root);
  SETCURRENT(trierep->root->either.firstchild);
  for (;;)
  {
    if (currentnodeisleaf)
    {
      assert (currentnode->either.startpos >= offset);
      if (lcpsubtab != NULL && nextfree > 0)
      {
        lcpsubtab[nextfree] = lcpnode->depth + offset;
      }
      suffixtable[nextfree++] = currentnode->either.startpos - offset;
      siblval = currentnode->rightsibling;
      if (siblval == NULL)
      {
        readyforpop = true;
        currentnodeisleaf = false; /* STATE 1 */
      } else
      {
        SETCURRENT (siblval);  /* current comes from brother */
        lcpnode = trierep->stack.spaceNodeptr[trierep->stack.nextfreeNodeptr-1];
      }
    } else
    {
      if (readyforpop)
      {
        if (trierep->stack.nextfreeNodeptr == 1UL)
        {
          break;
        }
        trierep->stack.nextfreeNodeptr--;
        siblval = trierep->stack.spaceNodeptr[
                           trierep->stack.nextfreeNodeptr]->rightsibling;
        if (siblval != NULL)
        {
          SETCURRENT (siblval);        /* current comes from brother */
          lcpnode = trierep->stack.spaceNodeptr[
                             trierep->stack.nextfreeNodeptr - 1];
          readyforpop = false;
        }
      } else
      {
        STOREINARRAY (&trierep->stack, Nodeptr, 128, currentnode);
        SETCURRENT (currentnode->either.firstchild);
      }
    }
  }
  return nextfree;
}

Blindtrierep *newBlindtrierep(unsigned long numofsuffixes,
                              const Encodedsequence *encseq,
                              Readmode readmode)
{
  Blindtrierep *trierep;

  ALLOCASSIGNSPACE(trierep,NULL,Blindtrierep,1);
  trierep->allocatedBlindtrienode = MULT2(numofsuffixes + 1) + 1;
  ALLOCASSIGNSPACE(trierep->spaceBlindtrienode,NULL,Blindtrienode,
                   trierep->allocatedBlindtrienode);
  trierep->nextfreeBlindtrienode = 0;
  trierep->encseq = encseq;
  trierep->readmode = readmode;
  trierep->root = NULL;
  trierep->esr1 = newEncodedsequencescanstate();
  trierep->esr2 = newEncodedsequencescanstate();
  trierep->totallength = getencseqtotallength(encseq);
  INITARRAY (&trierep->stack, Nodeptr);
  return trierep;
}

void freeBlindtrierep(Blindtrierep **trierep)
{
  STAMP;
  FREESPACE((*trierep)->spaceBlindtrienode);
  STAMP;
  FREEARRAY(&(*trierep)->stack, Nodeptr);
  STAMP;
  freeEncodedsequencescanstate(&((*trierep)->esr1));
  STAMP;
  freeEncodedsequencescanstate(&((*trierep)->esr2));
  STAMP;
  FREESPACE(*trierep);
  STAMP;
}

#ifdef SKDEBUG
static void checkcurrentblindtrie(Blindtrierep *trierep,
                                  Seqpos offset)
{
  Seqpos suffixtable[6];
  unsigned long idx, numofsuffixes;
  Seqpos maxcommon;
  int retval;

  numofsuffixes = enumeratetrieleaves (&suffixtable[0], NULL, trierep, offset);
  for (idx=1UL; idx < numofsuffixes; idx++)
  {
    maxcommon = 0;
    retval = comparetwostringsgeneric(trierep->encseq,
                                      ISDIRREVERSE(trierep->readmode)
                                        ? false : true,
                                      ISDIRCOMPLEMENT(trierep->readmode)
                                        ? true : false,
                                      &maxcommon,
                                      suffixtable[idx-1],
                                      suffixtable[idx],
                                      0);
    if (retval >= 0)
    {
      fprintf(stderr,"retval = %d, maxcommon = %u for idx = %lu\n",
              retval,maxcommon,idx);
      assert(retval < 0);
    }
  }
}

static void showleaf(const Blindtrierep *trierep,unsigned int level,
                     Nodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  assert(current != NULL);
  printf("Leaf(add=%lu,firstchar=%u,startpos=" FormatSeqpos
         ",rightsibling=%lu)\n",
         NODENUM(current),
         (unsigned int) current->firstchar,
         PRINTSeqposcast(current->either.startpos),
         NODENUM(current->rightsibling));
}

static void showintern(const Blindtrierep *trierep,unsigned int level,
                       Nodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  assert(current != NULL);
  printf("Intern(add=%lu,firstchar=%u,depth=" FormatSeqpos
         ",firstchild=%lu,rightsibling=%lu)\n",
          NODENUM(current),
          (unsigned int) current->firstchar,
          PRINTSeqposcast(current->depth),
          NODENUM(current->either.firstchild),
          NODENUM(current->rightsibling));
}

static void showblindtrie2(const Blindtrierep *trierep,
                           unsigned int level,
                           Nodeptr node)
{
  Nodeptr current;

  for (current = node->either.firstchild;
       current != NULL;
       current = current->rightsibling)
  {
    if (current->isleaf)
    {
      showleaf(trierep,level,current);
    } else
    {
      showintern(trierep,level,current);
      showblindtrie2(trierep,level+1,current);
    }
  }
}

static void showblindtrie(const Blindtrierep *trierep)
{
  showblindtrie2(trierep,0,trierep->root);
}
#endif

static int suffixcompare(const void *a, const void *b)
{
  assert(*((Seqpos *) a) != *((Seqpos *) b));
  if (*((Seqpos *) a) < *((Seqpos *) b))
  {
    return -1;
  }
  return 1;
}

void blindtriesuffixsort(Blindtrierep *trierep,
                         Seqpos *suffixtable,
                         Seqpos *lcpsubtab,
                         unsigned long numberofsuffixes,
                         Seqpos offset)
{
  unsigned long i, t;
  Nodeptr leafinsubtree, currentnode;
  Seqpos lcp;
  Uchar mm_oldsuffix, mm_newsuffix;

  qsort(suffixtable,(size_t) numberofsuffixes, sizeof (Seqpos), suffixcompare);
  trierep->nextfreeBlindtrienode = 0;
  trierep->root = makeroot(trierep,suffixtable[0] + offset);
#ifdef SKDEBUG
  printf("insert suffixes at offset " FormatSeqpos ":\n",
          PRINTSeqposcast(offset));
  for (i=0; i < numberofsuffixes; i++)
  {
    printf(FormatSeqpos " ",PRINTSeqposcast(suffixtable[i] + offset));
  }
  printf("\nstep 0\n");
  showblindtrie(trierep);
#endif
  for (i=1UL; i < numberofsuffixes; i++)
  {
    leafinsubtree = findcompanion(trierep, suffixtable[i] + offset);
    assert(leafinsubtree->isleaf);
    lcp = getlcp(&mm_oldsuffix,&mm_newsuffix,
                 trierep->encseq,
                 trierep->esr1,
                 trierep->esr2,
                 trierep->readmode,
                 trierep->totallength,
                 leafinsubtree->either.startpos,
                 suffixtable[i] + offset);
    currentnode = trierep->root;
    for (t=0;t<trierep->stack.nextfreeNodeptr;t++)
    {
      currentnode = trierep->stack.spaceNodeptr[t];
      if (currentnode->isleaf || currentnode->depth >= lcp)
      {
        break;
      }
    }
    insertsuffixintoblindtrie(trierep,
                              currentnode,
                              mm_oldsuffix,
                              lcp,
                              mm_newsuffix,
                              suffixtable[i] + offset);
#ifdef SKDEBUG
    printf("step %lu\n",i);
    showblindtrie(trierep);
    checkcurrentblindtrie(trierep,offset);
#endif
  }
  (void) enumeratetrieleaves (suffixtable, lcpsubtab, trierep, offset);
}
