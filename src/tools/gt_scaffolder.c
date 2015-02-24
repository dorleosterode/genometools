/*
  Copyright (c) 2015 Dorle Osterode <9osterod@informatik.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/unused_api.h"
#include "tools/gt_scaffolder.h"

typedef struct {
  /* options for gt scaffolder */
  GtUword min_contig_len;
  double rep_cp_cutoff;
  double rep_astat_cutoff;
  double p_cutoff;
  double cp_cutoff;
  GtUword overlap_cutoff;
  /* files for gt scaffolder */
  GtStr *contigs;
  GtStr *dist;
  GtStr *astat;
  GtStr *spm;
  /* options for bamparser */
  GtUword bam_min_qual;
  GtUword bam_min_nof_pairs;
  GtUword bam_min_ref_length;
  GtWord bam_min_dist;
  GtWord bam_max_dist;
  GtUword bam_min_align;
  /* files for bamparser */
  GtStr *bam;
} GtScaffolderArguments;

static void* gt_scaffolder_arguments_new(void)
{
  GtScaffolderArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->contigs = gt_str_new();
  arguments->dist = gt_str_new();
  arguments->astat = gt_str_new();
  arguments->spm = gt_str_new();
  arguments->bam = gt_str_new();
  return arguments;
}

static void gt_scaffolder_arguments_delete(void *tool_arguments)
{
  GtScaffolderArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->contigs);
    gt_str_delete(arguments->dist);
    gt_str_delete(arguments->astat);
    gt_str_delete(arguments->spm);
    gt_str_delete(arguments->bam);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_scaffolder_option_parser_new(void *tool_arguments)
{
  GtScaffolderArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  GtOption *dist;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]", /* XXX */
                            "Constructs scaffolds of the contigs."); /* XXX */

  /* options for gt scaffolder */

  /* - min contig length */
  option = gt_option_new_uword("min_contig_len", "minimal contig length for used contigs",
                               &arguments->min_contig_len, 200);
  gt_option_parser_add_option(op, option);


  /* - copy num cutoff for repeats */
  option = gt_option_new_double("rep_cp_cutoff", "minimal copy number for not repetitiv contigs",
                                &arguments->rep_cp_cutoff, 0.3);
  gt_option_parser_add_option(op, option);

  /* - astat cutoff */
  option = gt_option_new_double("rep_astat_cutoff", "minimal astatistic for not repetitiv contigs",
                               &arguments->rep_astat_cutoff, 20.0);
  gt_option_parser_add_option(op, option);

  /* - probability cutoff */
  option = gt_option_new_probability("p_cutoff", "probability cutoff for polymorphic vertices",
                                &arguments->p_cutoff, 0.01);
  gt_option_parser_add_option(op, option);

  /* - copy num cutoff */
  option = gt_option_new_double("cp_cutoff", "copy num cutoff for polymorphic vertices",
                                &arguments->cp_cutoff, 1.5);
  gt_option_parser_add_option(op, option);

  /* - overlap cutoff */
  option = gt_option_new_uword("overlap_cutoff", "overlap cutoff for inconsistent edges",
                               &arguments->overlap_cutoff, 400);
  gt_option_parser_add_option(op, option);

  /* files for gt scaffolder */

  /* - contigs file */
  option = gt_option_new_string("contigs", "contigs in FASTA format",
                            arguments->contigs, NULL);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /* - distance file */
  dist = gt_option_new_string("dist", "distanceinformation in ABySS .de format",
			      arguments->dist, NULL);
  gt_option_parser_add_option(op, dist);

  /* - astat file */
  option = gt_option_new_string("astat", "astatistic file in SGAs .astat format",
                            arguments->astat, NULL);
  gt_option_parser_add_option(op, option);

  /* - spm file */
  option = gt_option_new_string("spm", "basename of spm-representation of used string-graph",
                            arguments->spm, NULL);
  gt_option_parser_add_option(op, option);

  /* options for bamparser */

  /* - bam min qual */
  option = gt_option_new_uword("bam_min_qual", "minimal quality in bam file",
			       &arguments->bam_min_qual, 10);
  gt_option_parser_add_option(op, option);

  /* - bam min nof pairs */
  option = gt_option_new_uword("bam_min_nof_pairs", "minimal number of pairs",
			       &arguments->bam_min_nof_pairs, 10);
  gt_option_parser_add_option(op, option);

  /* - bam min ref length */
  option = gt_option_new_uword("bam_min_ref_length", "minimal reference length",
			       &arguments->bam_min_qual, 200);
  gt_option_parser_add_option(op, option);

  /* - bam min dist */
  option = gt_option_new_word("bam_min_dist", "minimal distance of contigs",
			      &arguments->bam_min_dist, -99);
  gt_option_parser_add_option(op, option);

  /* - bam max dist */
  option = gt_option_new_word("bam_max_dist", "maximal distance of contigs",
			      &arguments->bam_max_dist, GT_WORD_MAX);
  gt_option_parser_add_option(op, option);

  /* - bam min align */
  option = gt_option_new_uword("bam_min_align", "minimal alignment length",
			       &arguments->bam_min_align, 100);
  gt_option_parser_add_option(op, option);

  /* - bam file */
  option = gt_option_new_string("bam", "bam file containing alignments of paired reads to the contigs",
                                arguments->bam, NULL);
  gt_option_is_mandatory_either(dist, option);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_scaffolder_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtScaffolderArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  /* XXX: do some checking after the option have been parsed (usally this is not
     necessary and this function can be removed completely). */
  /* if (gt_str_length(arguments->str_option_scaffolder)) */
  /*   printf("%s\n", gt_str_get(arguments->str_option_scaffolder)); */

  return had_err;
}

static int gt_scaffolder_runner(int argc, const char **argv, int parsed_args,
                              void *tool_arguments, GT_UNUSED GtError *err)
{
  GtScaffolderArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* XXX */
  printf("argc=%d, parsed_args=%d\n", argc, parsed_args);
  printf("argv[0]=%s\n", argv[0]);

  return had_err;
}

GtTool* gt_scaffolder(void)
{
  return gt_tool_new(gt_scaffolder_arguments_new,
                     gt_scaffolder_arguments_delete,
                     gt_scaffolder_option_parser_new,
                     gt_scaffolder_arguments_check,
                     gt_scaffolder_runner);
}
