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

/* include the files of gt scaffolder located in match as symbolic links */
#include "match/gt_scaffolder_bamparser.h"
#include "match/gt_scaffolder_parser.h"
#include "match/gt_scaffolder_graph.h"
#include "match/gt_scaffolder_algorithms.h"

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
  option = gt_option_new_double("rep_astat_cutoff", "minimal A-statistics for not repetitiv contigs",
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
  option = gt_option_new_string("astat", "A-statistics file in SGAs .astat format",
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

static int gt_scaffolder_runner(GT_UNUSED int argc, GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args, void *tool_arguments,
                                GtError *err)
{
  DistRecords *dist;
  GtScaffolderGraph *graph = NULL;
  bool astat_is_annotated = true;
  GtScaffolderArguments *arguments = tool_arguments;
  char *dist_file;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* check if distance file or bam file is given */
  if (gt_str_length(arguments->bam) > 0) {
    /* parse bam file and construct distance information */
    dist_file = "gt_scaffolder_bamparser_distance_records.de";
    /* initialize distance records */
    dist = gt_scaffolder_bamparser_init_dist_records();

    /* read paired information from bam file */
    had_err = gt_scaffolder_bamparser_read_paired_information(dist,
                                        gt_str_get(arguments->bam),
                                        arguments->bam_min_dist,
                                        arguments->bam_max_dist,
                                        arguments->bam_min_qual,
                                        arguments->bam_min_nof_pairs,
                                        arguments->min_contig_len,
                                        arguments->bam_min_align,
                                        err);
    if (had_err == 0) {
      /* print distance records */
      had_err = gt_scaffolder_bamparser_print_dist_records(dist,
                                                           dist_file, err);
      /* delete distance records */
      gt_scaffolder_bamparser_delete_dist_records(dist);
    }
  }

  if (had_err == 0) {
    if (gt_str_length(arguments->dist) > 0)
      dist_file = gt_str_get(arguments->dist);

    if (gt_str_length(arguments->astat) > 0)
      astat_is_annotated = false;

    /* create graph */
    had_err = gt_scaffolder_graph_new_from_file(&graph,
                                                gt_str_get(arguments->contigs),
                                                arguments->min_contig_len,
                                                dist_file,
                                                astat_is_annotated,
                                                err);
  }

  if (had_err == 0) {
    if (!astat_is_annotated) {
      /* load astatistics and copy number from file */
      had_err = gt_scaffolder_graph_mark_repeats(gt_str_get(arguments->astat),
                                                 graph,
                                                 arguments->rep_cp_cutoff,
                                                 arguments->rep_astat_cutoff,
                                                 err);
    }
  }
  if (had_err == 0) {
    /* mark polymorphic vertices, edges and inconsistent edges */
    gt_scaffolder_graph_filter(graph,
                               arguments->p_cutoff,
                               arguments->cp_cutoff,
                               arguments->overlap_cutoff);

    gt_scaffolder_makescaffold(graph);
    had_err = gt_scaffolder_graph_print(graph,
                                        "gt_scaffolder_scaffolds.dot",
                                        err);
  }

  if (had_err == 0) {
    GtArray *recs;
    GtAssemblyStatsCalculator *scaf_stats =
      gt_assembly_stats_calculator_new();

    recs = gt_scaffolder_graph_iterate_scaffolds(graph, scaf_stats);

    had_err = gt_scaffolder_graph_write_scaffold(recs,
                                                 "gt_scaffolder_scaffolds.scaf",
                                                 err);
    if (had_err == 0) {
      GtUword i;
      GtScaffolderGraphRecord *rec;
      GtLogger *logger;

      /* if sequence generation is wanted, it would be done here */

      for (i = 0; i < gt_array_size(recs); i++) {
        rec = *(GtScaffolderGraphRecord **) gt_array_get(recs, i);
        gt_scaffolder_graph_record_delete(rec);
      }
      gt_array_delete(recs);

      logger = gt_logger_new(true, "[scaffolder] ", stderr);
      gt_assembly_stats_calculator_nstat(scaf_stats, 50);
      gt_assembly_stats_calculator_show(scaf_stats, logger);

      gt_logger_delete(logger);
      gt_assembly_stats_calculator_delete(scaf_stats);
    }
  }

  gt_scaffolder_graph_delete(graph);

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
