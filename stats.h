#ifndef STATS_H
#define STATS_H

#include <stdio.h>

typedef struct {
    long   n_genes;
    long   n_transcripts;
    long   n_exons;
    long   n_introns;
    double mean_exons_per_gene;
    double mean_introns_per_gene;
    long   n_single_exon_genes;
    double pct_single_exon_genes;
    long   n_single_exon_transcripts;
    double pct_single_exon_transcripts;
} GtfExonIntronStats;

/*
 * Read all exon features from fp and populate *stats.
 * Returns 0 on success, -1 on allocation failure.
 * n_introns is computed as sum of (exons_per_transcript - 1).
 * A single-exon gene is one where every transcript has exactly 1 exon.
 */
int  gtf_exon_intron_stats(FILE *fp, GtfExonIntronStats *stats);

void gtf_print_exon_intron_stats(const GtfExonIntronStats *stats);

#endif /* STATS_H */
