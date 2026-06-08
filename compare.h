#ifndef COMPARE_H
#define COMPARE_H

#include <stdio.h>
#include "gtf.h"

/* One intron: half-open on the GTF coordinate system (1-based closed input). */
typedef struct {
    long start;
    long end;
} IntronSpan;

/*
 * One transcript, derived from its exon features.
 * introns[] is a heap array of length n_introns, sorted by start.
 * src = 1 (GTF A) or 2 (GTF B), set by cmp_build_transcripts.
 * cmp_class is filled in by cmp_compare:
 *   '=' — exact intron chain match (all splice sites identical)
 *   'j' — partial match (at least one shared intron, but chains differ)
 *   'o' — positional overlap, no shared introns
 *   'u' — no overlap with any transcript from the other file
 */
typedef struct {
    char        seqname[GTF_FIELD_MAX];
    char        tx_id[GTF_FIELD_MAX];
    char        gene_id[GTF_FIELD_MAX];
    char        strand;
    long        start;       /* leftmost exon start  */
    long        end;         /* rightmost exon end   */
    int         n_exons;
    int         n_introns;
    IntronSpan *introns;     /* heap array [n_introns], sorted by start */
    int         src;
    char        cmp_class;
} Transcript;

/*
 * Comparison results — symmetric: both files get full classification counts.
 * Indices [1] and [2] correspond to GTF A and GTF B respectively.
 */
typedef struct {
    long n_tx[3];             /* [1]=GTF A, [2]=GTF B                  */
    long n_introns[3];        /* unique introns per file                */
    long n_shared_introns;    /* introns present in both files          */
    long n_exact[3];          /* (=) full intron chain match            */
    long n_partial[3];        /* (j) at least one shared intron         */
    long n_overlap[3];        /* (o) overlap, zero shared introns       */
    long n_unique[3];         /* (u) no overlap with opposite file      */
} CmpStats;

/*
 * Parse all exon features from fp and build one Transcript per tx_id.
 * src must be 1 or 2. *out_n is set to the number of transcripts.
 * Returns a heap array the caller must free with cmp_free_transcripts(),
 * or NULL on allocation failure (out_n is set to 0).
 */
Transcript *cmp_build_transcripts(FILE *fp, int src, size_t *out_n);

/* Free a transcript array and all nested intron arrays. */
void cmp_free_transcripts(Transcript *txs, size_t n);

/*
 * Compare two transcript arrays (built from different GTF files).
 * Both arrays are sorted in place by (seqname, strand, start, end).
 * Results are written into *stats.
 * Returns 0 on success, -1 on allocation failure.
 */
int cmp_compare(Transcript *a, size_t na,
                Transcript *b, size_t nb,
                CmpStats *stats);

/* Print a human-readable summary of *stats to stdout. */
void cmp_print_stats(const CmpStats *stats,
                     const char *label_a, const char *label_b);

#endif /* COMPARE_H */
