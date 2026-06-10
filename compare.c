#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gtf.h"
#include "compare.h"

/* ------------------------------------------------------------------ */
/* Internal types                                                       */
/* ------------------------------------------------------------------ */

/* One exon row collected during phase 1 of cmp_build_transcripts */
typedef struct {
    char seqname[GTF_FIELD_MAX];
    char tx_id[GTF_FIELD_MAX];
    char gene_id[GTF_FIELD_MAX];
    char strand;
    long start;
    long end;
} ExonRow;

/* Flattened intron key used for global deduplication */
typedef struct {
    char seqname[GTF_FIELD_MAX];
    char strand;
    long start;
    long end;
} IntronKey;

/* ------------------------------------------------------------------ */
/* Attribute extraction                                                 */
/* ------------------------------------------------------------------ */

static int attr_get(const char *attrs, const char *key,
                    char *val, size_t val_sz)
{
    size_t klen = strlen(key);
    const char *p = attrs;
    while (*p) {
        while (*p == ' ' || *p == '\t') p++;
        if (strncmp(p, key, klen) == 0 &&
            (p[klen] == ' ' || p[klen] == '\t')) {
            p += klen;
            while (*p == ' ' || *p == '\t') p++;
            if (*p == '"') {
                p++;
                size_t i = 0;
                while (*p && *p != '"' && i < val_sz - 1)
                    val[i++] = *p++;
                val[i] = '\0';
            } else {
                size_t i = 0;
                while (*p && *p != ';' && *p != ' ' && *p != '\t'
                       && i < val_sz - 1)
                    val[i++] = *p++;
                val[i] = '\0';
            }
            return 1;
        }
        while (*p && *p != ';') p++;
        if (*p == ';') p++;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Comparators                                                          */
/* ------------------------------------------------------------------ */

static int exon_cmp(const void *a, const void *b)
{
    const ExonRow *ea = (const ExonRow *)a;
    const ExonRow *eb = (const ExonRow *)b;
    int c = strcmp(ea->tx_id, eb->tx_id);
    if (c) return c;
    if (ea->start < eb->start) return -1;
    if (ea->start > eb->start) return  1;
    return 0;
}

/* Used to sort Transcript arrays in place */
static int tx_cmp(const void *a, const void *b)
{
    const Transcript *ta = (const Transcript *)a;
    const Transcript *tb = (const Transcript *)b;
    int c = strcmp(ta->seqname, tb->seqname);
    if (c) return c;
    if (ta->strand != tb->strand)
        return (ta->strand < tb->strand) ? -1 : 1;
    if (ta->start < tb->start) return -1;
    if (ta->start > tb->start) return  1;
    if (ta->end   < tb->end)   return -1;
    if (ta->end   > tb->end)   return  1;
    return 0;
}

/* Used for pointer array over the merged A+B transcripts in sweep */
static int tx_ptr_cmp(const void *a, const void *b)
{
    const Transcript *ta = *(const Transcript * const *)a;
    const Transcript *tb = *(const Transcript * const *)b;
    int c = strcmp(ta->seqname, tb->seqname);
    if (c) return c;
    if (ta->strand != tb->strand)
        return (ta->strand < tb->strand) ? -1 : 1;
    if (ta->start < tb->start) return -1;
    if (ta->start > tb->start) return  1;
    if (ta->end   < tb->end)   return -1;
    if (ta->end   > tb->end)   return  1;
    return 0;
}

static int intron_key_cmp(const void *a, const void *b)
{
    const IntronKey *ia = (const IntronKey *)a;
    const IntronKey *ib = (const IntronKey *)b;
    int c = strcmp(ia->seqname, ib->seqname);
    if (c) return c;
    if (ia->strand != ib->strand)
        return (ia->strand < ib->strand) ? -1 : 1;
    if (ia->start < ib->start) return -1;
    if (ia->start > ib->start) return  1;
    if (ia->end   < ib->end)   return -1;
    if (ia->end   > ib->end)   return  1;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Build transcript array from a GTF file                              */
/* ------------------------------------------------------------------ */

Transcript *cmp_build_transcripts(FILE *fp, int src, size_t *out_n)
{
    *out_n = 0;

    /* Phase 1 — collect exon rows */
    size_t ecap = 4096, en = 0;
    ExonRow *exons = malloc(ecap * sizeof(ExonRow));
    if (!exons) return NULL;

    char line[65536];
    GtfRawRecord rec;

    while (fgets(line, sizeof(line), fp)) {
        if (gtf_parse_line(line, &rec) != 0) continue;
        if (strcmp(rec.feature, "exon") != 0) continue;

        if (en >= ecap) {
            ecap *= 2;
            ExonRow *tmp = realloc(exons, ecap * sizeof(ExonRow));
            if (!tmp) { free(exons); return NULL; }
            exons = tmp;
        }

        ExonRow *e = &exons[en];
        e->strand = rec.strand;
        e->start  = rec.start;
        e->end    = rec.end;
        strncpy(e->seqname,  rec.seqname,  GTF_FIELD_MAX - 1);
        e->seqname[GTF_FIELD_MAX - 1]  = '\0';
        if (!attr_get(rec.attributes, "transcript_id", e->tx_id,  GTF_FIELD_MAX))
            e->tx_id[0] = '\0';
        if (!attr_get(rec.attributes, "gene_id",       e->gene_id, GTF_FIELD_MAX))
            e->gene_id[0] = '\0';
        en++;
    }

    if (en == 0) { free(exons); return NULL; }

    /* Phase 2 — sort by (tx_id, start), group into transcripts */
    qsort(exons, en, sizeof(ExonRow), exon_cmp);

    size_t tcap = 1024, tn = 0;
    Transcript *txs = malloc(tcap * sizeof(Transcript));
    if (!txs) { free(exons); return NULL; }

    for (size_t i = 0; i < en; ) {
        /* find the end of this transcript's exon group */
        size_t j = i + 1;
        while (j < en && strcmp(exons[j].tx_id, exons[i].tx_id) == 0) j++;

        if (tn >= tcap) {
            tcap *= 2;
            Transcript *tmp = realloc(txs, tcap * sizeof(Transcript));
            if (!tmp) {
                for (size_t k = 0; k < tn; k++) free(txs[k].introns);
                free(txs); free(exons);
                return NULL;
            }
            txs = tmp;
        }

        Transcript *t = &txs[tn];
        t->strand    = exons[i].strand;
        t->src       = src;
        t->cmp_class = 'u';
        t->n_exons   = (int)(j - i);
        t->n_introns = (int)(j - i) - 1;
        strncpy(t->seqname,  exons[i].seqname,  GTF_FIELD_MAX - 1);
        t->seqname[GTF_FIELD_MAX - 1]  = '\0';
        strncpy(t->tx_id,    exons[i].tx_id,    GTF_FIELD_MAX - 1);
        t->tx_id[GTF_FIELD_MAX - 1]    = '\0';
        strncpy(t->gene_id,  exons[i].gene_id,  GTF_FIELD_MAX - 1);
        t->gene_id[GTF_FIELD_MAX - 1]  = '\0';

        /* transcript span = min start / max end across all its exons */
        t->start = exons[i].start;
        t->end   = exons[i].end;
        for (size_t k = i + 1; k < j; k++) {
            if (exons[k].start < t->start) t->start = exons[k].start;
            if (exons[k].end   > t->end)   t->end   = exons[k].end;
        }

        /* intron array — gap between consecutive sorted exons */
        if (t->n_introns > 0) {
            t->introns = malloc(t->n_introns * sizeof(IntronSpan));
            if (!t->introns) {
                for (size_t k = 0; k < tn; k++) free(txs[k].introns);
                free(txs); free(exons);
                return NULL;
            }
            for (int k = 0; k < t->n_introns; k++) {
                t->introns[k].start = exons[i + k].end + 1;
                t->introns[k].end   = exons[i + k + 1].start - 1;
            }
        } else {
            t->introns = NULL;
        }

        tn++;
        i = j;
    }

    free(exons);
    *out_n = tn;
    return txs;
}

void cmp_free_transcripts(Transcript *txs, size_t n)
{
    if (!txs) return;
    for (size_t i = 0; i < n; i++) free(txs[i].introns);
    free(txs);
}

/* ------------------------------------------------------------------ */
/* Within-locus comparison helpers                                      */
/* ------------------------------------------------------------------ */

/*
 * Count introns that are identical (same start AND end) in both transcripts.
 * Both intron arrays are sorted by start, so this is a linear merge.
 */
static int count_shared_introns(const Transcript *a, const Transcript *b)
{
    int i = 0, j = 0, shared = 0;
    while (i < a->n_introns && j < b->n_introns) {
        long as = a->introns[i].start, ae = a->introns[i].end;
        long bs = b->introns[j].start, be = b->introns[j].end;
        if (as == bs && ae == be) {
            shared++; i++; j++;
        } else if (as < bs || (as == bs && ae < be)) {
            i++;
        } else {
            j++;
        }
    }
    return shared;
}

static int tx_overlap(const Transcript *a, const Transcript *b)
{
    return strcmp(a->seqname, b->seqname) == 0
        && a->strand == b->strand
        && a->start  <= b->end
        && b->start  <= a->end;
}

/* Returns 3 for '=', 2 for 'j', 1 for 'o', 0 for 'u'. */
static int class_rank(char c)
{
    if (c == '=') return 3;
    if (c == 'j') return 2;
    if (c == 'o') return 1;
    return 0;
}

/*
 * Classify every transcript in all[0..n).
 * For each transcript, scan all other-source transcripts in the locus,
 * pick the best class, and accumulate into stats.
 */
static void compare_locus(Transcript **all, size_t n, CmpStats *stats)
{
    for (size_t i = 0; i < n; i++) {
        Transcript *t = all[i];
        char best = 'u';

        for (size_t j = 0; j < n; j++) {
            Transcript *o = all[j];
            if (o->src == t->src) continue;
            if (!tx_overlap(t, o)) continue;

            int shared = count_shared_introns(t, o);
            char cls;
            if (t->n_introns > 0
                    && shared == t->n_introns
                    && shared == o->n_introns)
                cls = '=';
            else if (shared > 0)
                cls = 'j';
            else
                cls = 'o';

            if (class_rank(cls) > class_rank(best))
                best = cls;
        }

        t->cmp_class = best;
        int s = t->src;
        if      (best == '=') stats->n_exact[s]++;
        else if (best == 'j') stats->n_partial[s]++;
        else if (best == 'o') stats->n_overlap[s]++;
        else                  stats->n_unique[s]++;
    }
}

/* ------------------------------------------------------------------ */
/* Global intron deduplication for intron-level stats                  */
/* ------------------------------------------------------------------ */

/*
 * Flatten all introns from txs[0..n) into a sorted, deduplicated
 * IntronKey array.  *out_count is set to the number of unique introns.
 * Returns heap array (caller frees), or NULL on allocation failure
 * (*out_count set to -1).
 */
static IntronKey *dedup_introns(const Transcript *txs, size_t n,
                                 long *out_count)
{
    long total = 0;
    for (size_t i = 0; i < n; i++) total += txs[i].n_introns;

    if (total == 0) { *out_count = 0; return NULL; }

    IntronKey *arr = malloc(total * sizeof(IntronKey));
    if (!arr) { *out_count = -1; return NULL; }

    long k = 0;
    for (size_t i = 0; i < n; i++) {
        for (int m = 0; m < txs[i].n_introns; m++) {
            strncpy(arr[k].seqname, txs[i].seqname, GTF_FIELD_MAX - 1);
            arr[k].seqname[GTF_FIELD_MAX - 1] = '\0';
            arr[k].strand = txs[i].strand;
            arr[k].start  = txs[i].introns[m].start;
            arr[k].end    = txs[i].introns[m].end;
            k++;
        }
    }

    qsort(arr, total, sizeof(IntronKey), intron_key_cmp);

    /* deduplicate in place */
    long uniq = 0;
    for (long x = 0; x < total; x++) {
        if (x == 0 || intron_key_cmp(&arr[x], &arr[x - 1]) != 0)
            arr[uniq++] = arr[x];
    }

    *out_count = uniq;
    return arr;
}

/* ------------------------------------------------------------------ */
/* Main comparison                                                      */
/* ------------------------------------------------------------------ */

int cmp_compare(Transcript *a, size_t na,
                Transcript *b, size_t nb,
                CmpStats *stats)
{
    memset(stats, 0, sizeof(*stats));
    stats->n_tx[1] = (long)na;
    stats->n_tx[2] = (long)nb;

    /* --- intron-level stats --- */
    long n_uniq_a, n_uniq_b;
    IntronKey *ikeys_a = dedup_introns(a, na, &n_uniq_a);
    IntronKey *ikeys_b = dedup_introns(b, nb, &n_uniq_b);
    if (n_uniq_a < 0 || n_uniq_b < 0) {
        free(ikeys_a); free(ikeys_b); return -1;
    }
    stats->n_introns[1] = n_uniq_a;
    stats->n_introns[2] = n_uniq_b;

    /* count shared introns via merge scan on the two sorted unique arrays */
    for (long i = 0, j = 0; i < n_uniq_a && j < n_uniq_b; ) {
        int c = intron_key_cmp(&ikeys_a[i], &ikeys_b[j]);
        if      (c == 0) { stats->n_shared_introns++; i++; j++; }
        else if (c  < 0) { i++; }
        else             { j++; }
    }
    free(ikeys_a);
    free(ikeys_b);

    /* --- transcript-level sweep --- */

    /* Sort both arrays by genomic position for the sweep line */
    if (na) qsort(a, na, sizeof(Transcript), tx_cmp);
    if (nb) qsort(b, nb, sizeof(Transcript), tx_cmp);

    /*
     * Build a pointer array over all transcripts and sort it.
     * Sorting pointers avoids copying large structs and keeps the
     * original arrays intact for free later.
     */
    size_t ntotal = na + nb;
    if (ntotal == 0) return 0;

    Transcript **all = malloc(ntotal * sizeof(Transcript *));
    if (!all) return -1;
    for (size_t i = 0; i < na; i++) all[i]      = &a[i];
    for (size_t i = 0; i < nb; i++) all[na + i] = &b[i];
    qsort(all, ntotal, sizeof(Transcript *), tx_ptr_cmp);

    /*
     * Sweep line: advance through sorted transcripts extending a
     * locus as long as the next transcript overlaps (same seqname,
     * same strand, start <= current locus end).
     * When the locus closes, run compare_locus on it.
     */
    size_t lo = 0;
    while (lo < ntotal) {
        const char *lseq    = all[lo]->seqname;
        char        lstrand = all[lo]->strand;
        long        lend    = all[lo]->end;
        size_t hi = lo + 1;

        while (hi < ntotal
               && strcmp(all[hi]->seqname, lseq)  == 0
               && all[hi]->strand == lstrand
               && all[hi]->start  <= lend) {
            if (all[hi]->end > lend) lend = all[hi]->end;
            hi++;
        }

        compare_locus(all + lo, hi - lo, stats);
        lo = hi;
    }

    free(all);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Output                                                               */
/* ------------------------------------------------------------------ */

static void print_class_block(const CmpStats *stats, int s, const char *label)
{
    long n = stats->n_tx[s];
    if (n == 0) return;

    double pct_e = 100.0 * stats->n_exact[s]   / n;
    double pct_j = 100.0 * stats->n_partial[s] / n;
    double pct_o = 100.0 * stats->n_overlap[s] / n;
    double pct_u = 100.0 * stats->n_unique[s]  / n;

    int W = 26;
    printf("=== Transcript classification (%s) ===\n", label);
    printf("  %-*s: %8ld\n",              W, label,                  n);
    printf("  %-*s: %8ld  (%5.1f%%)\n",  W, "Exact match      (=)", stats->n_exact[s],   pct_e);
    printf("  %-*s: %8ld  (%5.1f%%)\n",  W, "Partial match    (j)", stats->n_partial[s], pct_j);
    printf("  %-*s: %8ld  (%5.1f%%)\n",  W, "Overlap only     (o)", stats->n_overlap[s], pct_o);
    printf("  %-*s: %8ld  (%5.1f%%)\n",  W, "Unique           (u)", stats->n_unique[s],  pct_u);
    printf("\n");
}

void cmp_print_stats(const CmpStats *stats,
                     const char *label_a, const char *label_b)
{
    printf("Structural GTF comparison\n");
    printf("  GTF1 = %s\n", label_a);
    printf("  GTF2 = %s\n\n", label_b);

    print_class_block(stats, 1, "GTF1");
    print_class_block(stats, 2, "GTF2");

    long ia = stats->n_introns[1];
    long ib = stats->n_introns[2];
    long sh = stats->n_shared_introns;

    printf("=== Intron-level ===\n");
    printf("  %-26s: %8ld\n", "GTF1 unique introns", ia);
    printf("  %-26s: %8ld\n", "GTF2 unique introns", ib);
    if (ia > 0 || ib > 0) {
        double pct_a = ia > 0 ? 100.0 * sh / ia : 0.0;
        double pct_b = ib > 0 ? 100.0 * sh / ib : 0.0;
        printf("  %-26s: %8ld  (%.1f%% of GTF1, %.1f%% of GTF2)\n",
               "Shared", sh, pct_a, pct_b);
    }
    printf("\n");

    /* Sensitivity / precision summary */
    long na = stats->n_tx[1], nb = stats->n_tx[2];
    if (na > 0 && nb > 0) {
        printf("  Sensitivity (GTF1 exact / GTF1 total) : %.1f%%\n",
               100.0 * stats->n_exact[1] / na);
        printf("  Precision   (GTF2 exact / GTF2 total) : %.1f%%\n",
               100.0 * stats->n_exact[2] / nb);
    }
}
