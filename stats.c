#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gtf.h"
#include "stats.h"

#define INIT_CAP 4096

typedef struct {
    char gene_id[GTF_FIELD_MAX];
    char tx_id[GTF_FIELD_MAX];
    long exon_count;
} TxEntry;

static int cmp_by_tx(const void *a, const void *b)
{
    return strcmp(((const TxEntry *)a)->tx_id, ((const TxEntry *)b)->tx_id);
}

static int cmp_by_gene(const void *a, const void *b)
{
    int c = strcmp(((const TxEntry *)a)->gene_id, ((const TxEntry *)b)->gene_id);
    if (c) return c;
    return strcmp(((const TxEntry *)a)->tx_id, ((const TxEntry *)b)->tx_id);
}

/*
 * Extract a quoted or unquoted GTF attribute value by key.
 * GTF format: key "value"; or key value;
 * Returns 1 if found, 0 otherwise.
 */
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
                while (*p && *p != ';' && *p != ' ' && *p != '\t' && i < val_sz - 1)
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

int gtf_exon_intron_stats(FILE *fp, GtfExonIntronStats *stats)
{
    memset(stats, 0, sizeof(*stats));

    size_t cap = INIT_CAP;
    TxEntry *rows = malloc(cap * sizeof(TxEntry));
    if (!rows) return -1;
    size_t n_rows = 0;

    char line[65536];
    GtfRecord rec;

    while (fgets(line, sizeof(line), fp)) {
        if (gtf_parse_line(line, &rec) != 0) continue;
        if (strcmp(rec.feature, "exon") != 0) continue;

        if (n_rows >= cap) {
            cap *= 2;
            TxEntry *tmp = realloc(rows, cap * sizeof(TxEntry));
            if (!tmp) { free(rows); return -1; }
            rows = tmp;
        }

        rows[n_rows].exon_count = 1;
        if (!attr_get(rec.attributes, "gene_id",
                      rows[n_rows].gene_id, GTF_FIELD_MAX))
            rows[n_rows].gene_id[0] = '\0';
        if (!attr_get(rec.attributes, "transcript_id",
                      rows[n_rows].tx_id, GTF_FIELD_MAX))
            rows[n_rows].tx_id[0] = '\0';
        n_rows++;
    }

    if (n_rows == 0) { free(rows); return 0; }

    /* --- merge into one TxEntry per transcript --- */
    qsort(rows, n_rows, sizeof(TxEntry), cmp_by_tx);

    TxEntry *txs = malloc(n_rows * sizeof(TxEntry));
    if (!txs) { free(rows); return -1; }
    size_t n_tx = 0;

    for (size_t i = 0; i < n_rows; ) {
        size_t j = i + 1;
        while (j < n_rows &&
               strcmp(rows[j].tx_id, rows[i].tx_id) == 0) j++;
        txs[n_tx]            = rows[i];
        txs[n_tx].exon_count = (long)(j - i);
        n_tx++;
        i = j;
    }
    free(rows);

    /* --- transcript-level stats --- */
    stats->n_transcripts = (long)n_tx;
    for (size_t i = 0; i < n_tx; i++) {
        stats->n_exons   += txs[i].exon_count;
        stats->n_introns += txs[i].exon_count - 1;
        if (txs[i].exon_count == 1)
            stats->n_single_exon_transcripts++;
    }

    /* --- gene-level stats --- */
    qsort(txs, n_tx, sizeof(TxEntry), cmp_by_gene);

    for (size_t i = 0; i < n_tx; ) {
        size_t j = i + 1;
        while (j < n_tx &&
               strcmp(txs[j].gene_id, txs[i].gene_id) == 0) j++;

        stats->n_genes++;

        int all_single = 1;
        for (size_t k = i; k < j; k++)
            if (txs[k].exon_count > 1) { all_single = 0; break; }
        if (all_single) stats->n_single_exon_genes++;

        i = j;
    }
    free(txs);

    if (stats->n_genes > 0) {
        stats->mean_exons_per_gene   = (double)stats->n_exons   / stats->n_genes;
        stats->mean_introns_per_gene = (double)stats->n_introns / stats->n_genes;
        stats->pct_single_exon_genes =
            100.0 * stats->n_single_exon_genes / stats->n_genes;
    }
    if (stats->n_transcripts > 0)
        stats->pct_single_exon_transcripts =
            100.0 * stats->n_single_exon_transcripts / stats->n_transcripts;

    return 0;
}

void gtf_print_exon_intron_stats(const GtfExonIntronStats *stats)
{
    printf("Genes                   : %ld\n",       stats->n_genes);
    printf("Transcripts             : %ld\n",       stats->n_transcripts);
    printf("Exons                   : %ld\n",       stats->n_exons);
    printf("Introns                 : %ld\n",       stats->n_introns);
    printf("Mean exons per gene     : %.2f\n",      stats->mean_exons_per_gene);
    printf("Mean introns per gene   : %.2f\n",      stats->mean_introns_per_gene);
    printf("Single-exon genes       : %ld (%.1f%%)\n",
           stats->n_single_exon_genes, stats->pct_single_exon_genes);
    printf("Single-exon transcripts : %ld (%.1f%%)\n",
           stats->n_single_exon_transcripts, stats->pct_single_exon_transcripts);
}
