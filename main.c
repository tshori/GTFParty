#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gtf.h"
#include "stats.h"
#include "compare.h"

static void usage(const char *prog)
{
    fprintf(stderr, "Usage: %s <command> [args]\n\n", prog);
    fprintf(stderr, "Commands:\n");
    fprintf(stderr, "  print   <file>                   print all records to stdout (GTF or GFF3)\n");
    fprintf(stderr, "  filter  <feature> <file>         print records matching feature type\n");
    fprintf(stderr, "  attrs   <key> <file>             extract one attribute value per record\n");
    fprintf(stderr, "  overlap <seqname:start-end> <file>  records overlapping a region\n");
    fprintf(stderr, "  stats   <file.gtf>               exon/intron statistics (GTF only)\n");
    fprintf(stderr, "  compare <file1> <file2>          structural comparison (GTF only)\n");
}

static int cmd_print(FILE *fp)
{
    GtfTable *t = gtf_table_load(fp);
    if (!t) {
        fprintf(stderr, "error: out of memory\n");
        return 1;
    }

    for (size_t i = 0; i < t->n; i++)
        gtf_print(&t->recs[i]);

    fprintf(stderr, "parsed=%zu  skipped=%zu  errors=%zu  format=%s\n",
            t->n, t->n_skipped, t->n_errors,
            t->format == GTF_FMT_GFF3 ? "gff3" : "gtf");

    int rc = t->n_errors > 0 ? 1 : 0;
    gtf_table_free(t);
    return rc;
}

static int cmd_attrs(FILE *fp, const char *key)
{
    GtfTable *t = gtf_table_load(fp);
    if (!t) {
        fprintf(stderr, "error: out of memory\n");
        return 1;
    }

    for (size_t i = 0; i < t->n; i++) {
        GtfAttrs *a = (t->format == GTF_FMT_GFF3)
                      ? gff3_attrs_parse(t->recs[i].attributes)
                      : gtf_attrs_parse(t->recs[i].attributes);
        if (!a) {
            fprintf(stderr, "error: out of memory\n");
            gtf_table_free(t);
            return 1;
        }
        const char *val = gtf_attrs_get(a, key);
        printf("%s\n", val ? val : ".");
        gtf_attrs_free(a);
    }

    gtf_table_free(t);
    return 0;
}

static int cmd_overlap(FILE *fp, const char *region)
{
    /* Parse region string: seqname:start-end  (1-based, fully closed) */
    const char *colon = strrchr(region, ':');
    if (!colon) {
        fprintf(stderr, "error: invalid region '%s' — expected seqname:start-end\n",
                region);
        return 1;
    }
    size_t seqlen = (size_t)(colon - region);
    if (seqlen == 0 || seqlen >= GTF_FIELD_MAX) {
        fprintf(stderr, "error: seqname in region '%s' is empty or too long\n", region);
        return 1;
    }
    char seqname[GTF_FIELD_MAX];
    memcpy(seqname, region, seqlen);
    seqname[seqlen] = '\0';

    long qstart, qend;
    if (sscanf(colon + 1, "%ld-%ld", &qstart, &qend) != 2 || qstart > qend) {
        fprintf(stderr, "error: invalid coordinates in region '%s'\n", region);
        return 1;
    }

    GtfTable *t = gtf_table_load(fp);
    if (!t) { fprintf(stderr, "error: out of memory\n"); return 1; }

    gtf_table_sort(t);

    size_t count;
    GtfRecord **hits = gtf_table_query(t, seqname, qstart, qend, &count);

    for (size_t i = 0; i < count; i++)
        gtf_print(hits[i]);

    fprintf(stderr, "hits=%zu  region=%s:%ld-%ld\n", count, seqname, qstart, qend);

    free(hits);
    gtf_table_free(t);
    return 0;
}

static int cmd_filter(FILE *fp, const char *feature)
{
    GtfTable *t = gtf_table_load(fp);
    if (!t) {
        fprintf(stderr, "error: out of memory\n");
        return 1;
    }

    size_t matched = 0;
    for (size_t i = 0; i < t->n; i++) {
        if (strcmp(t->recs[i].feature, feature) == 0) {
            gtf_print(&t->recs[i]);
            matched++;
        }
    }

    fprintf(stderr, "matched=%zu  total=%zu\n", matched, t->n);
    gtf_table_free(t);
    return 0;
}

static int cmd_stats(FILE *fp)
{
    GtfExonIntronStats stats;
    if (gtf_exon_intron_stats(fp, &stats) != 0) {
        fprintf(stderr, "error: stats computation failed (out of memory?)\n");
        return 1;
    }
    gtf_print_exon_intron_stats(&stats);
    return 0;
}

static int cmd_compare(const char *path_a, const char *path_b)
{
    FILE *fp_a = fopen(path_a, "r");
    if (!fp_a) { perror(path_a); return 1; }
    FILE *fp_b = fopen(path_b, "r");
    if (!fp_b) { perror(path_b); fclose(fp_a); return 1; }

    fprintf(stderr, "Loading %s ...\n", path_a);
    size_t na;
    Transcript *a = cmp_build_transcripts(fp_a, 1, &na);
    fclose(fp_a);

    fprintf(stderr, "Loading %s ...\n", path_b);
    size_t nb;
    Transcript *b = cmp_build_transcripts(fp_b, 2, &nb);
    fclose(fp_b);

    if (!a && na == 0) fprintf(stderr, "warning: no exon features in %s\n", path_a);
    if (!b && nb == 0) fprintf(stderr, "warning: no exon features in %s\n", path_b);

    fprintf(stderr, "Loaded %zu + %zu transcripts.  Comparing ...\n\n", na, nb);

    CmpStats stats;
    int rc = cmp_compare(a, na, b, nb, &stats);
    if (rc != 0) {
        fprintf(stderr, "error: comparison failed (out of memory?)\n");
        cmp_free_transcripts(a, na);
        cmp_free_transcripts(b, nb);
        return 1;
    }

    cmp_print_stats(&stats, path_a, path_b);

    cmp_free_transcripts(a, na);
    cmp_free_transcripts(b, nb);
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    const char *cmd = argv[1];

    /* compare takes two file arguments */
    if (strcmp(cmd, "compare") == 0) {
        if (argc < 4) {
            fprintf(stderr, "Usage: %s compare <file1.gtf> <file2.gtf>\n",
                    argv[0]);
            return 1;
        }
        return cmd_compare(argv[2], argv[3]);
    }

    /* filter, attrs, and overlap all take one string arg then a file */
    if (strcmp(cmd, "filter")  == 0 ||
        strcmp(cmd, "attrs")   == 0 ||
        strcmp(cmd, "overlap") == 0) {
        if (argc < 4) {
            fprintf(stderr, "Usage: %s %s <%s> <file>\n", argv[0], cmd,
                    strcmp(cmd, "overlap") == 0 ? "seqname:start-end"
                  : strcmp(cmd, "filter")  == 0 ? "feature" : "key");
            return 1;
        }
        FILE *fp = fopen(argv[3], "r");
        if (!fp) { perror(argv[3]); return 1; }
        int rc;
        if      (strcmp(cmd, "filter")  == 0) rc = cmd_filter(fp,  argv[2]);
        else if (strcmp(cmd, "attrs")   == 0) rc = cmd_attrs(fp,   argv[2]);
        else                                  rc = cmd_overlap(fp,  argv[2]);
        fclose(fp);
        return rc;
    }

    /* all other commands take one file argument */
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    FILE *fp = fopen(argv[2], "r");
    if (!fp) {
        perror(argv[2]);
        return 1;
    }

    int ret;
    if (strcmp(cmd, "print") == 0)
        ret = cmd_print(fp);
    else if (strcmp(cmd, "stats") == 0)
        ret = cmd_stats(fp);
    else {
        fprintf(stderr, "unknown command: %s\n\n", cmd);
        usage(argv[0]);
        ret = 1;
    }

    fclose(fp);
    return ret;
}
