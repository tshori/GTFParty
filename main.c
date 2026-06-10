#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gtf.h"
#include "stats.h"
#include "compare.h"

static void usage(const char *prog)
{
    fprintf(stderr, "Usage: %s <command> [args] [-o <file>]\n\n", prog);
    fprintf(stderr, "  -o, --output <file>  write output to <file> instead of stdout\n\n");
    fprintf(stderr, "Commands:\n");
    fprintf(stderr, "  print   [--chrom chr,...] [--source src,...] <file>\n");
    fprintf(stderr, "                                   print records, optionally filtered by chromosome or source\n");
    fprintf(stderr, "  filter  <feature> <file>         print records matching feature type\n");
    fprintf(stderr, "  merge   <feature> <file>         merge overlapping/adjacent intervals of a feature type\n");
    fprintf(stderr, "  attrs   <key> <file>             extract one attribute value per record\n");
    fprintf(stderr, "  overlap <seqname:start-end> <file>  records overlapping a region\n");
    fprintf(stderr, "  bed     [--name <attr>] <file>   BED6 output (0-based half-open coords)\n");
    fprintf(stderr, "  stats   <file.gtf>               exon/intron statistics (GTF only)\n");
    fprintf(stderr, "  compare <file1> <file2>          structural comparison (GTF only)\n");
}

/* Returns 1 if value appears in a comma-separated list, or if list is NULL. */
static int in_csv_list(const char *list, const char *value)
{
    if (!list) return 1;
    size_t vlen = strlen(value);
    const char *p = list;
    while (*p) {
        const char *end = strchr(p, ',');
        size_t len = end ? (size_t)(end - p) : strlen(p);
        if (len == vlen && strncmp(p, value, len) == 0) return 1;
        if (!end) break;
        p = end + 1;
    }
    return 0;
}

static int cmd_print(FILE *fp, const char *chrom_list, const char *source_list)
{
    GtfTable *t = gtf_table_load(fp);
    if (!t) {
        fprintf(stderr, "error: out of memory\n");
        return 1;
    }

    size_t printed = 0;
    for (size_t i = 0; i < t->n; i++) {
        if (!in_csv_list(chrom_list,  t->recs[i].seqname)) continue;
        if (!in_csv_list(source_list, t->recs[i].source))  continue;
        gtf_print(&t->recs[i]);
        printed++;
    }

    fprintf(stderr, "printed=%zu  parsed=%zu  skipped=%zu  errors=%zu  format=%s\n",
            printed, t->n, t->n_skipped, t->n_errors,
            t->format == GTF_FMT_GFF3 ? "gff3" : "gtf");

    int rc = t->n_errors > 0 ? 1 : 0;
    gtf_table_free(t);
    return rc;
}

static int cmd_bed(FILE *fp, const char *name_attr)
{
    GtfTable *t = gtf_table_load(fp);
    if (!t) {
        fprintf(stderr, "error: out of memory\n");
        return 1;
    }

    for (size_t i = 0; i < t->n; i++) {
        const GtfRecord *r = &t->recs[i];

        /* Score: NaN → 0; otherwise clamp to BED range [0, 1000] */
        int score = 0;
        if (!isnan(r->score)) {
            float s = r->score;
            score = (s < 0.0f) ? 0 : (s > 1000.0f) ? 1000 : (int)s;
        }

        /* Name from requested attribute, or "." if absent / not requested */
        const char *name = ".";
        GtfAttrs *attrs = NULL;
        if (name_attr) {
            attrs = (t->format == GTF_FMT_GFF3)
                    ? gff3_attrs_parse(r->attributes)
                    : gtf_attrs_parse(r->attributes);
            if (attrs) {
                const char *val = gtf_attrs_get(attrs, name_attr);
                if (val) name = val;
            }
        }

        /* GTF: 1-based closed [start, end]  →  BED: 0-based half-open */
        printf("%s\t%ld\t%ld\t%s\t%d\t%c\n",
               r->seqname,
               r->start - 1,
               r->end,
               name,
               score,
               r->strand);

        if (attrs) gtf_attrs_free(attrs);
    }

    fprintf(stderr, "records=%zu  format=%s\n", t->n,
            t->format == GTF_FMT_GFF3 ? "gff3" : "gtf");
    gtf_table_free(t);
    return 0;
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

    fprintf(stderr, "  sorting %zu records...\n", t->n);
    gtf_table_sort(t);
    fprintf(stderr, "  querying %s:%ld-%ld...\n", seqname, qstart, qend);

    size_t count;
    GtfRecord **hits = gtf_table_query(t, seqname, qstart, qend, &count);

    for (size_t i = 0; i < count; i++)
        gtf_print(hits[i]);

    fprintf(stderr, "hits=%zu  region=%s:%ld-%ld\n", count, seqname, qstart, qend);

    free(hits);
    gtf_table_free(t);
    return 0;
}

static int cmd_merge(FILE *fp, const char *feature)
{
    GtfTable *t = gtf_table_load(fp);
    if (!t) { fprintf(stderr, "error: out of memory\n"); return 1; }

    fprintf(stderr, "  sorting %zu records...\n", t->n);
    gtf_table_sort(t);

    size_t total = 0, n_out = 0;
    int    in_merge = 0;
    GtfRecord cur;
    size_t cur_count = 0;

    /*
     * merge_attrs lives for the whole function.
     * We temporarily point cur.attributes here before printing so that
     * gtf_print() emits the merged count instead of the original attributes.
     * No heap allocation needed: the pointer swap is safe because the buffer
     * outlives every gtf_print() call below.
     */
    char merge_attrs[64];

    for (size_t i = 0; i < t->n; i++) {
        GtfRecord *r = &t->recs[i];
        if (strcmp(r->feature, feature) != 0) continue;
        total++;

        if (!in_merge) {
            cur = *r;
            cur_count = 1;
            in_merge  = 1;
            continue;
        }

        /* same chromosome and overlapping or adjacent (1-based closed) */
        if (strcmp(r->seqname, cur.seqname) == 0 && r->start <= cur.end + 1) {
            if (r->end > cur.end) cur.end = r->end;
            cur_count++;
        } else {
            snprintf(merge_attrs, sizeof(merge_attrs),
                     "merged_count \"%zu\";", cur_count);
            cur.attributes = merge_attrs;
            gtf_print(&cur);
            n_out++;
            cur = *r;
            cur_count = 1;
        }
    }

    if (in_merge) {
        snprintf(merge_attrs, sizeof(merge_attrs),
                 "merged_count \"%zu\";", cur_count);
        cur.attributes = merge_attrs;
        gtf_print(&cur);
        n_out++;
    }

    fprintf(stderr, "input=%zu  merged_output=%zu\n", total, n_out);
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

    /*
     * Pre-pass: scan argv for --output / -o and compact it out.
     * We shift remaining elements left in-place and decrement argc so that
     * all downstream dispatch code sees a clean argv with no output flag.
     * freopen() replaces stdout's underlying file descriptor; every subsequent
     * printf / gtf_print call then writes to the requested file automatically.
     */
    const char *out_path = NULL;
    for (int i = 1; i < argc; ) {
        if ((strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0)
                && i + 1 < argc) {
            out_path = argv[i + 1];
            for (int j = i; j < argc - 2; j++)
                argv[j] = argv[j + 2];
            argc -= 2;
        } else if (strncmp(argv[i], "--output=", 9) == 0) {
            out_path = argv[i] + 9;
            for (int j = i; j < argc - 1; j++)
                argv[j] = argv[j + 1];
            argc--;
        } else {
            i++;
        }
    }

    if (out_path) {
        if (!freopen(out_path, "w", stdout)) {
            perror(out_path);
            return 1;
        }
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

    /* filter, merge, attrs, and overlap all take one string arg then a file */
    if (strcmp(cmd, "filter")  == 0 ||
        strcmp(cmd, "merge")   == 0 ||
        strcmp(cmd, "attrs")   == 0 ||
        strcmp(cmd, "overlap") == 0) {
        if (argc < 4) {
            fprintf(stderr, "Usage: %s %s <%s> <file>\n", argv[0], cmd,
                    strcmp(cmd, "overlap") == 0 ? "seqname:start-end"
                  : strcmp(cmd, "merge")   == 0 ? "feature"
                  : strcmp(cmd, "filter")  == 0 ? "feature" : "key");
            return 1;
        }
        FILE *fp = fopen(argv[3], "r");
        if (!fp) { perror(argv[3]); return 1; }
        int rc;
        if      (strcmp(cmd, "filter")  == 0) rc = cmd_filter(fp,  argv[2]);
        else if (strcmp(cmd, "merge")   == 0) rc = cmd_merge(fp,   argv[2]);
        else if (strcmp(cmd, "attrs")   == 0) rc = cmd_attrs(fp,   argv[2]);
        else                                  rc = cmd_overlap(fp,  argv[2]);
        fclose(fp);
        return rc;
    }

    /* print accepts optional --chrom and --source flags before the file */
    if (strcmp(cmd, "print") == 0) {
        const char *path        = NULL;
        const char *chrom_list  = NULL;
        const char *source_list = NULL;

        for (int i = 2; i < argc; i++) {
            if (strcmp(argv[i], "--chrom") == 0 && i + 1 < argc) {
                chrom_list = argv[++i];
            } else if (strncmp(argv[i], "--chrom=", 8) == 0) {
                chrom_list = argv[i] + 8;
            } else if (strcmp(argv[i], "--source") == 0 && i + 1 < argc) {
                source_list = argv[++i];
            } else if (strncmp(argv[i], "--source=", 9) == 0) {
                source_list = argv[i] + 9;
            } else if (argv[i][0] != '-') {
                path = argv[i];
            } else {
                fprintf(stderr, "error: unknown option '%s'\n", argv[i]);
                return 1;
            }
        }

        if (!path) {
            fprintf(stderr, "Usage: %s print [--chrom chr,...] [--source src,...] <file>\n",
                    argv[0]);
            return 1;
        }

        FILE *fp = fopen(path, "r");
        if (!fp) { perror(path); return 1; }
        int rc = cmd_print(fp, chrom_list, source_list);
        fclose(fp);
        return rc;
    }

    /* bed accepts optional --name flag before the file */
    if (strcmp(cmd, "bed") == 0) {
        const char *path      = NULL;
        const char *name_attr = NULL;

        for (int i = 2; i < argc; i++) {
            if (strcmp(argv[i], "--name") == 0 && i + 1 < argc) {
                name_attr = argv[++i];
            } else if (strncmp(argv[i], "--name=", 7) == 0) {
                name_attr = argv[i] + 7;
            } else if (argv[i][0] != '-') {
                path = argv[i];
            } else {
                fprintf(stderr, "error: unknown option '%s'\n", argv[i]);
                return 1;
            }
        }

        if (!path) {
            fprintf(stderr, "Usage: %s bed [--name <attr>] <file>\n", argv[0]);
            return 1;
        }

        FILE *fp = fopen(path, "r");
        if (!fp) { perror(path); return 1; }
        int rc = cmd_bed(fp, name_attr);
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
    if (strcmp(cmd, "stats") == 0)
        ret = cmd_stats(fp);
    else {
        fprintf(stderr, "unknown command: %s\n\n", cmd);
        usage(argv[0]);
        ret = 1;
    }

    fclose(fp);
    return ret;
}
