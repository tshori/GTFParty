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
    fprintf(stderr, "  print   <file.gtf>              print all records to stdout\n");
    fprintf(stderr, "  stats   <file.gtf>              exon/intron statistics\n");
    fprintf(stderr, "  compare <file1.gtf> <file2.gtf> structural comparison by coordinate\n");
}

static int cmd_print(FILE *fp)
{
    char line[65536];
    GtfRecord rec;
    int parsed = 0, skipped = 0, errors = 0;

    while (fgets(line, sizeof(line), fp)) {
        int ret = gtf_parse_line(line, &rec);
        if (ret == 1) { skipped++; continue; }
        if (ret == -1) { errors++; continue; }
        gtf_print(&rec);
        parsed++;
    }
    fprintf(stderr, "parsed=%d  skipped=%d  errors=%d\n",
            parsed, skipped, errors);
    return errors > 0 ? 1 : 0;
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
