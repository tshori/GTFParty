#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gtf.h"

/*
 * GTF has exactly 9 tab-delimited columns.  The 9th (attributes) may itself
 * contain spaces and semicolons, but never a tab — so splitting on '\t'
 * is safe and gives us clean field boundaries.
 *
 * strtok() modifies the string in place, so we copy `line` into a local
 * buffer first.  The copy also lets us strip the trailing newline without
 * touching the caller's memory.
 */
int gtf_parse_line(const char *line, GtfRecord *rec)
{
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\0')
        return 1;

    /* Local copy — stack allocation is fine; GTF lines are at most a few KB */
    char buf[65536];
    strncpy(buf, line, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';

    /* Strip trailing newline */
    size_t len = strlen(buf);
    if (len > 0 && buf[len - 1] == '\n')
        buf[len - 1] = '\0';

    /* Split on tabs, collect up to 9 tokens */
    char *fields[9];
    int n = 0;
    char *tok = strtok(buf, "\t");
    while (tok && n < 9) {
        fields[n++] = tok;
        tok = strtok(NULL, "\t");
    }

    if (n < 9)
        return -1;

    /* --- fixed string fields --- */
    strncpy(rec->seqname,    fields[0], GTF_FIELD_MAX - 1);
    rec->seqname[GTF_FIELD_MAX - 1] = '\0';

    strncpy(rec->source,     fields[1], GTF_FIELD_MAX - 1);
    rec->source[GTF_FIELD_MAX - 1] = '\0';

    strncpy(rec->feature,    fields[2], GTF_FIELD_MAX - 1);
    rec->feature[GTF_FIELD_MAX - 1] = '\0';

    /* --- numeric fields --- */
    rec->start = atol(fields[3]);
    rec->end   = atol(fields[4]);

    /* score: '.' → NAN so callers can test with isnan() */
    rec->score = (fields[5][0] == '.') ? NAN : (float)atof(fields[5]);

    /* strand: single character */
    rec->strand = fields[6][0];

    /* frame: '.' → -1 sentinel */
    rec->frame = (fields[7][0] == '.') ? -1 : atoi(fields[7]);

    /* --- attributes: everything in column 9 --- */
    strncpy(rec->attributes, fields[8], GTF_ATTR_MAX - 1);
    rec->attributes[GTF_ATTR_MAX - 1] = '\0';

    return 0;
}

void gtf_print(const GtfRecord *rec)
{
    printf("%s\t%s\t%s\t%ld\t%ld\t",
           rec->seqname, rec->source, rec->feature,
           rec->start, rec->end);

    if (isnan(rec->score))
        printf(".");
    else
        printf("%g", rec->score);

    printf("\t%c\t", rec->strand);

    if (rec->frame == -1)
        printf(".");
    else
        printf("%d", rec->frame);

    printf("\t%s\n", rec->attributes);
}
