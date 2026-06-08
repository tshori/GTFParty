#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gtf.h"

#define GTF_TABLE_INIT_CAP 4096

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

GtfTable *gtf_table_load(FILE *fp)
{
    GtfTable *t = malloc(sizeof(GtfTable));
    if (!t) return NULL;

    t->cap       = GTF_TABLE_INIT_CAP;
    t->n         = 0;
    t->n_skipped = 0;
    t->n_errors  = 0;
    t->recs      = malloc(t->cap * sizeof(GtfRecord));
    if (!t->recs) { free(t); return NULL; }

    char line[65536];
    GtfRecord rec;

    while (fgets(line, sizeof(line), fp)) {
        int ret = gtf_parse_line(line, &rec);
        if (ret == 1) { t->n_skipped++; continue; }
        if (ret == -1) { t->n_errors++;  continue; }

        if (t->n >= t->cap) {
            size_t new_cap = t->cap * 2;
            GtfRecord *tmp = realloc(t->recs, new_cap * sizeof(GtfRecord));
            if (!tmp) { gtf_table_free(t); return NULL; }
            t->recs = tmp;
            t->cap  = new_cap;
        }
        t->recs[t->n++] = rec;
    }

    return t;
}

void gtf_table_free(GtfTable *t)
{
    if (!t) return;
    free(t->recs);
    free(t);
}

/* ------------------------------------------------------------------ */
/* Attribute parsing                                                    */
/* ------------------------------------------------------------------ */

#define GTF_ATTRS_INIT_CAP 8

/*
 * Duplicate at most `len` bytes of `src` as a heap string.
 * Returns NULL on allocation failure.
 */
static char *strndup_safe(const char *src, size_t len)
{
    char *s = malloc(len + 1);
    if (!s) return NULL;
    memcpy(s, src, len);
    s[len] = '\0';
    return s;
}

GtfAttrs *gtf_attrs_parse(const char *raw)
{
    GtfAttrs *a = malloc(sizeof(GtfAttrs));
    if (!a) return NULL;

    a->n     = 0;
    size_t cap = GTF_ATTRS_INIT_CAP;
    a->pairs = malloc(cap * sizeof(GtfAttr));
    if (!a->pairs) { free(a); return NULL; }

    const char *p = raw;
    while (*p) {
        /* skip leading whitespace before key */
        while (*p == ' ' || *p == '\t') p++;
        if (!*p) break;

        /* extract key: runs to first whitespace or ';' */
        const char *key_start = p;
        while (*p && *p != ' ' && *p != '\t' && *p != ';') p++;
        size_t key_len = (size_t)(p - key_start);
        if (key_len == 0) {
            while (*p && *p != ';') p++;
            if (*p == ';') p++;
            continue;
        }

        char *key = strndup_safe(key_start, key_len);
        if (!key) { gtf_attrs_free(a); return NULL; }

        /* skip whitespace before value */
        while (*p == ' ' || *p == '\t') p++;

        /* extract value: quoted or unquoted */
        char *value = NULL;
        if (*p == '"') {
            p++;  /* skip opening quote */
            const char *val_start = p;
            while (*p && *p != '"') p++;
            value = strndup_safe(val_start, (size_t)(p - val_start));
            if (*p == '"') p++;  /* skip closing quote */
        } else {
            const char *val_start = p;
            while (*p && *p != ';' && *p != ' ' && *p != '\t') p++;
            value = strndup_safe(val_start, (size_t)(p - val_start));
        }
        if (!value) { free(key); gtf_attrs_free(a); return NULL; }

        /* grow pairs array if needed */
        if (a->n >= cap) {
            cap *= 2;
            GtfAttr *tmp = realloc(a->pairs, cap * sizeof(GtfAttr));
            if (!tmp) { free(key); free(value); gtf_attrs_free(a); return NULL; }
            a->pairs = tmp;
        }
        a->pairs[a->n].key   = key;
        a->pairs[a->n].value = value;
        a->n++;

        /* advance past the closing semicolon */
        while (*p && *p != ';') p++;
        if (*p == ';') p++;
    }

    return a;
}

const char *gtf_attrs_get(const GtfAttrs *a, const char *key)
{
    for (size_t i = 0; i < a->n; i++)
        if (strcmp(a->pairs[i].key, key) == 0)
            return a->pairs[i].value;
    return NULL;
}

void gtf_attrs_free(GtfAttrs *a)
{
    if (!a) return;
    for (size_t i = 0; i < a->n; i++) {
        free(a->pairs[i].key);
        free(a->pairs[i].value);
    }
    free(a->pairs);
    free(a);
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
