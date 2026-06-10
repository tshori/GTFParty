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
    t->format    = GTF_FMT_GTF;
    t->recs      = malloc(t->cap * sizeof(GtfRecord));
    if (!t->recs) { free(t); return NULL; }

    char line[65536];
    GtfRecord rec;

    while (fgets(line, sizeof(line), fp)) {
        /* GFF3: ##FASTA marks the start of embedded sequences — stop here */
        if (strncmp(line, "##FASTA", 7) == 0)
            break;
        /* GFF3: ##gff-version pragma identifies the format */
        if (strncmp(line, "##gff-version", 13) == 0) {
            t->format = GTF_FMT_GFF3;
            t->n_skipped++;
            continue;
        }

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
/* Sort and overlap query                                               */
/* ------------------------------------------------------------------ */

static int rec_cmp(const void *a, const void *b)
{
    const GtfRecord *ra = (const GtfRecord *)a;
    const GtfRecord *rb = (const GtfRecord *)b;
    int c = strcmp(ra->seqname, rb->seqname);
    if (c) return c;
    if (ra->start != rb->start) return (ra->start < rb->start) ? -1 : 1;
    return (ra->end < rb->end) ? -1 : (ra->end > rb->end) ? 1 : 0;
}

void gtf_table_sort(GtfTable *t)
{
    qsort(t->recs, t->n, sizeof(GtfRecord), rec_cmp);
}

/*
 * Return the index of the first record whose seqname is >= target.
 * This is the standard lower-bound binary search adapted for strings.
 */
static size_t lower_bound_seqname(const GtfTable *t, const char *target)
{
    size_t lo = 0, hi = t->n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (strcmp(t->recs[mid].seqname, target) < 0)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

GtfRecord **gtf_table_query(const GtfTable *t, const char *seqname,
                             long qstart, long qend, size_t *count)
{
    *count = 0;

    /* Jump to the first record whose seqname >= target */
    size_t lo = lower_bound_seqname(t, seqname);

    size_t cap = 16;
    GtfRecord **hits = malloc(cap * sizeof(GtfRecord *));
    if (!hits) return NULL;

    for (size_t i = lo; i < t->n; i++) {
        int cmp = strcmp(t->recs[i].seqname, seqname);
        if (cmp > 0) break;           /* past target seqname — done */
        if (t->recs[i].start > qend) break; /* sorted by start: no more overlaps */
        if (t->recs[i].end < qstart)  continue; /* starts in range but ends before query */

        if (*count >= cap) {
            cap *= 2;
            GtfRecord **tmp = realloc(hits, cap * sizeof(GtfRecord *));
            if (!tmp) { free(hits); *count = 0; return NULL; }
            hits = tmp;
        }
        hits[(*count)++] = &t->recs[i];
    }

    if (*count == 0) { free(hits); return NULL; }
    return hits;
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

/*
 * GFF3 attributes column: key=value;key=value;...
 * Keys and values are unquoted.  Values may be comma-separated lists;
 * the whole list is stored as a single string.
 * Percent-encoded characters (%XX) are left encoded — see README for details.
 */
GtfAttrs *gff3_attrs_parse(const char *raw)
{
    GtfAttrs *a = malloc(sizeof(GtfAttrs));
    if (!a) return NULL;

    a->n = 0;
    size_t cap = GTF_ATTRS_INIT_CAP;
    a->pairs = malloc(cap * sizeof(GtfAttr));
    if (!a->pairs) { free(a); return NULL; }

    const char *p = raw;
    while (*p) {
        while (*p == ' ' || *p == '\t') p++;
        if (!*p) break;

        /* locate the end of this key=value token (next ';' or end) */
        const char *tok_end = p;
        while (*tok_end && *tok_end != ';') tok_end++;

        /* locate the '=' separator within the token */
        const char *eq = p;
        while (eq < tok_end && *eq != '=') eq++;

        if (eq == tok_end) {
            /* no '=' — bare tag with no value; store with empty value */
            char *key = strndup_safe(p, (size_t)(tok_end - p));
            if (!key) { gtf_attrs_free(a); return NULL; }
            char *value = strndup_safe("", 0);
            if (!value) { free(key); gtf_attrs_free(a); return NULL; }

            if (a->n >= cap) {
                cap *= 2;
                GtfAttr *tmp = realloc(a->pairs, cap * sizeof(GtfAttr));
                if (!tmp) { free(key); free(value); gtf_attrs_free(a); return NULL; }
                a->pairs = tmp;
            }
            a->pairs[a->n].key   = key;
            a->pairs[a->n].value = value;
            a->n++;
        } else {
            char *key   = strndup_safe(p, (size_t)(eq - p));
            if (!key) { gtf_attrs_free(a); return NULL; }
            char *value = strndup_safe(eq + 1, (size_t)(tok_end - eq - 1));
            if (!value) { free(key); gtf_attrs_free(a); return NULL; }

            if (a->n >= cap) {
                cap *= 2;
                GtfAttr *tmp = realloc(a->pairs, cap * sizeof(GtfAttr));
                if (!tmp) { free(key); free(value); gtf_attrs_free(a); return NULL; }
                a->pairs = tmp;
            }
            a->pairs[a->n].key   = key;
            a->pairs[a->n].value = value;
            a->n++;
        }

        p = tok_end;
        if (*p == ';') p++;
    }

    return a;
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
