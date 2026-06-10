#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gtf.h"

#define GTF_TABLE_INIT_CAP 4096

/* Duplicate at most `len` bytes of `src` as a heap string. */
static char *strndup_safe(const char *src, size_t len)
{
    char *s = malloc(len + 1);
    if (!s) return NULL;
    memcpy(s, src, len);
    s[len] = '\0';
    return s;
}

/* ------------------------------------------------------------------ */
/* String intern pool                                                    */
/* ------------------------------------------------------------------ */

/*
 * Return a pointer to an interned copy of s.
 * If s is already in the pool the existing pointer is returned — same string,
 * same address.  Otherwise a new heap copy is added and its address returned.
 *
 * Linear scan is fine here: there are at most a few hundred unique
 * seqname/source/feature strings even in large files.
 */
static const char *str_pool_intern(StrPool *pool, const char *s)
{
    for (size_t i = 0; i < pool->n; i++)
        if (strcmp(pool->data[i], s) == 0) return pool->data[i];

    if (pool->n >= pool->cap) {
        size_t new_cap = pool->cap ? pool->cap * 2 : 16;
        char **tmp = realloc(pool->data, new_cap * sizeof(char *));
        if (!tmp) return NULL;
        pool->data = tmp;
        pool->cap  = new_cap;
    }
    char *dup = strndup_safe(s, strlen(s));
    if (!dup) return NULL;
    pool->data[pool->n++] = dup;
    return dup;
}

/* ------------------------------------------------------------------ */
/* Open-addressing hash table  (const char * key → size_t value)       */
/* ------------------------------------------------------------------ */

/*
 * FNV-1a 64-bit hash — fast, good avalanche, no dependencies.
 * https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
 */
static size_t fnv1a(const char *s)
{
    size_t h = 14695981039346656037ULL;
    for (; *s; s++) {
        h ^= (unsigned char)*s;
        h *= 1099511628211ULL;
    }
    return h;
}

static HTable *htable_new(size_t init_cap)
{
    HTable *ht = malloc(sizeof(HTable));
    if (!ht) return NULL;
    ht->slots = calloc(init_cap, sizeof(HEntry));
    if (!ht->slots) { free(ht); return NULL; }
    ht->cap = init_cap;
    ht->n   = 0;
    return ht;
}

void htable_free(HTable *ht)
{
    if (!ht) return;
    free(ht->slots);
    free(ht);
}

/* Rebuild into a larger slot array; called when load > 0.75. */
static int htable_resize(HTable *ht, size_t new_cap)
{
    HEntry *new_slots = calloc(new_cap, sizeof(HEntry));
    if (!new_slots) return -1;
    for (size_t i = 0; i < ht->cap; i++) {
        if (!ht->slots[i].key) continue;
        size_t h = fnv1a(ht->slots[i].key) & (new_cap - 1);
        while (new_slots[h].key)
            h = (h + 1) & (new_cap - 1);
        new_slots[h] = ht->slots[i];
    }
    free(ht->slots);
    ht->slots = new_slots;
    ht->cap   = new_cap;
    return 0;
}

static int htable_insert(HTable *ht, const char *key, size_t val)
{
    if (ht->n * 4 > ht->cap * 3) {
        if (htable_resize(ht, ht->cap * 2) < 0) return -1;
    }
    size_t h = fnv1a(key) & (ht->cap - 1);
    while (ht->slots[h].key) {
        if (ht->slots[h].key == key) { ht->slots[h].val = val; return 0; }
        h = (h + 1) & (ht->cap - 1);
    }
    ht->slots[h].key = key;
    ht->slots[h].val = val;
    ht->n++;
    return 0;
}

/*
 * Look up key by string content (strcmp).
 * The stored keys are interned pointers, but the caller's key may be any
 * string — a stack buffer from cmd_overlap, for instance.
 * After the slot is found, the interned pointer (slot->key) is used for the
 * fast pointer-equality scan inside gtf_table_query.
 */
static size_t *htable_get(const HTable *ht, const char *key)
{
    size_t h = fnv1a(key) & (ht->cap - 1);
    while (ht->slots[h].key) {
        if (strcmp(ht->slots[h].key, key) == 0) return &ht->slots[h].val;
        h = (h + 1) & (ht->cap - 1);
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/* GTF/GFF3 line parser                                                 */
/* ------------------------------------------------------------------ */

/*
 * Parse one GTF/GFF3 line into a GtfRawRecord scratch buffer.
 * The raw record uses fixed-size arrays — it is meant to live on the stack
 * for the duration of one loop iteration, never stored long-term.
 *
 * Returns 0 on success, 1 to skip (comment/blank), -1 on parse error.
 */
int gtf_parse_line(const char *line, GtfRawRecord *rec)
{
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\0')
        return 1;

    char buf[65536];
    strncpy(buf, line, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';

    size_t len = strlen(buf);
    if (len > 0 && buf[len - 1] == '\n')
        buf[len - 1] = '\0';

    char *fields[9];
    int n = 0;
    char *tok = strtok(buf, "\t");
    while (tok && n < 9) {
        fields[n++] = tok;
        tok = strtok(NULL, "\t");
    }
    if (n < 9) return -1;

    strncpy(rec->seqname,    fields[0], GTF_FIELD_MAX - 1);
    rec->seqname[GTF_FIELD_MAX - 1] = '\0';
    strncpy(rec->source,     fields[1], GTF_FIELD_MAX - 1);
    rec->source[GTF_FIELD_MAX - 1]  = '\0';
    strncpy(rec->feature,    fields[2], GTF_FIELD_MAX - 1);
    rec->feature[GTF_FIELD_MAX - 1] = '\0';

    rec->start  = atol(fields[3]);
    rec->end    = atol(fields[4]);
    rec->score  = (fields[5][0] == '.') ? NAN : (float)atof(fields[5]);
    rec->strand = fields[6][0];
    rec->frame  = (fields[7][0] == '.') ? -1 : atoi(fields[7]);

    strncpy(rec->attributes, fields[8], GTF_ATTR_MAX - 1);
    rec->attributes[GTF_ATTR_MAX - 1] = '\0';

    return 0;
}

/* ------------------------------------------------------------------ */
/* GtfTable: load, free, sort, query                                    */
/* ------------------------------------------------------------------ */

/*
 * Load all records from fp.
 *
 * Each unique seqname/source/feature string is interned in t->pool so that
 * records sharing the same chromosome name share exactly one heap allocation
 * — and their pointers compare equal, enabling O(1) string equality checks.
 *
 * attributes is heap-allocated per record via strdup (it is always unique).
 */
GtfTable *gtf_table_load(FILE *fp)
{
    GtfTable *t = malloc(sizeof(GtfTable));
    if (!t) return NULL;

    t->cap       = GTF_TABLE_INIT_CAP;
    t->n         = 0;
    t->n_skipped = 0;
    t->n_errors  = 0;
    t->format    = GTF_FMT_GTF;
    t->index     = NULL;
    t->pool.data = NULL;
    t->pool.n    = 0;
    t->pool.cap  = 0;

    t->recs = malloc(t->cap * sizeof(GtfRecord));
    if (!t->recs) { free(t); return NULL; }

    char line[65536];
    GtfRawRecord raw;

    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "##FASTA", 7) == 0) break;
        if (strncmp(line, "##gff-version", 13) == 0) {
            t->format = GTF_FMT_GFF3;
            t->n_skipped++;
            continue;
        }

        int ret = gtf_parse_line(line, &raw);
        if (ret == 1) { t->n_skipped++; continue; }
        if (ret == -1) { t->n_errors++;  continue; }

        if (t->n >= t->cap) {
            size_t new_cap = t->cap * 2;
            GtfRecord *tmp = realloc(t->recs, new_cap * sizeof(GtfRecord));
            if (!tmp) { gtf_table_free(t); return NULL; }
            t->recs = tmp;
            t->cap  = new_cap;
        }

        GtfRecord *r = &t->recs[t->n];
        r->seqname    = str_pool_intern(&t->pool, raw.seqname);
        r->source     = str_pool_intern(&t->pool, raw.source);
        r->feature    = str_pool_intern(&t->pool, raw.feature);
        r->start      = raw.start;
        r->end        = raw.end;
        r->score      = raw.score;
        r->strand     = raw.strand;
        r->frame      = raw.frame;
        r->attributes = strndup_safe(raw.attributes, strlen(raw.attributes));

        if (!r->seqname || !r->source || !r->feature || !r->attributes) {
            gtf_table_free(t);
            return NULL;
        }

        t->n++;
        if (t->n % 100000 == 0)
            fprintf(stderr, "  loading: %zu records...\n", t->n);
    }

    return t;
}

void gtf_table_free(GtfTable *t)
{
    if (!t) return;
    for (size_t i = 0; i < t->n; i++)
        free(t->recs[i].attributes);
    free(t->recs);
    for (size_t i = 0; i < t->pool.n; i++)
        free(t->pool.data[i]);
    free(t->pool.data);
    htable_free(t->index);
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

/*
 * Sort records by (seqname, start, end), then build the seqname index.
 *
 * The index maps each unique seqname to the index of its first record in the
 * sorted array.  Because seqnames are interned, block boundaries are detected
 * with pointer comparison (O(1)) rather than strcmp (O(n)).
 */
void gtf_table_sort(GtfTable *t)
{
    qsort(t->recs, t->n, sizeof(GtfRecord), rec_cmp);

    htable_free(t->index);
    t->index = htable_new(64);
    if (!t->index) return;

    for (size_t i = 0; i < t->n; i++) {
        /* Pointer inequality is sufficient because seqnames are interned */
        if (i == 0 || t->recs[i].seqname != t->recs[i - 1].seqname)
            htable_insert(t->index, t->recs[i].seqname, i);
    }
}

/*
 * Fallback used when the index is not available (table not yet sorted).
 * Standard lower-bound binary search.
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

    const char *interned = NULL;
    size_t lo;

    if (t->index) {
        /*
         * Hash lookup gives us the first index and the interned pointer.
         * The interned pointer lets us use pointer comparison (not strcmp)
         * to detect the end of the seqname block inside the scan loop.
         */
        size_t *idx = htable_get(t->index, seqname);
        if (!idx) return NULL;
        lo       = *idx;
        interned = t->recs[lo].seqname;
    } else {
        lo = lower_bound_seqname(t, seqname);
    }

    size_t cap = 16;
    GtfRecord **hits = malloc(cap * sizeof(GtfRecord *));
    if (!hits) return NULL;

    for (size_t i = lo; i < t->n; i++) {
        if (interned) {
            if (t->recs[i].seqname != interned) break;
        } else {
            int cmp = strcmp(t->recs[i].seqname, seqname);
            if (cmp > 0) break;
        }
        if (t->recs[i].start > qend) break;
        if (t->recs[i].end   < qstart) continue;

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
        while (*p == ' ' || *p == '\t') p++;
        if (!*p) break;

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

        while (*p == ' ' || *p == '\t') p++;

        char *value = NULL;
        if (*p == '"') {
            p++;
            const char *val_start = p;
            while (*p && *p != '"') p++;
            value = strndup_safe(val_start, (size_t)(p - val_start));
            if (*p == '"') p++;
        } else {
            const char *val_start = p;
            while (*p && *p != ';' && *p != ' ' && *p != '\t') p++;
            value = strndup_safe(val_start, (size_t)(p - val_start));
        }
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

        const char *tok_end = p;
        while (*tok_end && *tok_end != ';') tok_end++;

        const char *eq = p;
        while (eq < tok_end && *eq != '=') eq++;

        char *key, *value;
        if (eq == tok_end) {
            key   = strndup_safe(p, (size_t)(tok_end - p));
            value = strndup_safe("", 0);
        } else {
            key   = strndup_safe(p,      (size_t)(eq - p));
            value = strndup_safe(eq + 1, (size_t)(tok_end - eq - 1));
        }
        if (!key) { gtf_attrs_free(a); return NULL; }
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

        p = tok_end;
        if (*p == ';') p++;
    }

    return a;
}

/* ------------------------------------------------------------------ */
/* Printer                                                              */
/* ------------------------------------------------------------------ */

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
