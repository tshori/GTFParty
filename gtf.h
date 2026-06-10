#ifndef GTF_H
#define GTF_H

#include <stddef.h>

typedef enum { GTF_FMT_GTF = 0, GTF_FMT_GFF3 = 1 } GtfFormat;

#define GTF_FIELD_MAX  256
#define GTF_ATTR_MAX  8192

/*
 * Raw parse buffer filled by gtf_parse_line().
 * Fixed-size arrays — lives on the stack during parsing only; never stored
 * long-term.  stats.c and compare.c use this type directly since they parse
 * line-by-line without building a GtfTable.
 */
typedef struct {
    char  seqname[GTF_FIELD_MAX];
    char  source[GTF_FIELD_MAX];
    char  feature[GTF_FIELD_MAX];
    long  start;
    long  end;
    float score;
    char  strand;
    int   frame;
    char  attributes[GTF_ATTR_MAX];
} GtfRawRecord;

int  gtf_parse_line(const char *line, GtfRawRecord *rec);

/*
 * Compact record stored in GtfTable.
 *   seqname / source / feature  — interned in the table's StrPool
 *                                  (pointer-owned by the pool, not this struct)
 *   attributes                  — heap-allocated string owned by this record
 *
 * Interning means equal strings share one allocation, so pointer equality
 * implies string equality: rec_a->seqname == rec_b->seqname is O(1).
 *
 * GTF coordinates: 1-based, fully closed [start, end].
 * score == NAN   means the field was '.' in the file.
 * frame == -1    means the field was '.' in the file.
 */
typedef struct {
    const char *seqname;
    const char *source;
    const char *feature;
    long  start;
    long  end;
    float score;
    char  strand;
    int   frame;
    char *attributes;
} GtfRecord;

void gtf_print(const GtfRecord *rec);

/* ------------------------------------------------------------------ */
/* String intern pool                                                    */
/* ------------------------------------------------------------------ */

/*
 * Dynamic array of heap-allocated strings.  Stores one copy of each unique
 * string and returns the same pointer for every duplicate.
 * Owned and freed by the containing GtfTable.
 */
typedef struct {
    char  **data;
    size_t  n;
    size_t  cap;
} StrPool;

/* ------------------------------------------------------------------ */
/* Open-addressing hash table  (const char * key → size_t value)       */
/* ------------------------------------------------------------------ */

/*
 * Each slot holds one key/value pair.  key == NULL means the slot is empty.
 * cap is always a power of 2 so (hash & (cap-1)) replaces modulo.
 *
 * Because keys are interned strings, equality in the inner scan loop can
 * use pointer comparison instead of strcmp — the pointer IS the identity.
 */
typedef struct { const char *key; size_t val; } HEntry;

typedef struct {
    HEntry *slots;
    size_t  cap;
    size_t  n;
} HTable;

/* ------------------------------------------------------------------ */
/* GtfTable: dynamic array of GtfRecord                                */
/* ------------------------------------------------------------------ */

typedef struct {
    GtfRecord *recs;
    size_t     n;
    size_t     cap;
    size_t     n_skipped;
    size_t     n_errors;
    GtfFormat  format;
    StrPool    pool;    /* owns all interned seqname/source/feature strings */
    HTable    *index;   /* seqname → first sorted index; NULL until sorted */
} GtfTable;

GtfTable   *gtf_table_load(FILE *fp);
void        gtf_table_free(GtfTable *t);

/* Sort records in place by (seqname, start, end) and build the seqname index. */
void        gtf_table_sort(GtfTable *t);

/*
 * Find all records overlapping [qstart, qend] on seqname.
 * Table must be sorted first.  Overlap: rec.start <= qend && rec.end >= qstart.
 * Returns heap array of pointers into t->recs (*count entries), or NULL when
 * count == 0.  Caller frees the array; the pointed-to records belong to the table.
 */
GtfRecord **gtf_table_query(const GtfTable *t, const char *seqname,
                             long qstart, long qend, size_t *count);

/* ------------------------------------------------------------------ */
/* Attribute parsing                                                    */
/* ------------------------------------------------------------------ */

typedef struct { char *key; char *value; } GtfAttr;
typedef struct { GtfAttr *pairs; size_t n; } GtfAttrs;

GtfAttrs   *gtf_attrs_parse(const char *raw);
GtfAttrs   *gff3_attrs_parse(const char *raw);
const char *gtf_attrs_get(const GtfAttrs *a, const char *key);
void        gtf_attrs_free(GtfAttrs *a);

#endif /* GTF_H */
