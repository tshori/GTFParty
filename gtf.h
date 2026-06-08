#ifndef GTF_H
#define GTF_H

#include <stddef.h>

#define GTF_FIELD_MAX  256
#define GTF_ATTR_MAX  8192

/*
 * One record from a GTF file.
 * GTF coordinates are 1-based, fully closed [start, end].
 * score == NAN means the field was '.' in the file.
 * frame == -1  means the field was '.' in the file.
 */
typedef struct {
    char  seqname[GTF_FIELD_MAX];
    char  source[GTF_FIELD_MAX];
    char  feature[GTF_FIELD_MAX];
    long  start;
    long  end;
    float score;
    char  strand;    /* '+', '-', or '.' */
    int   frame;     /* 0, 1, 2, or -1  */
    char  attributes[GTF_ATTR_MAX];
} GtfRecord;

/*
 * Parse one line of a GTF file into *rec.
 *
 * Returns:
 *   0   success
 *   1   line is a comment or blank — skip it, *rec is untouched
 *  -1   parse error
 */
int  gtf_parse_line(const char *line, GtfRecord *rec);

/* Print a record back as a GTF line (tab-separated, newline at end) */
void gtf_print(const GtfRecord *rec);

/*
 * Heap-allocated dynamic array of GtfRecord values.
 * n_skipped counts comment/blank lines; n_errors counts malformed lines.
 */
typedef struct {
    GtfRecord *recs;
    size_t     n;
    size_t     cap;
    size_t     n_skipped;
    size_t     n_errors;
} GtfTable;

/*
 * Read all records from fp into a newly allocated GtfTable.
 * Returns NULL on allocation failure.
 * Caller must free with gtf_table_free().
 */
GtfTable *gtf_table_load(FILE *fp);

void      gtf_table_free(GtfTable *t);

/* ------------------------------------------------------------------ */
/* Attribute parsing                                                    */
/* ------------------------------------------------------------------ */

/* One key/value pair from the attributes column.  Both strings are    */
/* heap-allocated and owned by the containing GtfAttrs.                */
typedef struct {
    char *key;
    char *value;
} GtfAttr;

/* Parsed attributes column: a dynamic array of GtfAttr pairs.        */
typedef struct {
    GtfAttr *pairs;
    size_t   n;
} GtfAttrs;

/*
 * Parse the raw GTF attributes string into a heap-allocated GtfAttrs.
 * Handles both quoted values (key "value";) and unquoted (key value;).
 * Returns NULL on allocation failure.
 * Caller must free with gtf_attrs_free().
 */
GtfAttrs   *gtf_attrs_parse(const char *raw);

/*
 * Look up a key in a parsed GtfAttrs.
 * Returns a pointer to the value string, or NULL if the key is absent.
 * The pointer is valid for the lifetime of the GtfAttrs.
 */
const char *gtf_attrs_get(const GtfAttrs *a, const char *key);

void        gtf_attrs_free(GtfAttrs *a);

#endif /* GTF_H */
