#ifndef GTF_H
#define GTF_H

#include <stddef.h>

/* File format — detected automatically from the ##gff-version pragma. */
typedef enum { GTF_FMT_GTF = 0, GTF_FMT_GFF3 = 1 } GtfFormat;

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
    GtfFormat  format;    /* GTF_FMT_GTF or GTF_FMT_GFF3 */
} GtfTable;

/*
 * Read all records from fp into a newly allocated GtfTable.
 * Returns NULL on allocation failure.
 * Caller must free with gtf_table_free().
 */
GtfTable *gtf_table_load(FILE *fp);

void      gtf_table_free(GtfTable *t);

/* Sort all records in place by (seqname, start, end). */
void gtf_table_sort(GtfTable *t);

/*
 * Find all records that overlap [qstart, qend] on seqname.
 * The table must be sorted with gtf_table_sort() before calling this.
 * Overlap condition (1-based closed): rec.start <= qend && rec.end >= qstart.
 *
 * Returns a heap-allocated array of pointers into t->recs, with *count
 * set to the number of hits.  Returns NULL when *count == 0.
 * The caller must free the array itself (not the records it points to).
 * The pointers are only valid while the GtfTable lives.
 */
GtfRecord **gtf_table_query(const GtfTable *t, const char *seqname,
                             long qstart, long qend, size_t *count);

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
 * Parse GTF attributes  (key "value"; key "value"; ...)  into a GtfAttrs.
 * Parse GFF3 attributes (key=value;key=value;...)         into a GtfAttrs.
 * Both return the same GtfAttrs type; free with gtf_attrs_free().
 * Return NULL on allocation failure.
 */
GtfAttrs   *gtf_attrs_parse(const char *raw);
GtfAttrs   *gff3_attrs_parse(const char *raw);

/*
 * Look up a key in a parsed GtfAttrs.
 * Returns a pointer to the value string, or NULL if the key is absent.
 * The pointer is valid for the lifetime of the GtfAttrs.
 */
const char *gtf_attrs_get(const GtfAttrs *a, const char *key);

void        gtf_attrs_free(GtfAttrs *a);

#endif /* GTF_H */
