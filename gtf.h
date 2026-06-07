#ifndef GTF_H
#define GTF_H

/* Fixed-size buffers for now — Step 2 will introduce heap allocation */
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

#endif /* GTF_H */
