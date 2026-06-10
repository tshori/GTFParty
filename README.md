# GTFing — GTF/GFF Utilities in C

A set of command-line utilities for processing GTF and GFF genome annotation
files, built step by step as a C learning project.

## Current tools

| Tool                  | Language | Description |
|-----------------------|----------|-------------|
| `gtfparse print`      | C        | Parse a GTF/GFF3 file and print it back out (round-trip validator) |
| `gtfparse filter`     | C        | Print only records matching a given feature type |
| `gtfparse attrs`      | C        | Extract one attribute value per record (GTF and GFF3) |
| `gtfparse overlap`    | C        | Find all records overlapping a genomic region |
| `gtfparse bed`        | C        | Convert GTF/GFF3 to BED6 (0-based half-open coordinates) |
| `gtfparse sort`       | C        | Sort records by (seqname, start, end) |
| `gtfparse keys`       | C        | List all unique attribute keys in a file |
| `gtfparse validate`   | C        | Check coordinate sanity, strand, frame, and required attributes |
| `gtfparse stats`      | C        | Exon and intron statistics: counts, per-gene averages, single-exon rates |
| `gtfparse compare`    | C        | Structural comparison of two GTFs by genomic coordinate and intron chain |
| `compare_gtf.py`      | Python   | Compare two GTFs by gene/transcript ID; write filtered or merged output; Venn diagram |

## Build

**C tools** — requires GCC and GNU Make:

```sh
make          # builds gtfparse
make clean    # removes object files and binaries
```

**Python script** — no build step; optional dependency for `--venn`:

```sh
pip install matplotlib matplotlib-venn
```

## Usage

```sh
./gtfparse <command> [args] [-o <file>]
```

### Global flags

| Flag | Short | Effect |
|------|-------|--------|
| `--output <file>` | `-o <file>` | Write output to `<file>` instead of stdout |

The `-o` flag works with every command and can appear anywhere in the argument list.
Diagnostic messages always go to stderr regardless of `-o`, so progress and
summary lines remain visible on the terminal while output goes to the file.

```sh
# Equivalent to shell redirection, but built in
./gtfparse filter gene file.gtf -o genes.gtf
./gtfparse bed --name gene_id file.gtf -o genes.bed
./gtfparse print --chrom chr1 -o chr1.gtf file.gtf

# --output=<file> form also accepted
./gtfparse filter exon file.gtf --output=exons.gtf
```

### print

Parses every record and writes it back to stdout.
Accepts both GTF and GFF3 — format is detected automatically from the
`##gff-version 3` pragma.  Diagnostic counts (including a `printed=` tally)
go to stderr.  Loading progress is printed to stderr every 100,000 records.

Optional flags narrow output to one or more chromosomes or sources:

| Flag | Value | Effect |
|------|-------|--------|
| `--chrom` | comma-separated seqnames | keep only records on these chromosomes |
| `--source` | comma-separated source names | keep only records from these sources |

Both flags can be combined (AND logic: record must pass both).

```sh
./gtfparse print file.gtf
./gtfparse print file.gff3

# Single chromosome
./gtfparse print --chrom chr1 file.gtf

# Multiple chromosomes
./gtfparse print --chrom chr1,chr2,chrX file.gtf > subset.gtf

# Filter by source/origin
./gtfparse print --source RefSeq file.gtf

# Combine: chrX records from Gnomon only
./gtfparse print --chrom chrX --source Gnomon file.gtf

# Round-trip test — diff should be empty
./gtfparse print file.gtf > out.gtf
diff <(grep -v '^#' file.gtf) out.gtf

# Count records by feature type
./gtfparse print file.gff3 | awk '{print $3}' | sort | uniq -c | sort -rn
```

### filter

Prints all records whose feature column matches the given type.
Works with both GTF and GFF3.  Diagnostic counts go to stderr.

```sh
./gtfparse filter exon file.gtf
./gtfparse filter CDS  file.gff3
./gtfparse filter gene file.gtf > genes.gtf
```

### overlap

Finds all records whose `[start, end]` interval overlaps a query region.
Accepts `seqname:start-end` in 1-based fully-closed coordinates (the GTF/GFF3
native coordinate system).  Loads the file into memory, sorts by position,
then binary-searches to the target chromosome before scanning.

```sh
./gtfparse overlap chr1:1000000-2000000 file.gtf
./gtfparse overlap NC_092344.1:31000-32000 file.gff3

# How many features overlap a region?
./gtfparse overlap chr1:1-5000000 file.gtf 2>&1 | grep ^hits
```

### attrs

Parses the attributes column of every record and prints the value of the
requested key, one line per record.  Prints `.` when the key is absent.
Automatically uses GTF (`key "value";`) or GFF3 (`key=value;`) parsing
based on the detected format.

```sh
# GTF: list all unique gene IDs
./gtfparse attrs gene_id file.gtf | sort -u

# GFF3: list all unique feature IDs
./gtfparse attrs ID file.gff3 | sort -u

# Extract Parent IDs from exon records
./gtfparse filter exon file.gff3 | ./gtfparse attrs Parent /dev/stdin

# Count records per gene
./gtfparse attrs gene_id file.gtf | sort | uniq -c | sort -rn | head
```

**Note on GFF3 multi-value attributes:** GFF3 allows comma-separated values
(e.g. `Parent=mRNA:t1,mRNA:t2`).  The entire comma-separated string is
returned as a single value; splitting is left to the caller.

**Note on `stats` and `compare`:** these commands look for `gene_id` and
`transcript_id` attributes in GTF format and will not produce meaningful
results on GFF3 files.

### bed

Converts every record to BED6 format using 0-based half-open coordinates
(the standard expected by bedtools, IGV, and other genome tools).

| Column | Value |
|--------|-------|
| chrom  | seqname as-is |
| start  | GTF `start − 1` (0-based) |
| end    | GTF `end` (unchanged; now exclusive) |
| name   | value of `--name` attribute, or `.` |
| score  | GTF score (0–1000, integer), or `0` when `.` in file |
| strand | `+`, `-`, or `.` |

```sh
# BED6 with no name column
./gtfparse bed file.gtf > file.bed

# Use gene_id as the name (GTF)
./gtfparse bed --name gene_id file.gtf > genes.bed

# Use ID as the name (GFF3)
./gtfparse bed --name ID file.gff3 > features.bed

# Gene records only, named by gene_id — common pipeline step
./gtfparse filter gene file.gtf | ./gtfparse bed --name gene_id /dev/stdin > genes.bed

# Check with bedtools
bedtools sort -i genes.bed | head
```

Coordinate conversion example: a GTF record `chr1 . gene 1001 2000` becomes
BED `chr1 1000 2000 . 0 .`.

### sort

Sorts all records by `(seqname, start, end)` and writes them back out in the
original format.  The sort order is lexicographic on seqname, then numeric on
start and end.

Sorting is required before `overlap` queries, and is useful any time
downstream tools expect a canonical record order.  Diagnostic counts go to
stderr; sorted records go to stdout.

```sh
./gtfparse sort file.gtf > sorted.gtf
./gtfparse sort file.gff3 -o sorted.gff3

# Sort then query — overlap sorts internally, but explicit sort lets you
# inspect the sorted file or feed it to other tools
./gtfparse sort file.gtf | ./gtfparse overlap chr1:1000-5000 /dev/stdin
```

### keys

Scans every record in a file, parses the attributes column, and prints all
unique attribute key names — one per line, sorted alphabetically.  Works with
both GTF (`key "value";`) and GFF3 (`key=value;`).

Useful for discovering what metadata is available before querying with `attrs`,
or for auditing two annotation files to see if they share the same attribute
vocabulary.

```sh
./gtfparse keys file.gtf
./gtfparse keys file.gff3

# Compare attribute keys between two files
diff <(./gtfparse keys file1.gtf 2>/dev/null) <(./gtfparse keys file2.gtf 2>/dev/null)

# Save key list to a file
./gtfparse keys file.gtf -o keys.txt
```

### validate

Checks every record for common errors and prints each violation to stdout,
one per line.  Exits with code `0` if no errors are found, `1` otherwise —
making it usable directly in shell scripts.

**Checks performed:**

| Check | What is tested |
|-------|---------------|
| `coord` | `start >= 1`; `end >= start` |
| `strand` | value is `+`, `-`, or `.` |
| `frame` | CDS / start_codon / stop_codon must have a numeric frame (0, 1, 2) |
| `attr` (GTF) | every record has `gene_id`; exon/CDS/codon features also have `transcript_id` |
| `attr` (GFF3) | gene/mRNA/transcript features have `ID`; mRNA/exon/CDS features have `Parent` |

Format (GTF vs GFF3) is detected automatically.  Violations go to stdout so
they can be saved with `-o` or piped; the summary line goes to stderr.

```sh
./gtfparse validate file.gtf
./gtfparse validate file.gff3

# Save violations to a file, still see summary on terminal
./gtfparse validate file.gtf -o errors.txt

# Use exit code in a script
if ./gtfparse validate file.gtf > /dev/null 2>&1; then
    echo "file is valid"
else
    echo "validation failed"
fi

# Count violations by check type
./gtfparse validate file.gtf 2>/dev/null | awk -F': ' '{print $2}' | sort | uniq -c
```

Example output for a file with problems:
```
record 3: coord: start (0) must be >= 1
record 7: strand: invalid value 'q'
record 12: frame: feature 'CDS' requires a frame (0, 1, or 2)
record 41: attr: missing required attribute 'gene_id'
```

### stats

Reads all `exon` features and computes:

- Total genes, transcripts, exons, and introns
- Mean exons and introns per gene
- Number and percentage of single-exon genes
- Number and percentage of single-exon transcripts

Introns are inferred as the gaps between exons within each transcript
(`exons − 1` per transcript); `intron` feature lines are not used.
A single-exon gene is one where **every** transcript has exactly 1 exon.

```sh
./gtfparse stats file.gtf
```

Example output:
```
Genes                   : 14234
Transcripts             : 32891
Exons                   : 187442
Introns                 : 154551
Mean exons per gene     : 13.17
Mean introns per gene   : 10.86
Single-exon genes       : 1823 (12.8%)
Single-exon transcripts : 4201 (12.8%)
```

### compare

Compares two GTF files by **genomic coordinate and intron chain structure**,
not by ID strings. This is the right approach when comparing two distinct
gene predictors (e.g. AUGUSTUS vs StringTie) where IDs will differ.

Each transcript in GTF1 is classified against every overlapping transcript
in GTF2:

| Code | Meaning |
|------|---------|
| `=`  | Exact intron chain match — all splice sites identical |
| `j`  | Partial match — at least one shared intron, chains differ |
| `o`  | Positional overlap, no shared introns (e.g. two single-exon genes) |
| `u`  | No overlap with any transcript in the other file |

```sh
./gtfparse compare pred1.gtf pred2.gtf
```

Example output:
```
=== Transcript classification (GTF1) ===
  GTF1                      :    32891
  Exact match      (=)      :    18432  ( 56.0%)
  Partial match    (j)      :     6201  ( 18.9%)
  Overlap only     (o)      :     1834  (  5.6%)
  Unique           (u)      :     6424  ( 19.5%)

=== Transcript classification (GTF2) ===
  ...

=== Intron-level ===
  GTF1 unique introns       :   154551
  GTF2 unique introns       :   141203
  Shared                    :   112834  (73.0% of GTF1, 79.9% of GTF2)

  Sensitivity (GTF1 exact / GTF1 total) : 56.0%
  Precision   (GTF2 exact / GTF2 total) : 62.3%
```

## `compare_gtf.py` usage

```sh
python3 compare_gtf.py <gtf1> <gtf2> [options]
```

Comparison is by ID string (`gene_id` / `transcript_id`), not coordinates.
Summary is always written to stderr; stdout is reserved for GTF output.

```sh
# Summary only
python3 compare_gtf.py v110.gtf v111.gtf

# Summary + write genes present in v110 but absent from v111
python3 compare_gtf.py v110.gtf v111.gtf --output only-in-1 --out-file lost_genes.gtf

# Summary + write merged annotation (v110 + unique v111 genes)
python3 compare_gtf.py v110.gtf v111.gtf --output merge --out-file merged.gtf

# Two-panel Venn diagram (genes + transcripts)
python3 compare_gtf.py v110.gtf v111.gtf --venn venn.png

# Custom labels in the output
python3 compare_gtf.py v110.gtf v111.gtf --label1 "Ensembl v110" --label2 "Ensembl v111" --venn venn.png
```

Output modes:

| `--output`    | Records written |
|---------------|-----------------|
| `merge`       | All of GTF1, plus any gene from GTF2 absent in GTF1 |
| `only-in-1`   | GTF1 records whose `gene_id` is absent from GTF2 |
| `only-in-2`   | GTF2 records whose `gene_id` is absent from GTF1 |

## Project structure

```
.
├── gtf.h            # GtfRecord, GtfTable, GtfAttrs structs and declarations
├── gtf.c            # Parser, printer, dynamic record table, attribute parsing
├── stats.h          # GtfExonIntronStats struct and function declarations
├── stats.c          # Exon/intron statistics implementation
├── compare.h        # Transcript struct, CmpStats, function declarations
├── compare.c        # Structural comparison implementation
├── main.c           # Entry point — subcommand dispatch
├── Makefile         # Build rules
├── compare_gtf.py   # ID-based comparison script (Python)
├── requirements.txt # Python dependencies for --venn
├── README.md        # This file
└── MANUAL.md        # Detailed C code walkthrough
```

## Tested on

- `GCF_963676685.1_xbMytEdul2.2_genomic.gtf` — NCBI RefSeq GTF for
  *Mytilus edulis* (blue mussel), 1,725,022 records, 0 parse errors.

## Roadmap

1. **Step 1** ✓ — Parse a GTF record into a struct, round-trip printer;
   exon/intron stats (heap allocation, `qsort`, attribute parsing);
   ID-based GTF comparison with Venn diagram (`compare_gtf.py`);
   coordinate-based structural comparison (sweep line, intron chain matching)
2. **Step 2** ✓ — Dynamic array of all records (`GtfTable`: heap-allocated, doubling-growth array; `gtf_table_load` / `gtf_table_free`)
3. **Step 3** ✓ — Filter by feature type (`gtfparse filter <feature> <file.gtf>`)
4. **Step 4** ✓ — Parse the attributes column into key/value pairs (`GtfAttr`/`GtfAttrs`; `gtf_attrs_parse` / `gtf_attrs_get` / `gtf_attrs_free`)
5. **Step 5** ✓ — Sort records by position (`gtf_table_sort`: in-place `qsort` by seqname/start/end)
6. **Step 6** ✓ — GFF3 support (`gff3_attrs_parse`; auto-detect via `##gff-version`; stop at `##FASTA`)
7. **Step 7** ✓ — Interval overlap query (`gtf_table_query`: binary search to seqname block + forward scan; `gtfparse overlap seqname:start-end`)
8. **Step 8** ✓ — BED6 export (`gtfparse bed [--name attr]`: GTF→BED coordinate conversion, optional attribute-based name extraction)
9. **Step 9** ✓ — Sort command (`gtfparse sort`: expose `gtf_table_sort` as a user-facing command; library/command separation)
10. **Step 10** ✓ — Attribute key discovery (`gtfparse keys`: dynamic string array, deduplication by linear scan, `qsort` on `char **`, `goto` for OOM cleanup)
11. **Step 11** ✓ — Validation (`gtfparse validate`: per-record checks, NULL-terminated feature sets, meaningful exit codes)
