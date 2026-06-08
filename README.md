# GTFing — GTF/GFF Utilities in C

A set of command-line utilities for processing GTF and GFF genome annotation
files, built step by step as a C learning project.

## Current tools

| Tool                  | Language | Description |
|-----------------------|----------|-------------|
| `gtfparse print`      | C        | Parse a GTF file and print it back out (round-trip validator) |
| `gtfparse filter`     | C        | Print only records matching a given feature type |
| `gtfparse attrs`      | C        | Extract one attribute value per record |
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
./gtfparse <command> [file.gtf ...]
```

### print

Parses every record and writes it back to stdout.
Diagnostic counts (parsed / skipped / errors) go to stderr.

```sh
./gtfparse print file.gtf

# Round-trip test — diff should be empty
./gtfparse print file.gtf > out.gtf
diff <(grep -v '^#' file.gtf) out.gtf

# Count records by feature type
./gtfparse print file.gtf | awk '{print $3}' | sort | uniq -c | sort -rn
```

### filter

Prints all records whose feature column matches the given type.
Diagnostic counts (matched / total) go to stderr.

```sh
./gtfparse filter exon file.gtf
./gtfparse filter CDS  file.gtf
./gtfparse filter gene file.gtf > genes.gtf
```

### attrs

Parses the attributes column of every record and prints the value of the
requested key, one line per record.  Prints `.` when the key is absent.

```sh
# List all unique gene IDs
./gtfparse attrs gene_id file.gtf | sort -u

# Extract transcript IDs from exon records only
./gtfparse filter exon file.gtf | ./gtfparse attrs transcript_id /dev/stdin

# Count records per gene
./gtfparse attrs gene_id file.gtf | sort | uniq -c | sort -rn | head
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
5. **Step 5** — Sort records by position
6. **Step 6** — GFF3 support
7. **Step 7** — Interval overlap query
