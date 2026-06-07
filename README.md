# GTFing — GTF/GFF Utilities in C

A set of command-line utilities for processing GTF and GFF genome annotation
files, built step by step as a C learning project.

## Current tools

| Binary     | Description                          |
|------------|--------------------------------------|
| `gtfparse` | Parse a GTF file and print it back out (round-trip validator) |

## Build

Requires GCC and GNU Make.

```sh
make          # builds gtfparse
make clean    # removes object files and binaries
```

## Usage

```sh
./gtfparse <file.gtf>
```

Reads the GTF file, parses every record, and writes it back to stdout.
Diagnostic counts (parsed / skipped / errors) are written to stderr so they
do not pollute redirected output.

```sh
# Round-trip test — diff should be empty
./gtfparse file.gtf > out.gtf
diff <(grep -v '^#' file.gtf) out.gtf

# Count records by feature type
./gtfparse file.gtf | awk '{print $3}' | sort | uniq -c | sort -rn
```

## Project structure

```
.
├── gtf.h       # GtfRecord struct and function declarations
├── gtf.c       # Parser and printer implementation
├── main.c      # Entry point — file I/O loop
├── Makefile    # Build rules
├── README.md   # This file
└── MANUAL.md   # Detailed code walkthrough
```

## Tested on

- `GCF_963676685.1_xbMytEdul2.2_genomic.gtf` — NCBI RefSeq GTF for
  *Mytilus edulis* (blue mussel), 1,725,022 records, 0 parse errors.

## Roadmap

1. **Step 1** ✓ — Parse a GTF record into a struct, round-trip printer
2. **Step 2** — Dynamic array of records (heap allocation)
3. **Step 3** — Filter by feature type (CLI flags)
4. **Step 4** — Parse the attributes column into key/value pairs
5. **Step 5** — Sort records by position
6. **Step 6** — GFF3 support
7. **Step 7** — Interval overlap query
