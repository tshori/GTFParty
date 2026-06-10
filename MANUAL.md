# GTFing — Code Manual

A detailed walkthrough of every design decision and C concept in the codebase.
Written alongside the code as a learning companion.

---

## Table of contents

1. [The GTF format](#1-the-gtf-format)
2. [Project layout and the build system](#2-project-layout-and-the-build-system)
3. [gtf.h — the header file](#3-gtfh--the-header-file)
4. [gtf.c — the implementation](#4-gtfc--the-implementation)
5. [main.c — the entry point](#5-mainc--the-entry-point)
6. [stats.h and stats.c — the analytics layer](#6-statsh-and-statsc--the-analytics-layer)
7. [compare.h and compare.c — the comparison layer](#7-compareh-and-comparec--the-comparison-layer)
8. [Pointers — a complete picture](#8-pointers--a-complete-picture)
9. [Memory lifetimes and dangling pointers](#9-memory-lifetimes-and-dangling-pointers)
10. [Current limitations](#10-current-limitations)
11. [GtfTable — a dynamic array of records](#11-gtftable--a-dynamic-array-of-records)
12. [filter — using GtfTable in a command](#12-filter--using-gtftable-in-a-command)
13. [GtfAttrs — parsing the attributes column](#13-gtfattrs--parsing-the-attributes-column)
14. [GFF3 support — a second file format](#14-gff3-support--a-second-file-format)
15. [gtf_table_sort — sorting records by position](#15-gtf_table_sort--sorting-records-by-position)
16. [gtf_table_query — interval overlap queries](#16-gtf_table_query--interval-overlap-queries)
17. [Interval merge — a scan on sorted data](#17-interval-merge--a-scan-on-sorted-data)
18. [String interning — one allocation per unique string](#18-string-interning--one-allocation-per-unique-string)
19. [Hash table — O(1) seqname lookup](#19-hash-table--o1-seqname-lookup)
20. [BED output — coordinate conversion and format interop](#20-bed-output--coordinate-conversion-and-format-interop)
21. [Output redirection — argv compaction and freopen](#21-output-redirection--argv-compaction-and-freopen)
22. [sort — exposing library code as a command](#22-sort--exposing-library-code-as-a-command)
23. [keys — dynamic string deduplication and goto](#23-keys--dynamic-string-deduplication-and-goto)
24. [validate — per-record checks, feature sets, and exit codes](#24-validate--per-record-checks-feature-sets-and-exit-codes)

---

## 1. The GTF format

GTF (Gene Transfer Format) is a tab-delimited text format for genome
annotations. Each data line has exactly 9 columns:

```
seqname  source  feature  start  end  score  strand  frame  attributes
```

Example:
```
NC_092344.1  Gnomon  exon  31352  31529  .  -  .  gene_id "LOC139519001"; transcript_id "XM_071310724.1";
```

Key rules:
- Lines beginning with `#` are comments and must be skipped.
- Coordinates are **1-based** and **fully closed** — both start and end are
  included in the feature. This differs from BED format, which is 0-based
  half-open.
- `score` and `frame` use `.` to mean "not applicable."
- The `attributes` column contains `key "value";` pairs separated by
  semicolons (GTF2 format). Values are always double-quoted.
- Fields are separated by tabs. The `attributes` field may contain spaces
  but never tabs, so splitting on `\t` cleanly isolates all 9 columns.

---

## 2. Project layout and the build system

### Four-layer architecture

The project is split into four layers, each in its own pair of files:

```
main.c                  ← CLI layer: reads arguments, opens files, dispatches
stats.c   / stats.h     ← analytics layer: per-gene exon/intron statistics
compare.c / compare.h   ← comparison layer: structural GTF comparison
gtf.c     / gtf.h       ← parser layer: parses and prints individual records
```

Each layer only depends on the layers below it. `main.c` knows about
`stats.h`, `compare.h`, and `gtf.h`. Both `stats.c` and `compare.c` know
about `gtf.h` (they use `GtfRecord` and `gtf_parse_line`). `gtf.c` knows
nothing about any of the others.

This separation has a practical consequence: you can change how parsing works
inside `gtf.c` without touching `main.c` or `stats.c`, as long as the
function signatures in `gtf.h` stay the same. The compiler enforces this
boundary — if `main.c` tries to call a function that isn't declared in a
header it includes, it will refuse to compile.

### The Makefile

```makefile
CC      = gcc
CFLAGS  = -Wall -Wextra -std=c11 -g
LDFLAGS = -lm
```

- `CC` — the compiler. Assigning it to a variable means you can override it
  from the command line: `make CC=clang`.
- `-Wall -Wextra` — enable almost all compiler warnings. Treat these as
  errors in practice; they catch real bugs.
- `-std=c11` — compile to the C11 standard. This ensures consistent behaviour
  across compilers and machines.
- `-g` — include debug symbols in the output. This lets debuggers (gdb,
  lldb) show you source line numbers and variable names instead of raw
  addresses.
- `-lm` — link the math library. Required for `NAN` and `isnan()`. The `-l`
  flag must come **after** the object files that use it because the linker
  resolves symbols left to right.

### Two-stage compilation

```makefile
gtfparse: main.o gtf.o stats.o compare.o
    $(CC) $(CFLAGS) -o gtfparse main.o gtf.o stats.o compare.o $(LDFLAGS)

main.o: main.c gtf.h stats.h compare.h
    $(CC) $(CFLAGS) -c main.c

gtf.o: gtf.c gtf.h
    $(CC) $(CFLAGS) -c gtf.c

stats.o: stats.c gtf.h stats.h
    $(CC) $(CFLAGS) -c stats.c

compare.o: compare.c gtf.h compare.h
    $(CC) $(CFLAGS) -c compare.c
```

**Stage 1 — compile** (`-c` flag): each `.c` file is translated into machine
code and written to a `.o` (object) file. At this stage, calls to functions
in other files are left as unresolved placeholders.

**Stage 2 — link**: the linker combines the `.o` files, resolves every
placeholder, and writes the final executable.

The benefit of this split: `make` tracks which files have changed. If you
edit only `gtf.c`, only `gtf.o` is recompiled. `main.o` and `stats.o` are
reused as-is. On large projects with hundreds of files this matters
enormously.

Each rule also lists its header dependencies. `main.o: main.c gtf.h stats.h`
tells `make` that if either header changes, `main.c` must be recompiled.
Without this, you could edit a struct definition in a header and end up with
object files built from different versions of that struct — a very subtle bug.

You can inspect what symbols a `.o` file exports and imports:

```sh
nm stats.o
```

Output:
```
T gtf_exon_intron_stats       # T = defined (exported) by this file
T gtf_print_exon_intron_stats
U gtf_parse_line              # U = undefined (imported from elsewhere)
U malloc
U realloc
U free
U qsort
U strcmp
```

The `static` helper functions (`attr_get`, `cmp_by_tx`, `cmp_by_gene`) do
not appear in this list because they are invisible outside the file.
See section 6 for a detailed discussion.

---

## 3. `gtf.h` — the header file

### Include guards

```c
#ifndef GTF_H
#define GTF_H
// ...
#endif
```

If two `.c` files both `#include "gtf.h"`, the preprocessor would process the
header twice. The second pass would try to re-declare the struct and functions,
causing a compilation error. The include guard prevents this:

- First encounter: `GTF_H` is not defined → enter the block, define it,
  process the contents.
- Second encounter: `GTF_H` is already defined → skip everything to `#endif`.

### Preprocessor constants

```c
#define GTF_FIELD_MAX  256
#define GTF_ATTR_MAX  8192
```

`#define` is a preprocessor directive, not a C statement. Before the compiler
sees your code, the preprocessor performs a global text substitution:
every occurrence of `GTF_FIELD_MAX` is replaced with `256`. These constants
have no type, no address, and consume no memory. The advantage over raw
numbers ("magic numbers") is that you can change the buffer size in one
place and all uses update automatically.

### Two record types

The header defines two related structs — one for parsing, one for storage.

**`GtfRawRecord`** is the scratch buffer filled by `gtf_parse_line`. It uses
fixed-size character arrays and lives on the stack for one loop iteration:

```c
typedef struct {
    char  seqname[GTF_FIELD_MAX];   /* 256 bytes */
    char  source[GTF_FIELD_MAX];    /* 256 bytes */
    char  feature[GTF_FIELD_MAX];   /* 256 bytes */
    long  start;
    long  end;
    float score;
    char  strand;
    int   frame;
    char  attributes[GTF_ATTR_MAX]; /* 8192 bytes */
} GtfRawRecord;
```

Any string longer than the array minus one will be silently truncated.
In practice, NCBI RefSeq files fit within these limits.

**`GtfRecord`** is stored inside `GtfTable`. Its string fields are pointers
rather than embedded arrays:

```c
typedef struct {
    const char *seqname;   /* interned — owned by the table's pool */
    const char *source;    /* interned — owned by the table's pool */
    const char *feature;   /* interned — owned by the table's pool */
    long  start;
    long  end;
    float score;
    char  strand;
    int   frame;
    char *attributes;      /* heap-allocated — owned by this record */
} GtfRecord;
```

`seqname`, `source`, and `feature` point into a shared *string intern pool*
(§18): every record on chromosome 1 shares the same single allocation of
`"chr1"`, so pointer equality is string equality. `attributes` is unique per
record and allocated with `strndup_safe` during loading.

The size of a `GtfRecord` on disk is ~40 bytes of fixed fields plus a pointer.
The `attributes` string averages ~100 bytes. This is a **60× reduction** from
the ~9 KB that the original fixed-array layout required per record.

A `struct` groups related values into a single named unit. Without `typedef`
you would write `struct GtfRecord rec` everywhere. The `typedef` lets you
write just `GtfRecord rec`.

**Sentinel values** encode "not applicable" without a separate boolean flag:

| Field    | Type    | Normal values | Sentinel meaning "not applicable" |
|----------|---------|---------------|-----------------------------------|
| `score`  | `float` | any number    | `NAN` (IEEE 754 not-a-number)     |
| `frame`  | `int`   | 0, 1, 2       | `-1`                              |
| `strand` | `char`  | `'+'`, `'-'`  | `'.'`                             |

`NAN` is the standard C choice for a missing float because you cannot
accidentally compare it equal to anything — `NAN == NAN` is always false
(IEEE 754 rule). Callers must use `isnan()` to detect it.

### Function declarations

```c
int  gtf_parse_line(const char *line, GtfRawRecord *rec);
void gtf_print(const GtfRecord *rec);
```

These are **declarations** (also called prototypes), not definitions. They
tell the compiler what types the functions accept and return, so it can
type-check call sites in `main.c` before it has seen `gtf.c`. The actual
code lives in `gtf.c`.

`const char *line` — a pointer to char, with `const` on the data. The
function promises not to modify the string it receives.

`GtfRawRecord *rec` — a pointer to the caller's scratch buffer. The
function writes tab-parsed fields into the fixed arrays. Using
`GtfRawRecord` (not `GtfRecord`) here is intentional: the parser just
fills arrays; the conversion to interned pointer fields happens in
`gtf_table_load`. `stats.c` and `compare.c` use `GtfRawRecord` directly
and never build a `GtfTable`.

### GtfTable and related types

`gtf.h` declares several higher-level types built on top of `GtfRecord`:

**`StrPool`** is a dynamic array of `char *` — the intern pool that owns
all `seqname`/`source`/`feature` strings stored in a `GtfTable`. See §18.

**`HEntry` / `HTable`** form an open-addressing hash table mapping
interned `const char *` keys to `size_t` values. The table's `index` field
uses this to map each seqname to the first record index in the sorted
array, giving O(1) chromosome lookup. See §19.

**`GtfTable`** holds a heap-allocated, dynamically-grown array of every
record in a file, plus the intern pool and the hash index. See §11.

**`GtfAttr` / `GtfAttrs`** represent the parsed attributes column: a
dynamic array of `{char *key; char *value;}` pairs allocated on the heap
on demand. Rather than scanning the raw attribute string each time you need
a value, you parse once and look up by key. See §13.

---

## 4. `gtf.c` — the implementation

### `gtf_parse_line`

#### Skipping non-data lines

```c
if (line[0] == '#' || line[0] == '\n' || line[0] == '\0')
    return 1;
```

`line[0]` reads the first byte of the string. In C, array indexing and
pointer arithmetic are the same thing: `line[0]` is identical to `*(line + 0)`.
Character literals like `'#'` are just integers (their ASCII value).
Returning `1` is the "skip this line" signal to the caller, distinct from
`0` (success) and `-1` (error).

#### Copying into a local buffer

```c
char buf[65536];
strncpy(buf, line, sizeof(buf) - 1);
buf[sizeof(buf) - 1] = '\0';
```

`strtok` (used below) modifies the string it operates on — it replaces
delimiter characters with null bytes in place. If we called it directly on
`line`, we would modify the caller's memory, which is both unexpected and
dangerous (e.g. if `line` pointed to read-only memory, it would crash).
Copying into our own `buf` means we own the memory and can do whatever we
like with it.

`strncpy(dst, src, n)` copies at most `n` bytes. Critical gotcha: if `src`
is longer than `n` bytes, `strncpy` fills exactly `n` bytes and does **not**
write a null terminator. The manual null assignment on the next line is the
standard defensive pattern to guarantee the buffer is always null-terminated.

`sizeof(buf)` gives the total byte size of the array — `65536`. For a
pointer (`char *p`), `sizeof(p)` would give the size of the pointer itself
(8 bytes on a 64-bit system), not the size of whatever it points to. This
distinction is important.

#### Stripping the trailing newline

```c
size_t len = strlen(buf);
if (len > 0 && buf[len - 1] == '\n')
    buf[len - 1] = '\0';
```

`fgets` always includes the newline character in the buffer if one was
present. Overwriting it with `'\0'` shortens the string by one character.
`strlen` returns the number of characters before the null terminator, so
`buf[len - 1]` is the last real character.

#### Splitting on tabs with `strtok`

```c
char *fields[9];
int n = 0;
char *tok = strtok(buf, "\t");
while (tok && n < 9) {
    fields[n++] = tok;
    tok = strtok(NULL, "\t");
}
```

`strtok(str, delimiters)` works in two modes:

- **First call** — pass the string to search. `strtok` finds the first
  delimiter, replaces it with `'\0'`, and returns a pointer to the start of
  the token.
- **Subsequent calls** — pass `NULL`. `strtok` resumes from where it left off
  using an internal static pointer.

After the loop, `buf` looks like this (tabs replaced by null bytes):

```
Before:  [NC_092344.1\tGnomon\texon\t31352\t...]
After:   [NC_092344.1\0Gnomon\0exon\031352\0...]
                      ↑       ↑
                 fields[1]  fields[2]
```

`fields[i]` points into `buf` at the start of each field. No new memory is
allocated — this is efficient but means `fields` is only valid while `buf`
exists. See section 9 for a detailed discussion of this.

`strtok` returns `NULL` when there are no more tokens, which is why
`while (tok && ...)` terminates correctly. `NULL` is address zero, and zero
is falsy in C.

`n++` is post-increment: it evaluates to the current value of `n` (used as
the index), then increments `n`. So `fields[n++] = tok` stores to
`fields[0]`, `fields[1]`, etc.

One important limitation: `strtok` is **not thread-safe** because it stores
state in a static (global) variable. If two threads called it simultaneously
they would corrupt each other's progress. The thread-safe version is
`strtok_r`, which takes an explicit `char **saveptr` argument instead of
using hidden state.

#### Copying fields into the scratch buffer

```c
strncpy(rec->seqname, fields[0], GTF_FIELD_MAX - 1);
rec->seqname[GTF_FIELD_MAX - 1] = '\0';
```

`rec` is a pointer to the caller's `GtfRawRecord`. The `->` operator is
shorthand: `rec->seqname` means `(*rec).seqname` — dereference the pointer
to get the struct, then access the field. Writing to `rec->seqname` writes
directly into the caller's memory; no copy of the struct is made.

We copy the bytes that `fields[0]` points to into the struct's fixed array.
This is the crucial step that makes our data safe — once copied, the struct
fields no longer depend on `buf`, which disappears when the function returns.

`gtf_parse_line` fills only a `GtfRawRecord` (fixed arrays). The conversion
to a `GtfRecord` with interned pointer fields happens in `gtf_table_load`:
it calls `str_pool_intern` for the three short fields and `strndup_safe` for
`attributes`. See §18.

#### Parsing numeric and special fields

```c
rec->start = atol(fields[3]);
rec->end   = atol(fields[4]);
```

`atol` (ASCII to long) converts a string like `"31352"` to the integer
`31352`. It stops at the first non-numeric character and returns 0 on
failure (with no error reporting — `strtol` with error checking is the
robust alternative, but `atol` is fine here since we trust the GTF format).

```c
rec->score = (fields[5][0] == '.') ? NAN : (float)atof(fields[5]);
```

`fields[5][0]` — `fields[5]` is a `char *`, so `[0]` gives the first
character. If it's `'.'`, we use `NAN`. Otherwise `atof` parses the float,
and `(float)` casts the `double` result down to match the struct field type.

```c
rec->frame = (fields[7][0] == '.') ? -1 : atoi(fields[7]);
```

Same pattern. `-1` is the sentinel for "frame not applicable."

### `gtf_print`

```c
void gtf_print(const GtfRecord *rec)
{
    printf("%s\t%s\t%s\t%ld\t%ld\t", ...);
    if (isnan(rec->score)) printf("."); else printf("%g", rec->score);
    printf("\t%c\t", rec->strand);
    if (rec->frame == -1) printf("."); else printf("%d", rec->frame);
    printf("\t%s\n", rec->attributes);
}
```

`const GtfRecord *rec` — pointer to a record we promise not to modify.
Using `const` here lets the compiler catch any accidental writes and signals
intent to the reader.

`%g` prints a float in the shortest representation — no trailing zeros,
uses scientific notation only when necessary. We reconstruct the original
`.` from our sentinels, which is why the round-trip diff against the
original file is empty.

Diagnostics always go to `stderr`. Actual output goes to `stdout`. This
separation means that `./gtfparse print file.gtf > out.gtf` sends records to
the file while counts still appear in the terminal.

---

## 5. `main.c` — the entry point

### Subcommand dispatch

```c
int main(int argc, char *argv[])
{
    if (argc < 2) { usage(argv[0]); return 1; }

    const char *cmd = argv[1];

    /* compare: two file arguments */
    if (strcmp(cmd, "compare") == 0) { ... }

    /* print: own block — accepts optional --chrom / --source flags */
    if (strcmp(cmd, "print") == 0) {
        /* scan argv for --chrom=..., --source=..., and the file path */
        ...
        int rc = cmd_print(fp, chrom_list, source_list);
        ...
    }

    /* filter, merge, attrs, overlap: one string arg then a file */
    if (strcmp(cmd, "filter")  == 0 || strcmp(cmd, "merge")   == 0 ||
        strcmp(cmd, "attrs")   == 0 || strcmp(cmd, "overlap") == 0) {
        if (argc < 4) { ... }
        FILE *fp = fopen(argv[3], "r");
        int rc;
        if      (strcmp(cmd, "filter")  == 0) rc = cmd_filter(fp,  argv[2]);
        else if (strcmp(cmd, "merge")   == 0) rc = cmd_merge(fp,   argv[2]);
        else if (strcmp(cmd, "attrs")   == 0) rc = cmd_attrs(fp,   argv[2]);
        else                                  rc = cmd_overlap(fp,  argv[2]);
        fclose(fp);
        return rc;
    }

    /* stats: one file argument */
    FILE *fp = fopen(argv[2], "r");
    if (strcmp(cmd, "stats") == 0) ret = cmd_stats(fp);
    else { fprintf(stderr, "unknown command: %s\n", cmd); ... }
}
```

`argc` (argument count) is the number of strings on the command line,
including the program name itself. `argv` (argument vector) is an array of
pointers to those strings:

```
./gtfparse compare a.gtf b.gtf

argv[0] → "./gtfparse\0"
argv[1] → "compare\0"
argv[2] → "a.gtf\0"
argv[3] → "b.gtf\0"
argv[4] → NULL
argc    == 4
```

The dispatch pattern — `strcmp` the first argument against each known
command name — is the standard approach used by tools like `git`, `samtools`,
and `bedtools`. The key advantage over flags (`-p`, `-s`) is that each
subcommand can grow its own options later without polluting a shared flag
namespace.

**Four argument shapes.** Commands fall into four groups checked from
most-specific to least-specific:

| Shape | Commands | `argv` layout |
|-------|----------|---------------|
| Two files | `compare` | `argv[2]`, `argv[3]` are both files |
| Flags + file | `print` | scanned for `--chrom`/`--source`; file is last positional arg |
| String + file | `filter`, `merge`, `attrs`, `overlap` | `argv[2]` is a string, `argv[3]` is a file |
| One file | `stats` | `argv[2]` is a file |

Checking the special cases first is a common C idiom: each block either
returns or falls through to the next, so the final `fopen(argv[2])` is only
reached by commands that take exactly one file argument.

**`print` gets its own block.** The `--chrom` and `--source` options are
optional and combinable, which does not fit the simple `argv[2]`/`argv[3]`
pattern. The `print` block scans all remaining arguments: anything starting
with `--` is a flag; anything else is the file path. `cmd_print` receives
the comma-separated filter lists (or `NULL` for no filter):

```c
static int cmd_print(FILE *fp, const char *chrom_list, const char *source_list)
```

**Sharing a dispatch case.** `filter`, `merge`, `attrs`, and `overlap` all
take the same shape (`<string> <file>`), so they share one `if` block with
a chain of comparisons to pick the right function. The `fopen`/`fclose`
pair appears once regardless of how many commands share the shape.

### Splitting the commands into static functions

```c
static int cmd_print(FILE *fp) { ... }
static int cmd_stats(FILE *fp) { ... }
```

Each command is its own function. This keeps `main` short and readable —
it only does dispatch. The `static` keyword here means these functions are
not visible outside of `main.c`. They are implementation details; nothing
else needs to call them.

### File I/O

```c
FILE *fp = fopen(path, "r");
if (!fp) {
    perror(path);
    return 1;
}
```

`FILE` is an opaque struct managed by the C standard library. You never
access its fields directly — you pass the pointer to `fgets`, `fclose`, etc.
`fopen` returns `NULL` on failure, which is why we check immediately after
every call. `perror` prints the filename plus the OS error message (from the
global `errno` variable that `fopen` set on failure).

```c
while (fgets(line, sizeof(line), fp)) {
```

`fgets(buf, size, stream)` reads at most `size - 1` bytes from `stream`,
stops at a newline (including it in the buffer), and always null-terminates.
It returns `buf` on success, `NULL` at end-of-file or error. Using the
return value directly as the loop condition is the standard C idiom for
line-by-line reading.

### Return codes

```c
return errors > 0 ? 1 : 0;
```

By convention, a Unix program returns `0` on success and a non-zero value
on failure. Shell scripts check `$?` after running a command. Returning `1`
when there were parse errors lets callers detect failure:

```sh
./gtfparse print bad.gtf > /dev/null || echo "parse failed"
```

---

## 6. `stats.h` and `stats.c` — the analytics layer

This section introduces several new C concepts: `static` functions, heap
memory allocation, `qsort` and function pointers, and pointer-walking string
parsing.

### The stats struct

```c
typedef struct {
    long   n_genes;
    long   n_transcripts;
    long   n_exons;
    long   n_introns;
    double mean_exons_per_gene;
    double mean_introns_per_gene;
    long   n_single_exon_genes;
    double pct_single_exon_genes;
    long   n_single_exon_transcripts;
    double pct_single_exon_transcripts;
} GtfExonIntronStats;
```

`double` (64-bit floating point) is preferred over `float` (32-bit) for
computed statistics. The extra precision costs very little — doubles are the
natural floating-point type in C and most library functions use them by
default.

### `static` — file-scope visibility

```c
static int cmp_by_tx(const void *a, const void *b) { ... }
static int cmp_by_gene(const void *a, const void *b) { ... }
static int attr_get(const char *attrs, const char *key,
                    char *val, size_t val_sz) { ... }
```

`static` on a function declaration limits the function's **linkage** to the
current translation unit (the `.c` file). The linker cannot see it; no other
`.c` file can call it. This is the C mechanism for "private" functions.

Without `static`, all three functions would be globally visible. That creates
two problems: naming collisions (another file could define its own `attr_get`)
and a misleading public API (callers might depend on implementation details
that you later want to change). As a rule of thumb, any function that is not
declared in the header should be `static`.

### `attr_get` — walking a string with a pointer

GTF attributes look like this:

```
gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "3";
```

We need to extract a value by key name. The function advances a pointer `p`
through the attribute string one character at a time:

```c
static int attr_get(const char *attrs, const char *key,
                    char *val, size_t val_sz)
{
    size_t klen = strlen(key);
    const char *p = attrs;

    while (*p) {
        while (*p == ' ' || *p == '\t') p++;

        if (strncmp(p, key, klen) == 0 &&
            (p[klen] == ' ' || p[klen] == '\t')) {
            p += klen;
            while (*p == ' ' || *p == '\t') p++;

            if (*p == '"') {
                p++;  /* skip opening quote */
                size_t i = 0;
                while (*p && *p != '"' && i < val_sz - 1)
                    val[i++] = *p++;
                val[i] = '\0';
            }
            return 1;
        }

        while (*p && *p != ';') p++;   /* skip to next semicolon */
        if (*p == ';') p++;
    }
    return 0;
}
```

**`*p` — the current character.** `p` is a `const char *`: a pointer to a
character. `*p` dereferences it to read the byte at that address. The loop
condition `while (*p)` continues until `*p` is `'\0'` (the null terminator
at the end of the string). This is the standard idiom for iterating over a
C string without knowing its length up front.

**`p++` — advancing the pointer.** Incrementing a pointer moves it forward
by `sizeof(*p)` bytes — one byte for `char *`. `p++` is post-increment: the
expression evaluates to the current value of `p` (used here for the copy
`val[i++] = *p++`), then `p` is incremented. So both `i` and `p` advance
together, copying one character per iteration.

**`strncmp` vs `strcmp`.** `strcmp(a, b)` compares entire null-terminated
strings. Here we want to check if the string at position `p` *starts with*
`key`, not that it equals it. `strncmp(p, key, klen)` compares exactly
`klen` bytes, stopping early. After confirming the key matches, we check
`p[klen]` to ensure the match is a complete word (the character after the
key must be whitespace, not, say, the start of a longer key name like
`gene_id_version`).

**`size_t val_sz` — the size parameter.** This function writes into a
caller-provided buffer `val`. Passing the buffer size as a parameter is the
standard safe pattern in C: we check `i < val_sz - 1` before every write,
so we never write past the end. The `- 1` reserves space for the null
terminator we write after the loop.

### Dynamic arrays — `malloc`, `realloc`, `free`

The parser encounters exon records one line at a time and does not know in
advance how many there are. We need a container that can grow.

```c
#define INIT_CAP 4096

size_t cap = INIT_CAP;
TxEntry *rows = malloc(cap * sizeof(TxEntry));
```

`malloc(n)` allocates `n` bytes on the **heap** and returns a `void *`
pointer to the start. The heap is a large pool of memory managed by the
runtime; allocations on it persist until explicitly freed, regardless of
which function called `malloc`. Unlike stack arrays, the size does not need
to be a compile-time constant.

`cap * sizeof(TxEntry)` is the number of bytes needed to hold `cap` elements
of type `TxEntry`. `sizeof` an struct gives its total byte size including any
internal padding the compiler adds for alignment.

`void *` is a "pointer to anything." C allows implicit conversion from
`void *` to any other pointer type, so no cast is needed when assigning to
`TxEntry *rows`.

**Growing the array:**

```c
if (n_rows >= cap) {
    cap *= 2;
    TxEntry *tmp = realloc(rows, cap * sizeof(TxEntry));
    if (!tmp) { free(rows); return -1; }
    rows = tmp;
}
```

`realloc(ptr, new_size)` requests a resized allocation. It may:
- Extend the existing block in place and return the same pointer.
- Allocate a new block, copy the old data, free the old block, and return
  the new pointer.
- Fail and return `NULL`, leaving the original allocation untouched.

Because `realloc` might return a different address, we assign to a temporary
`tmp` first. If we wrote `rows = realloc(rows, ...)` and `realloc` returned
`NULL`, we would overwrite `rows` with `NULL` and lose our only pointer to
the original allocation — a memory leak. The temporary pointer avoids this.

Doubling `cap` on each resize is the standard strategy. An element that
triggers the n-th resize costs O(n) to copy, but that cost is amortized over
the 2^(n-1) elements added since the last resize, giving O(1) amortized cost
per insertion. This is why `std::vector` in C++ uses the same strategy.

**Freeing memory:**

```c
free(rows);
```

Every `malloc` must be matched by exactly one `free`. Failing to free is a
memory leak. Freeing twice ("double free") corrupts the allocator's internal
state and typically crashes the program or creates a security vulnerability.

Note the error-handling paths: every early `return -1` in `gtf_exon_intron_stats`
is preceded by `free(rows)` (or `free(txs)`). If you forget to free on an
error path, the leak only happens when things go wrong — which makes it easy
to miss in testing.

### `qsort` and function pointers

After collecting all exon rows, we sort them to group entries by the same
transcript together:

```c
qsort(rows, n_rows, sizeof(TxEntry), cmp_by_tx);
```

`qsort` is the C standard library's general-purpose sort. Its signature is:

```c
void qsort(void *base, size_t nmemb, size_t size,
           int (*compar)(const void *, const void *));
```

- `base` — pointer to the start of the array.
- `nmemb` — number of elements.
- `size` — byte size of each element (`sizeof(TxEntry)`).
- `compar` — a **function pointer**: the address of a function that `qsort`
  calls to compare any two elements.

`int (*compar)(const void *, const void *)` is the type of a function pointer.
Reading it inside-out: `compar` is a pointer (`*compar`) to a function that
takes two `const void *` arguments and returns `int`. The `&` before the
function name is optional — in C, using a function name in an expression
without calling it automatically decays to a pointer to that function.

**The comparator:**

```c
static int cmp_by_tx(const void *a, const void *b)
{
    return strcmp(((const TxEntry *)a)->tx_id,
                  ((const TxEntry *)b)->tx_id);
}
```

`qsort` passes element addresses as `void *` because it doesn't know the
concrete type — it was designed before C had generics. The comparator must
cast these pointers to the concrete type before using them.

`(const TxEntry *)a` is a **cast**: tell the compiler "treat this `void *`
as a pointer to `const TxEntry`." After the cast, `->tx_id` accesses the
field normally.

The comparator's contract: return a negative number if `a < b`, zero if
`a == b`, positive if `a > b`. `strcmp` already follows this contract for
strings, so we return its result directly.

**Multi-level sort:** `cmp_by_gene` sorts first by `gene_id`, then by
`tx_id` as a tiebreaker. This ensures entries for the same gene are
contiguous, and within a gene, entries are in a consistent order.

```c
static int cmp_by_gene(const void *a, const void *b)
{
    int c = strcmp(((const TxEntry *)a)->gene_id,
                   ((const TxEntry *)b)->gene_id);
    if (c) return c;
    return strcmp(((const TxEntry *)a)->tx_id,
                  ((const TxEntry *)b)->tx_id);
}
```

### The two-pass group-by algorithm

The core of `gtf_exon_intron_stats` is a two-pass approach:

**Pass 1:** collect one row per exon, sort by `tx_id`, merge into one row
per transcript with the exon count.

**Pass 2:** sort the per-transcript rows by `gene_id`, walk through them
gene by gene.

The merge loop after pass 1:

```c
for (size_t i = 0; i < n_rows; ) {
    size_t j = i + 1;
    while (j < n_rows &&
           strcmp(rows[j].tx_id, rows[i].tx_id) == 0) j++;

    txs[n_tx]            = rows[i];
    txs[n_tx].exon_count = (long)(j - i);
    n_tx++;
    i = j;
}
```

After sorting, all rows with the same `tx_id` are adjacent. The inner
`while` advances `j` to the first row with a *different* `tx_id`. The
distance `j - i` is the number of exon rows for this transcript — its
exon count.

`txs[n_tx] = rows[i]` copies the entire struct (gene_id, tx_id, and all
fields) from `rows[i]` into `txs[n_tx]`. C struct assignment copies all
fields by value. We then overwrite `exon_count` with the merged value.

`(long)(j - i)` — `j` and `i` are both `size_t` (unsigned). Their
difference is also `size_t`. The cast to `long` is needed because
`exon_count` is declared `long`. On 64-bit systems both types are 8 bytes,
but mixing signed and unsigned without a cast produces a compiler warning.

The gene-level pass is the same structure: `j` scans forward while `gene_id`
matches, then we process the group `[i, j)` as a single gene before moving
`i = j` to the next group.

### Zeroing a struct with `memset`

```c
memset(stats, 0, sizeof(*stats));
```

`memset(ptr, byte, n)` fills `n` bytes starting at `ptr` with `byte`. Here
we fill every byte of the `GtfExonIntronStats` struct with zero. For integer
fields, all-zero bytes means 0. For `double` fields, all-zero bytes is 0.0
(IEEE 754 guarantees this). This is the standard way to zero-initialize a
struct when you cannot use a literal.

`sizeof(*stats)` — note the dereference. `stats` is a pointer, so
`sizeof(stats)` would give 8 (the pointer size). `sizeof(*stats)` dereferences
the pointer type and gives the size of the struct it points to. `sizeof`
does not actually dereference the pointer at runtime — it only inspects the
type — so this is safe even before the struct contains valid data.

---

## 7. `compare.h` and `compare.c` — the comparison layer

This section introduces the largest collection of new C concepts in the
project: structs that own heap memory, cascading frees, sorting a pointer
array, the sweep line algorithm, and the linear merge for counting shared
elements. It also explains the biological reasoning behind the design.

### Why coordinate-based comparison?

`compare_gtf.py` compares by ID: if both files have a feature named
`ENSG00000001`, they match. That is the right approach when comparing two
versions of the same annotation database.

`gtfparse compare` solves a different problem: two distinct gene predictors
will assign completely different IDs to the same genomic feature. The only
shared currency is position. To decide whether transcript A from predictor 1
corresponds to transcript B from predictor 2, we compare their **intron
chains** — the list of (start, end) coordinates of the gaps between exons.

Two transcripts with identical intron chains use exactly the same splice
sites, which means they represent the same underlying gene structure
regardless of what name each tool assigned. This is the standard criterion
used by tools like GFFCompare.

The four classification codes follow from this logic:

| Code | Condition |
|------|-----------|
| `=`  | All introns match (count and coordinates identical) |
| `j`  | At least one shared intron; chains differ (novel isoform candidate) |
| `o`  | Positional overlap, zero shared introns |
| `u`  | No overlap with any transcript in the other file |

### A struct that owns heap memory

`GtfRecord` (section 3) stores all its data in fixed-size arrays embedded
inside the struct — no pointers, no heap. `Transcript` is different:

```c
typedef struct {
    char        seqname[GTF_FIELD_MAX];
    char        tx_id[GTF_FIELD_MAX];
    char        gene_id[GTF_FIELD_MAX];
    char        strand;
    long        start;
    long        end;
    int         n_exons;
    int         n_introns;
    IntronSpan *introns;     /* heap array [n_introns] */
    int         src;
    char        cmp_class;
} Transcript;
```

`IntronSpan *introns` is a pointer to a heap-allocated array. The struct
itself does not contain the intron data — it contains only the address of
where that data lives. This means `sizeof(Transcript)` is the same regardless
of how many introns the transcript has.

This design creates a two-level ownership structure:

```
Transcript array (heap)
├── txs[0].introns ──→ IntronSpan[n_introns_0]  (heap)
├── txs[1].introns ──→ IntronSpan[n_introns_1]  (heap)
├── txs[2].introns ──→ NULL  (single-exon, no introns)
└── txs[3].introns ──→ IntronSpan[n_introns_3]  (heap)
```

The consequence: freeing the transcript array alone is not enough. Every
non-NULL `introns` pointer must be freed first, then the array itself.
`cmp_free_transcripts` handles this:

```c
void cmp_free_transcripts(Transcript *txs, size_t n)
{
    if (!txs) return;
    for (size_t i = 0; i < n; i++) free(txs[i].introns);
    free(txs);
}
```

Forgetting the inner loop is a classic memory leak — `free(txs)` would
release the array of Transcript structs but leave every intron array
orphaned on the heap with no pointer left to free them.

### Building transcripts from exon rows

`cmp_build_transcripts` uses the same two-phase sort-and-group pattern
introduced in `stats.c`, but adds a third step: computing introns.

**Phase 1** collects raw exon records into an `ExonRow` dynamic array.

**Phase 2** sorts by `(tx_id, start)` and groups into `Transcript` structs.
Within each group of exons for the same transcript, exons are now in
start-order, so introns can be computed as the gaps between consecutive
exons:

```c
for (int k = 0; k < t->n_introns; k++) {
    t->introns[k].start = exons[i + k].end   + 1;
    t->introns[k].end   = exons[i + k + 1].start - 1;
}
```

GTF coordinates are 1-based closed: an exon ending at position 1200 and
the next exon starting at 1500 implies an intron spanning positions
1201–1499. Adding 1 and subtracting 1 converts from the exon boundaries
to the intron boundaries within the same coordinate system.

**Error paths with nested allocations** are tricky. If `malloc` for
`t->introns` fails midway through building the transcript array, we must
free all previously allocated intron arrays before returning:

```c
t->introns = malloc(t->n_introns * sizeof(IntronSpan));
if (!t->introns) {
    for (size_t k = 0; k < tn; k++) free(txs[k].introns);
    free(txs); free(exons);
    return NULL;
}
```

This is the correct pattern: clean up everything you own before returning
an error. If you return NULL without freeing, the caller has no way to
release that memory.

### Sort-then-deduplicate

To count shared introns globally (not just within a locus), we need the set
of unique introns in each file. C has no built-in set type. The standard
idiom is **sort then scan**:

```c
qsort(arr, total, sizeof(IntronKey), intron_key_cmp);

long uniq = 0;
for (long x = 0; x < total; x++) {
    if (x == 0 || intron_key_cmp(&arr[x], &arr[x - 1]) != 0)
        arr[uniq++] = arr[x];
}
```

After sorting, duplicates are adjacent. A single forward scan copies only
elements that differ from their predecessor. The result is a deduplicated
array in the same memory. `uniq` is the number of unique elements.

This pattern replaces a hash set. It costs O(n log n) for the sort, whereas
a hash set would give O(n) expected. But it needs no data structure
implementation — just `qsort` and a loop — which is why it is the dominant
idiom in C genomics code.

### Sorting a pointer array

After sorting each file's transcripts in place, `cmp_compare` builds a
pointer array over all of them:

```c
Transcript **all = malloc(ntotal * sizeof(Transcript *));
for (size_t i = 0; i < na; i++) all[i]      = &a[i];
for (size_t i = 0; i < nb; i++) all[na + i] = &b[i];
qsort(all, ntotal, sizeof(Transcript *), tx_ptr_cmp);
```

`all` is an array of *pointers to* Transcript. Each element is 8 bytes
(the size of a pointer on a 64-bit system). Sorting this pointer array
instead of copying the Transcript structs themselves matters: a Transcript
is several hundred bytes because of its fixed-size string fields. Copying
tens of thousands of them during a sort would be expensive. Sorting 8-byte
pointers is much cheaper.

The comparator for a pointer array requires one extra level of indirection:

```c
static int tx_ptr_cmp(const void *a, const void *b)
{
    const Transcript *ta = *(const Transcript * const *)a;
    const Transcript *tb = *(const Transcript * const *)b;
    // ... compare ta and tb
}
```

`qsort` passes the *address of each element* as `void *`. Each element in
`all` is a `Transcript *`, so each element's address is a `Transcript **`.
The cast `(const Transcript * const *)a` tells the compiler "treat this
`void *` as a pointer to a `const Transcript *`". Dereferencing it with `*`
gives the `Transcript *` that was stored in that slot of `all`.

Reading the type inside-out: `const Transcript * const *` is a pointer to
a const pointer to a const Transcript. The outer `const` (between `*` and
`*`) prevents the comparator from modifying the pointer stored in the array.
The inner `const` (before `Transcript`) prevents it from modifying the
struct itself. Both are good defensive practice for a comparator.

Contrast with the comparators in `stats.c`:

```c
// stats.c: sorting an array of TxEntry structs directly
static int cmp_by_tx(const void *a, const void *b)
{
    const TxEntry *ea = (const TxEntry *)a;   // one level of indirection
    ...
}

// compare.c: sorting an array of Transcript * pointers
static int tx_ptr_cmp(const void *a, const void *b)
{
    const Transcript *ta = *(const Transcript * const *)a;  // two levels
    ...
}
```

The extra dereference (`*`) is the only difference. One level of
indirection when the array contains structs; two levels when it contains
pointers.

### The sweep line algorithm

The sweep line is the fundamental algorithm for grouping overlapping
genomic intervals. It underlies virtually all coordinate-aware genomics
tools: bedtools intersect, samtools sort, STAR's junction detection.

The idea: **if you sort intervals by start position, overlapping intervals
are always adjacent.** A single forward scan can find all connected groups
(loci) in O(n) after the O(n log n) sort.

```c
size_t lo = 0;
while (lo < ntotal) {
    const char *lseq    = all[lo]->seqname;
    char        lstrand = all[lo]->strand;
    long        lend    = all[lo]->end;
    size_t hi = lo + 1;

    while (hi < ntotal
           && strcmp(all[hi]->seqname, lseq)  == 0
           && all[hi]->strand == lstrand
           && all[hi]->start  <= lend) {
        if (all[hi]->end > lend) lend = all[hi]->end;
        hi++;
    }

    compare_locus(all + lo, hi - lo, stats);
    lo = hi;
}
```

`lo` is the index of the first transcript in the current locus. `lend`
tracks the rightmost end seen so far across all transcripts added to the
locus. The inner loop extends the locus as long as the next transcript's
start is within `lend` (i.e. it overlaps with *something* already in the
locus). When it does not, the locus is closed and `compare_locus` is called.

Two conditions break the locus regardless of position: a different `seqname`
(different chromosome) or a different `strand`. A transcript on chr1+ cannot
overlap one on chr1−. Sorting by `(seqname, strand, start)` ensures these
boundaries are respected automatically — when `seqname` or `strand` changes,
the `strcmp`/`!=` checks fail and the while loop stops.

Note that `lend` must be updated with `max(lend, all[hi]->end)` as each
transcript is added. Without this, a long transcript that extends further
right than the next one would fail to capture transcripts that start after
the next transcript's end but still within the long one's span.

### Linear merge for intron chain comparison

Within a locus, for each pair of transcripts from different files,
`count_shared_introns` counts how many introns are identical:

```c
static int count_shared_introns(const Transcript *a, const Transcript *b)
{
    int i = 0, j = 0, shared = 0;
    while (i < a->n_introns && j < b->n_introns) {
        long as = a->introns[i].start, ae = a->introns[i].end;
        long bs = b->introns[j].start, be = b->introns[j].end;
        if (as == bs && ae == be) {
            shared++; i++; j++;
        } else if (as < bs || (as == bs && ae < be)) {
            i++;
        } else {
            j++;
        }
    }
    return shared;
}
```

This is the **merge scan** idiom, the same pattern used to count shared
introns globally in `cmp_compare`. It works because both intron arrays are
already sorted by start position (they were built from exons sorted by start
in `cmp_build_transcripts`).

The rule is: always advance the pointer into whichever array currently has
the smaller element. If both are equal, count the match and advance both.
This ensures every element is visited exactly once — O(m + n) total — rather
than the O(m × n) of a naive nested loop.

Visualised on a small example:

```
A introns: [1201-1499]  [1701-1999]
B introns: [1201-1499]  [2001-2399]

Step 1: A[0]==B[0] → shared=1, i=1, j=1
Step 2: A[1].start(1701) < B[1].start(2001) → advance i
Step 3: i==2, loop ends. shared=1 → class 'j'
```

### Classifying transcripts

Once `count_shared_introns` returns, the classification logic is:

```c
if (t->n_introns > 0
        && shared == t->n_introns
        && shared == o->n_introns)
    cls = '=';
else if (shared > 0)
    cls = 'j';
else
    cls = 'o';
```

`=` requires that *all* introns on *both* sides match. If A has 2 introns
and B has 3, but the 2 match, that is `j` (a partial match), not `=`. The
extra intron in B means B has an additional exon that A does not — the
structures differ even though they share some splice sites.

The `n_introns > 0` guard ensures that two single-exon transcripts with
positional overlap never receive `=`: with no introns there is nothing to
match, so `o` is the correct class.

`class_rank` assigns a numeric priority (`=`=3, `j`=2, `o`=1, `u`=0) so
that the best class across all opponents is kept with a simple `if
(class_rank(cls) > class_rank(best))` comparison. This pattern — encoding
an ordered preference as integers — is cleaner than a chain of if/else
comparisons and easier to extend.

---

## 8. Pointers — a complete picture

A pointer is a variable that holds a **memory address**. The type of the
pointer tells the compiler how to interpret the bytes at that address.

```
  int x = 42;          int *p = &x;

  Address  Value        Address  Value
  0x1000   42    ◀───   0x2000   0x1000
    x                     p
```

`&x` is the address-of operator — it gives you the address of `x`.
`*p` is the dereference operator — it follows the pointer and gives you the
value at that address. `p->field` is shorthand for `(*p).field`.

### The five pointers in this program

**1. `const char *line`** — read-only view of caller's data

The function receives the address of a string it may read but not modify.
`const` on the data (not the pointer) enforces this at compile time.

**2. `GtfRecord *rec`** — output parameter

The caller passes `&rec` (address of its local struct). The function writes
results directly into the caller's memory. This is how C achieves
"multiple return values."

**3. `char *fields[9]`** — array of pointers into `buf`

After `strtok` runs, each `fields[i]` points into a different position in
`buf`. No new memory is allocated; this is pointer arithmetic into an
existing buffer.

**4. `char *tok`** — strtok's return value

`strtok` returns a pointer into the string it was given, or `NULL` when
exhausted. Checking `while (tok && ...)` exploits the fact that `NULL`
(address zero) is falsy.

**5. `FILE *fp`** — opaque library handle

`FILE` is a struct you never look inside. The library owns it; you own only
the pointer. Always pass it back to `fclose` when done.

**6. `TxEntry *rows`** — heap-allocated array

`malloc` returns a pointer to newly allocated heap memory. Unlike the five
pointers above, which point into stack memory or library-managed memory,
this pointer represents an allocation you own: you must call `free(rows)`
when you are done with it.

**7. `Transcript **all`** — array of pointers to structs

A pointer to a pointer. Each element of `all` is itself a pointer to a
`Transcript`. This indirection is what allows sorting by pointer (cheap,
8 bytes per element) without copying the full struct (expensive, hundreds
of bytes). Section 7 explains this in detail.

---

## 9. Memory lifetimes and dangling pointers

Every piece of memory in C has a **lifetime** — a period during which it is
valid to read and write through pointers to it.

### Stack memory

Local variables live on the **stack** — a region of memory the CPU manages
automatically. When a function is called, a stack frame is pushed. When the
function returns, the frame is popped and that memory is reclaimed.

```c
int gtf_parse_line(const char *line, GtfRecord *rec)
{
    char buf[65536];      // lives on gtf_parse_line's stack frame
    char *fields[9];      // also on this stack frame

    // fields[0..8] point INTO buf — valid here

    return 0;
    // stack frame is destroyed — buf and fields are gone
    // any pointer that still points into buf is now DANGLING
}
```

A **dangling pointer** is a pointer whose target has been destroyed. Reading
through it is undefined behaviour: you might get stale data, garbage, or a
crash. The compiler is not required to warn you.

### Why our code is safe

We copy all data out of `buf` before the function returns:

```c
strncpy(rec->seqname, fields[0], GTF_FIELD_MAX - 1);
```

`fields[0]` points into `buf`. We copy the bytes into `rec->seqname`, which
lives in `main`'s stack frame and outlives the call. By the time `buf`
disappears, we have everything we need in the struct.

### What would be broken

If we returned a pointer into `buf` instead of copying, the caller would
hold a dangling pointer:

```c
// BROKEN — never do this
char *broken_get_feature(const char *line) {
    char buf[65536];
    strncpy(buf, line, sizeof(buf) - 1);
    char *tok = strtok(buf, "\t");
    tok = strtok(NULL, "\t");
    tok = strtok(NULL, "\t");
    return tok;  // points into buf, which is about to be destroyed
}
```

### Heap memory

`malloc(n)` allocates `n` bytes on the **heap** — memory that persists until
explicitly freed with `free()`, regardless of which function allocated it.
This removes the stack lifetime problem.

`stats.c` already uses heap memory: `rows` and `txs` are heap arrays. They
are allocated in `gtf_exon_intron_stats`, built up over the course of the
function, and freed before the function returns. The data we care about
(the computed statistics) is written into the caller's `GtfExonIntronStats`
struct — which lives on the caller's stack — before we free the heap arrays.

The pattern is the same as with `buf` in the parser, just with a longer
lifetime: heap memory lives until `free()`, not until the function returns.

### The three ownership rules

These rules prevent the two most common memory bugs (leaks and dangling pointers):

1. **Every `malloc` must have exactly one `free`.**
2. **Never use a pointer after the memory it points to has been freed.**
3. **Never return a pointer to a local variable.**

`stats.c` follows all three: every allocation is freed on every code path
(including error paths), freed pointers are never used again, and no local
addresses are returned.

---

## 10. Current limitations

This section describes the known constraints of the current code: what has
been addressed and what remains.

### Buffer size caps — partially addressed

**Stored records** — addressed in §18. `GtfRecord` now stores
`seqname`/`source`/`feature` as interned `const char *` pointers and
`attributes` as a `strndup`'d `char *`. There is no fixed cap on what
the table holds; the exact string is preserved.

**Parsing** — still present. `GtfRawRecord` (the scratch buffer used by
`gtf_parse_line`) still has `char seqname[256]` etc. A value longer than
255 characters is silently truncated *before* it reaches the pool. Removing
this last cap would require replacing `fgets` + `strtok` with `getline`
and scanning the raw buffer without a fixed intermediate copy.

### Line length cap — still present

`buf[65536]` in `gtf_parse_line` limits lines to 65,535 characters. GTF
attributes columns can be long on heavily annotated genomes, but 64 KB
is generous in practice. The POSIX `getline()` function allocates a buffer
large enough for each line dynamically and would remove this limit.

### `strtok` is not thread-safe — still present

`strtok` stores its position in a hidden static variable. Calling it from
two threads simultaneously corrupts both. `strtok_r` (POSIX) takes an
explicit `char **saveptr` instead and is safe. This matters only if the
parser is ever called from multiple threads.

### Dynamic record array — addressed in Step 2

`GtfTable` (see section 11) provides a heap-allocated, doubling-growth
array that holds all records from a file in memory. Commands that need
to make multiple passes over the data (`filter`, `attrs`) load into a
`GtfTable` rather than streaming line-by-line.

### Attribute parsing — addressed in Step 4

`stats.c` had a private `attr_get` function that scanned the raw attribute
string character by character each time it needed a value. `GtfAttrs`
(see section 13) provides a public API in the `gtf.c` layer: parse once,
look up by key, free when done. The `attrs` command demonstrates it.

### GFF3 support — addressed in Step 6

`gtf_table_load` now detects the `##gff-version` pragma and sets
`t->format = GTF_FMT_GFF3`. It also stops at `##FASTA` rather than trying
to parse sequence data as annotation records. `gff3_attrs_parse` handles
the `key=value;` attribute syntax. See section 14.

### Sorting and interval queries — addressed in Steps 5, 7, and further refined

`gtf_table_sort` sorts the table in place by `(seqname, start, end)` and
then builds a hash index (§19) mapping each seqname to its first sorted
index. `gtf_table_query` uses the hash index for O(1) chromosome lookup,
then scans forward using pointer comparison for seqname-block detection,
giving O(1 + k) per query after sorting. See §15, §16, §19.

### `stats.c` loads all exon records into memory

`gtf_exon_intron_stats` reads all exon lines before computing anything.
For a genome with tens of millions of exons this may require several hundred
megabytes. A streaming algorithm that processes records in gene-sorted order
could reduce peak memory to O(transcripts per gene), but it would require
the input to be sorted — a requirement the current code does not impose.

---

## 11. `GtfTable` — a dynamic array of records

### Why load everything into memory?

The original `cmd_print` streamed the file one line at a time: parse a
record, print it, discard it, repeat. Streaming uses O(1) memory but limits
you to a single forward pass. Commands that need to look at all records
more than once — filtering, sorting, interval queries — require the whole
file to be in memory first.

`GtfTable` is that in-memory representation. It is the same doubling-growth
dynamic array pattern introduced in `stats.c` (section 6), but generalised
to hold `GtfRecord` values instead of private `TxEntry` structs.

### The struct

```c
typedef struct {
    GtfRecord *recs;
    size_t     n;
    size_t     cap;
    size_t     n_skipped;
    size_t     n_errors;
    GtfFormat  format;   /* GTF_FMT_GTF or GTF_FMT_GFF3 */
    StrPool    pool;     /* owns all interned seqname/source/feature strings */
    HTable    *index;    /* seqname → first sorted index; NULL until sorted */
} GtfTable;
```

- `recs` — pointer to the heap-allocated array of `GtfRecord` values.
- `n` — the number of valid records currently stored.
- `cap` — the allocated capacity (always `>= n`).
- `n_skipped` — count of comment and blank lines while loading.
- `n_errors` — count of lines that failed to parse.
- `format` — detected file format (`GTF_FMT_GTF` or `GTF_FMT_GFF3`).
- `pool` — the string intern pool that owns every `seqname`, `source`, and
  `feature` string stored in `recs`. See §18.
- `index` — hash table built by `gtf_table_sort`; maps each seqname to the
  first record index on that chromosome. `NULL` until `gtf_table_sort` is
  called. See §19.

The diagnostic fields are written to `stderr` by commands that care
(`cmd_print`) and ignored by commands that do not.

### `gtf_table_load`

```c
GtfTable *gtf_table_load(FILE *fp)
{
    GtfTable *t = malloc(sizeof(GtfTable));
    // ... initialise fields including pool and index = NULL ...
    t->recs = malloc(t->cap * sizeof(GtfRecord));

    char line[65536];
    GtfRawRecord raw;   /* stack scratch buffer, re-used each iteration */

    while (fgets(line, sizeof(line), fp)) {
        int ret = gtf_parse_line(line, &raw);
        if (ret == 1) { t->n_skipped++; continue; }
        if (ret == -1) { t->n_errors++;  continue; }

        // grow recs array if needed (same doubling pattern as before) ...

        GtfRecord *r  = &t->recs[t->n];
        r->seqname    = str_pool_intern(&t->pool, raw.seqname);
        r->source     = str_pool_intern(&t->pool, raw.source);
        r->feature    = str_pool_intern(&t->pool, raw.feature);
        r->start      = raw.start;
        r->end        = raw.end;
        r->score      = raw.score;
        r->strand     = raw.strand;
        r->frame      = raw.frame;
        r->attributes = strndup_safe(raw.attributes, strlen(raw.attributes));
        t->n++;
    }
    return t;
}
```

**`GtfRawRecord raw` — a stack scratch buffer.** Every iteration overwrites
the same stack memory. `gtf_parse_line` fills the fixed arrays, then we
convert into the table slot by calling `str_pool_intern` (which either
returns an existing pooled pointer or allocates a new one) and
`strndup_safe` (which copies the exact attribute string to the heap).

**Three-way return from `gtf_parse_line`.** Returns `1` for comments/blanks
(increment `n_skipped`), `-1` for errors (increment `n_errors`), and `0`
for a valid record. Only `0` causes a record to be added.

**Grow before insert.** The capacity check comes *before* storing the new
record, maintaining the invariant `t->n < t->cap` at the point of assignment.

**`realloc` via a temporary.** Assigning `realloc`'s return value directly
to `t->recs` would leak the old allocation if `realloc` returned `NULL`.
The temporary `tmp` avoids this.

**Field-by-field fill, not struct assignment.** Unlike the old
`t->recs[t->n++] = rec`, the new code fills one field at a time using
`str_pool_intern` and `strndup_safe`. The struct cannot be bulk-copied
because each pointer field needs its own allocation decision.

**Failure path calls `gtf_table_free`.** If any allocation fails mid-load,
`gtf_table_free(t)` frees the partially filled table. The free function
handles partial initialisation safely — see below.

### `gtf_table_free`

```c
void gtf_table_free(GtfTable *t)
{
    if (!t) return;
    for (size_t i = 0; i < t->n; i++)
        free(t->recs[i].attributes);    /* each record owns its attributes */
    free(t->recs);
    for (size_t i = 0; i < t->pool.n; i++)
        free(t->pool.data[i]);          /* pool owns the interned strings */
    free(t->pool.data);
    htable_free(t->index);              /* frees slots and HTable struct */
    free(t);
}
```

After string interning, `GtfRecord` is no longer a plain value type — it
contains `char *attributes` (heap-allocated per record) and `const char *`
pointer fields (pool-owned). The free cascade must handle three ownership
layers before freeing the container:

1. **Per-record**: free each `attributes` string individually.
2. **Pool**: free each interned string in `pool.data`, then `pool.data` itself.
3. **Index**: `htable_free` frees the slot array and `HTable` struct.
   It does *not* free the key strings — those are in the pool (layer 2).

Compare with the original one-liner `free(t->recs); free(t)`: that worked
because the old `GtfRecord` had only fixed-size embedded arrays and no
separately allocated memory. The switch to pointer fields trades simplicity
of freeing for a 60× reduction in memory usage.

Compare also with `Transcript` in `compare.c`, which has `IntronSpan *introns`
in each element — it also requires a loop to free inner allocations before
the outer array. The pattern is the same: owner frees contents, then container.

---

## 12. `filter` — using `GtfTable` in a command

`cmd_filter` is the simplest consumer of `GtfTable`:

```c
static int cmd_filter(FILE *fp, const char *feature)
{
    GtfTable *t = gtf_table_load(fp);
    if (!t) { fprintf(stderr, "error: out of memory\n"); return 1; }

    size_t matched = 0;
    for (size_t i = 0; i < t->n; i++) {
        if (strcmp(t->recs[i].feature, feature) == 0) {
            gtf_print(&t->recs[i]);
            matched++;
        }
    }

    fprintf(stderr, "matched=%zu  total=%zu\n", matched, t->n);
    gtf_table_free(t);
    return 0;
}
```

There is not much new C here — the interest is in the design pattern. The
function does not need to know how many records the file has. It loads once,
iterates once, and prints matching records to `stdout` while reporting
counts to `stderr`.

**Why not stream?** For filtering alone, streaming would work and use less
memory. The reason to use `GtfTable` here is consistency: every command in
this tool works from an in-memory table, which means any command can later
be extended to make a second pass without changing its file-reading logic.
The `overlap` command (section 16) shows this payoff: it sorts the table in
place after loading, then binary-searches it — a step that would be impossible
without the whole file in memory.

**`strcmp` on a struct field.** `t->recs[i].feature` is a `const char *`
— an interned pointer stored directly in the `GtfRecord`. `strcmp` takes
two `const char *` arguments and works identically whether the string came
from a fixed array or a heap allocation. No array-to-pointer decay is
needed here; the field is already a pointer.

(In earlier versions `feature` was `char feature[GTF_FIELD_MAX]`, an
embedded array. In that case, passing `t->recs[i].feature` to `strcmp`
relied on array-to-pointer decay: an array name, in most expression
contexts, becomes a pointer to its first element. The call was syntactically
the same but the underlying type was different.)

---

## 13. `GtfAttrs` — parsing the attributes column

### The problem with raw string scanning

`stats.c` contains a private `attr_get` function that finds a key's value
by walking the raw attribute string character by character. It works, but
it re-scans the whole string for every key you look up. If you need three
attributes from one record, you scan the string three times.

More importantly, `attr_get` is `static` — private to `stats.c`. Any other
caller that needs attribute access must either duplicate the function or
work without it. The `GtfAttrs` API moves attribute parsing into the shared
`gtf.c` layer where any caller can use it.

### The two structs

```c
typedef struct {
    char *key;
    char *value;
} GtfAttr;

typedef struct {
    GtfAttr *pairs;
    size_t   n;
} GtfAttrs;
```

`GtfAttr` holds one key/value pair. Both strings are heap-allocated and
owned by the pair. `GtfAttrs` holds a dynamic array of pairs. The two types
form a small two-level ownership tree, similar to the `Transcript`/`IntronSpan`
relationship in `compare.c`:

```
GtfAttrs (heap)
├── pairs[0].key   → "gene_id\0"        (heap)
├── pairs[0].value → "ENSG00000001\0"   (heap)
├── pairs[1].key   → "transcript_id\0"  (heap)
├── pairs[1].value → "ENST00000001\0"   (heap)
└── ...
```

Freeing a `GtfAttrs` requires three levels: free each key, free each value,
free the pairs array, free the struct.

### `strndup_safe` — copying a substring to the heap

The parser needs to copy substrings out of the raw attributes string into
individual heap allocations. The standard library has `strdup` (duplicate an
entire string) but not `strndup` on all platforms. We define our own:

```c
static char *strndup_safe(const char *src, size_t len)
{
    char *s = malloc(len + 1);
    if (!s) return NULL;
    memcpy(s, src, len);
    s[len] = '\0';
    return s;
}
```

`malloc(len + 1)` allocates exactly enough bytes for the substring plus its
null terminator. `memcpy` copies the bytes (it does not stop at `'\0'`, which
is what we want — we are copying a substring, not a null-terminated string).
The explicit `s[len] = '\0'` adds the terminator.

This is `static` (private to `gtf.c`) because it is an implementation detail
with no reason to be visible outside the file.

### `gtf_attrs_parse` — pointer walking

```c
GtfAttrs *gtf_attrs_parse(const char *raw)
{
    GtfAttrs *a = malloc(sizeof(GtfAttrs));
    // ...
    a->n = 0;
    a->pairs = malloc(GTF_ATTRS_INIT_CAP * sizeof(GtfAttr));  /* 8 pairs */

    const char *p = raw;
    while (*p) {
        while (*p == ' ' || *p == '\t') p++;  /* skip whitespace */
        if (!*p) break;

        /* extract key */
        const char *key_start = p;
        while (*p && *p != ' ' && *p != '\t' && *p != ';') p++;
        char *key = strndup_safe(key_start, p - key_start);

        while (*p == ' ' || *p == '\t') p++;  /* skip whitespace */

        /* extract value — quoted or bare */
        if (*p == '"') {
            p++;  /* skip opening quote */
            const char *val_start = p;
            while (*p && *p != '"') p++;
            value = strndup_safe(val_start, p - val_start);
            if (*p == '"') p++;  /* skip closing quote */
        } else { ... }

        /* store the pair, grow array if needed */
        a->pairs[a->n].key   = key;
        a->pairs[a->n].value = value;
        a->n++;

        while (*p && *p != ';') p++;  /* skip to next semicolon */
        if (*p == ';') p++;
    }
    return a;
}
```

`p` starts at the beginning of the raw string and advances forward — it
never moves backward. This is the same pointer-walking style used in
`attr_get` in `stats.c`, but instead of stopping when one key is found, it
continues until the entire string has been consumed, producing all pairs.

**`p - key_start` — pointer subtraction.** Both `p` and `key_start` point
into the same string. Subtracting two pointers of the same type gives the
number of elements between them — here, the number of characters in the key.
Pointer subtraction yields a value of type `ptrdiff_t` (a signed integer),
but `strndup_safe` takes `size_t` (unsigned); the cast is safe because `p`
is always at or ahead of `key_start`.

**Error paths with partial construction.** If any `malloc` inside the loop
fails, we have a partially built `GtfAttrs` — some pairs have been stored,
some have not, and the just-allocated `key` or `value` that triggered the
failure is dangling. The cleanup order is:
1. Free the just-allocated string that triggered failure (if any).
2. Call `gtf_attrs_free(a)` — it iterates over `a->n` already-stored pairs
   and frees their keys, values, the pairs array, and the struct itself.

This works because `a->n` always reflects only the pairs that have been
*fully* stored (both key and value set). A pair that failed mid-construction
has not yet been added to `a->n`, so `gtf_attrs_free` never sees it.

### `gtf_attrs_get` — linear scan

```c
const char *gtf_attrs_get(const GtfAttrs *a, const char *key)
{
    for (size_t i = 0; i < a->n; i++)
        if (strcmp(a->pairs[i].key, key) == 0)
            return a->pairs[i].value;
    return NULL;
}
```

A linear scan over the pairs array. A typical GTF record has 5–20 attributes,
so this is fast in practice. A hash table would give O(1) lookup, but adds
significant implementation complexity for little practical gain at this scale.

The function returns `const char *` — a read-only pointer to the value
string owned by the `GtfAttrs`. The caller must not free or modify it, and
must not use it after calling `gtf_attrs_free`.

### `gtf_attrs_free` — cascading free

```c
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
```

The `if (!a) return` guard makes the function safe to call with `NULL`,
which is the standard convention for free-like functions in C (`free(NULL)`
is also defined as a no-op). This matters in error paths where `gtf_attrs_parse`
may return `NULL` — the caller can call `gtf_attrs_free(a)` without
checking first.

The loop frees the two heap strings inside each pair. After the loop,
`free(a->pairs)` releases the array of pairs. Finally `free(a)` releases the
`GtfAttrs` struct itself. Reversing this order — freeing the struct or the
pairs array before the strings — would not cause a crash immediately, but
the string allocations would become unreachable: a memory leak.

### Why not embed `GtfAttrs` inside `GtfRecord`?

The simpler design might be to parse attributes automatically in
`gtf_parse_line` and store the result directly in `GtfRecord`:

```c
typedef struct {
    // ... existing fields ...
    char     attributes[GTF_ATTR_MAX];  /* raw, for printing */
    GtfAttr *attrs;                     /* parsed, for lookup */
    size_t   n_attrs;
} GtfRecord;
```

This would be convenient for callers — no separate parse step. But it comes
with a significant cost: `GtfRecord` would no longer be a **value type**.

Currently you can write:

```c
GtfRecord rec;
gtf_parse_line(line, &rec);
// use rec
// no cleanup needed — rec lives on the stack and disappears automatically
```

If `GtfRecord` contained a heap-allocated `attrs` array, every `gtf_parse_line`
call would allocate memory that must eventually be freed. The streaming loops
in `stats.c` and `compare.c` — which declare a single `GtfRecord rec` and
reuse it each iteration — would suddenly leak memory on every line. The
`GtfTable` array would require a cascading free loop over all records, as
`cmp_free_transcripts` does for `Transcript`.

Keeping `GtfAttrs` separate preserves `GtfRecord`'s value-type semantics.
Callers that need attribute access parse on demand and free when done.
Callers that do not (like `stats.c`, which uses its own internal parser)
pay nothing.

### `cmd_attrs` — tying it together

```c
static int cmd_attrs(FILE *fp, const char *key)
{
    GtfTable *t = gtf_table_load(fp);
    if (!t) { ... return 1; }

    for (size_t i = 0; i < t->n; i++) {
        GtfAttrs *a = gtf_attrs_parse(t->recs[i].attributes);
        if (!a) { gtf_table_free(t); return 1; }
        const char *val = gtf_attrs_get(a, key);
        printf("%s\n", val ? val : ".");
        gtf_attrs_free(a);
    }

    gtf_table_free(t);
    return 0;
}
```

`GtfAttrs *a` is created, used, and freed inside a single loop iteration.
It is a **temporary** — its lifetime is shorter than the loop variable `i`.
This is a useful mental model: when you see `malloc` and `free` bracketing
a use, the allocation is scoped to that block even if C does not enforce it.

`val ? val : "."` — the ternary chooses between the found value and the
sentinel string `"."`. Printing `"."` for absent keys (the GTF convention
for missing fields) means the output has exactly one line per input record,
making it safe to `paste` with other per-record output or process with `awk`.

`printf("%s\n", val ? val : ".")` — the `%s` format specifier requires a
non-NULL pointer. If `gtf_attrs_get` returned `NULL` and we printed it
directly, the behaviour would be undefined (and typically a crash on most
platforms). The ternary ensures a valid pointer is always passed.

---

## 14. GFF3 support — a second file format

### GTF vs GFF3

GTF and GFF3 share the same nine-column tab-delimited structure. The
differences that matter for parsing are entirely in the attributes column
and in a handful of pragma lines that appear at the top of GFF3 files.

| Aspect | GTF | GFF3 |
|--------|-----|------|
| Attribute separator | `;` between pairs | `;` between pairs |
| Key/value separator | space: `key "value"` | equals sign: `key=value` |
| Value quoting | always double-quoted | never quoted |
| Pragma lines | none | `##gff-version`, `##FASTA`, etc. |
| Embedded sequences | never | optional `##FASTA` section at end |

Because the nine columns are identical, `gtf_parse_line` works for GFF3
lines without any changes — the raw attribute string is stored as-is
regardless of format.

### The `GtfFormat` enum

```c
typedef enum { GTF_FMT_GTF = 0, GTF_FMT_GFF3 = 1 } GtfFormat;
```

An `enum` (enumeration) assigns names to integer constants. `GTF_FMT_GTF`
has value `0` and `GTF_FMT_GFF3` has value `1`. Using named constants
instead of raw integers makes code self-documenting and lets the compiler
warn if a switch statement omits a case.

The `= 0` initialiser on `GTF_FMT_GTF` means that a zero-initialised
`GtfTable` (e.g. from `calloc` or a `= {0}` initialiser) defaults to GTF
format — the safer assumption when no pragma has been seen.

`GtfTable` gains a `format` field of this type. Every command that parses
attributes (`cmd_attrs`) reads it to choose the right parser.

### Format detection in `gtf_table_load`

```c
while (fgets(line, sizeof(line), fp)) {
    if (strncmp(line, "##FASTA", 7) == 0)
        break;
    if (strncmp(line, "##gff-version", 13) == 0) {
        t->format = GTF_FMT_GFF3;
        t->n_skipped++;
        continue;
    }
    int ret = gtf_parse_line(line, &rec);
    // ...
}
```

Each line is inspected *before* being passed to `gtf_parse_line`. Two
GFF3-specific cases are handled first:

**`##FASTA` — break.** GFF3 files may embed the reference sequences after
a `##FASTA` pragma. Everything after that line is FASTA, not annotation.
`break` exits the loop immediately; `fgets` is never called again for this
file. Without this check, `gtf_parse_line` would try to parse FASTA lines
as nine-column records and report a flood of parse errors.

**`##gff-version` — set format, continue.** The pragma is counted as a
skipped line (it carries no record data) and the format flag is flipped.
All subsequent attribute parsing will use `gff3_attrs_parse`.

`strncmp(line, prefix, len)` compares only the first `len` characters,
so it matches `##gff-version 3`, `##gff-version 3.1.25`, and any other
version string equally well. The full pragma content is irrelevant — the
presence of the prefix is what identifies GFF3.

### `gff3_attrs_parse`

```c
GtfAttrs *gff3_attrs_parse(const char *raw)
{
    // ...
    const char *p = raw;
    while (*p) {
        while (*p == ' ' || *p == '\t') p++;
        if (!*p) break;

        /* locate end of this token at the next ';' */
        const char *tok_end = p;
        while (*tok_end && *tok_end != ';') tok_end++;

        /* locate the '=' separator within the token */
        const char *eq = p;
        while (eq < tok_end && *eq != '=') eq++;

        if (eq == tok_end) {
            /* bare tag with no value — store with empty value */
        } else {
            char *key   = strndup_safe(p,      eq - p);
            char *value = strndup_safe(eq + 1, tok_end - eq - 1);
            // store pair, grow array...
        }

        p = tok_end;
        if (*p == ';') p++;
    }
    return a;
}
```

The key structural difference from `gtf_attrs_parse` is how the token
boundaries are found. GTF parsing walks character-by-character and extracts
the value based on quote delimiters. GFF3 parsing instead locates the entire
token boundary first (`tok_end` at the next `;`), then splits that token at
the `=` sign.

**Two-pointer scanning.** `tok_end` finds the end of the current
`key=value` token. `eq` finds the `=` within it. Both scan forward from
`p` and both are bounded: `tok_end` is bounded by the null terminator,
`eq` is bounded by `tok_end`. Neither can overrun the string.

**Bare tags.** Some GFF3 files include tags without values (e.g.
`Is_circular`). When `eq == tok_end`, no `=` was found. Rather than
skipping the tag silently, it is stored with an empty string value so
callers can test for its presence with `gtf_attrs_get`.

**Multi-value attributes.** GFF3 allows comma-separated values for one
key: `Parent=mRNA:t1,mRNA:t2`. The entire comma-separated string is stored
as a single value; splitting is intentionally left to the caller.

**Percent-encoding.** The GFF3 spec requires certain characters (`;`, `=`,
`%`, `&`, `,`) to be percent-encoded when they appear inside values. This
implementation does not decode them. For files from standard sources
(NCBI, Ensembl) this is not an issue in practice.

`gff3_attrs_parse` returns the same `GtfAttrs` type as `gtf_attrs_parse`.
The same `gtf_attrs_get` and `gtf_attrs_free` functions work on both.

---

## 15. `gtf_table_sort` — sorting records by position

### Why sort?

An unsorted table supports only linear scan: to find all records overlapping
a region you must examine every record. Sorting by position enables binary
search: jump directly to the chromosome of interest, then scan only the
records that could possibly overlap. For a file with a million records and
a query that spans 100 of them, the difference is roughly 1,000,000
comparisons vs ~20 (binary search) + 100 (scan).

Sorting also groups all records from the same chromosome together, which is
a prerequisite for the sweep-line algorithm used in `compare.c`.

### The comparator

```c
static int rec_cmp(const void *a, const void *b)
{
    const GtfRecord *ra = (const GtfRecord *)a;
    const GtfRecord *rb = (const GtfRecord *)b;
    int c = strcmp(ra->seqname, rb->seqname);
    if (c) return c;
    if (ra->start != rb->start) return (ra->start < rb->start) ? -1 : 1;
    return (ra->end < rb->end) ? -1 : (ra->end > rb->end) ? 1 : 0;
}
```

`qsort` requires a comparator that returns negative, zero, or positive.
The three keys are applied in priority order:

1. **`seqname`** — `strcmp` already returns the right sign. Records on
   different chromosomes are separated first.
2. **`start`** — within a chromosome, earlier-starting records come first.
   The ternary `(a < b) ? -1 : 1` is used instead of the tempting
   subtraction `a - b` because subtraction can overflow: if `a` is a large
   positive `long` and `b` is a large negative `long`, `a - b` wraps around
   to a negative result, inverting the sort order silently.
3. **`end`** — tiebreaker for records that start at the same position.
   Shorter records sort before longer ones.

**Lexicographic chromosome order.** `strcmp` sorts strings
lexicographically, which means chromosome names sort as `chr1, chr10,
chr11, chr2` rather than `chr1, chr2, chr10, chr11`. This is a known
limitation of naive sorting. Tools like `bedtools` use a reference
genome's chromosome order instead. For our purposes, what matters is
consistency: both the sort order and the hash index use the same string
identity, so a query for `"chr1"` will always find the correct block.

### `gtf_table_sort`

```c
void gtf_table_sort(GtfTable *t)
{
    qsort(t->recs, t->n, sizeof(GtfRecord), rec_cmp);

    htable_free(t->index);
    t->index = htable_new(64);
    for (size_t i = 0; i < t->n; i++) {
        if (i == 0 || t->recs[i].seqname != t->recs[i - 1].seqname)
            htable_insert(t->index, t->recs[i].seqname, i);
    }
}
```

`qsort` sorts the array **in place** — `GtfRecord` values physically move
within `t->recs`. Any pointer taken before sorting may point to a different
record afterwards; this is why `gtf_table_query` is always called *after*
`gtf_table_sort`.

After sorting, the function immediately builds the seqname index. Because
seqnames are interned (§18), consecutive records with the same chromosome
share the same pointer. The boundary check `seqname != prev_seqname` is a
pointer comparison — O(1) and touching no string data. Each unique seqname
is inserted into the hash table with the index of its first record. See §19
for the full hash table discussion.

`qsort` is not guaranteed to be stable. For our purposes this does not
matter — records at identical `(seqname, start, end)` coordinates are
duplicates from the annotation perspective.

---

## 16. `gtf_table_query` — interval overlap queries

### The overlap condition

Two intervals are defined to overlap when they share at least one point.
For 1-based, fully-closed intervals `[a, b]` and `[c, d]`:

```
overlaps  iff  a <= d  AND  c <= b
no overlap iff  a > d  OR   c > b
```

In our notation: record `r` overlaps query `[qstart, qend]` when
`r.start <= qend AND r.end >= qstart`. This is the test applied inside
`gtf_table_query`.

### Finding the start of the seqname block

There are two paths, chosen by whether the hash index has been built:

**Primary path — hash index (O(1)).** After `gtf_table_sort`, `t->index`
maps each seqname to the first record on that chromosome. `htable_get`
finds the slot by hashing and comparing with `strcmp`, then returns a
pointer to the stored index:

```c
if (t->index) {
    size_t *idx = htable_get(t->index, seqname);
    if (!idx) return NULL;          /* chromosome not in table at all */
    lo       = *idx;
    interned = t->recs[lo].seqname; /* canonical pointer for this seqname */
}
```

If `htable_get` returns `NULL`, the chromosome is simply absent — no
scan is needed, return immediately with zero hits.

**Fallback — lower-bound binary search.** Used when the table has been
loaded but not yet sorted (i.e. `t->index == NULL`):

```c
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
```

This is the **lower bound** pattern: the first index where
`seqname >= target`. Unlike `bsearch`, it never returns `NULL` — it returns
the position where the target *would* be inserted. The invariant: every
index `< lo` has `seqname < target`; every index `>= hi` has
`seqname >= target`. When `lo == hi`, that position is the answer.

**`lo + (hi - lo) / 2` vs `(lo + hi) / 2`.** Both compute the midpoint,
but `(lo + hi)` can overflow when both are near `SIZE_MAX`. Subtracting
first keeps the intermediate value in bounds — a classic binary search
correctness detail.

### `gtf_table_query`

```c
GtfRecord **gtf_table_query(const GtfTable *t, const char *seqname,
                             long qstart, long qend, size_t *count)
{
    *count = 0;
    const char *interned = NULL;
    size_t lo;

    if (t->index) {
        size_t *idx = htable_get(t->index, seqname);
        if (!idx) return NULL;
        lo = *idx;
        interned = t->recs[lo].seqname;   /* interned pointer for this seqname */
    } else {
        lo = lower_bound_seqname(t, seqname);
    }

    GtfRecord **hits = malloc(16 * sizeof(GtfRecord *));

    for (size_t i = lo; i < t->n; i++) {
        if (interned) {
            if (t->recs[i].seqname != interned) break; /* (A) pointer compare */
        } else {
            int cmp = strcmp(t->recs[i].seqname, seqname);
            if (cmp > 0) break;                        /* (A) strcmp fallback */
        }
        if (t->recs[i].start > qend) break;            /* (B) past query end */
        if (t->recs[i].end < qstart) continue;         /* (C) before query start */

        hits[(*count)++] = &t->recs[i];
        // grow hits if needed...
    }

    if (*count == 0) { free(hits); return NULL; }
    return hits;
}
```

The loop starts at `lo` and uses three conditions to process each record:

**(A) Seqname boundary — break.** When the index is available, this is a
**pointer comparison**: `t->recs[i].seqname != interned`. Because seqnames
are interned, every record on the target chromosome holds exactly the same
pointer; the first record with a different pointer is on a different
chromosome. When the fallback binary search is used, this falls back to
`strcmp > 0`. The pointer version is O(1) and touches no string data.

**(B) `start > qend` — break.** Within the target chromosome, records are
sorted by `start`. Once a record starts after the query ends, no subsequent
record can overlap. Stop immediately.

**(C) `end < qstart` — continue.** This record starts within the query's
coordinate range (`start <= qend`, since condition B did not fire), but it
ends before the query begins. This seems contradictory, but consider a
long query `[100, 9000]` and a tiny record `[150, 160]`: the record
starts within range but ends long before the query starts... wait, that
cannot happen because `qstart=100 <= end=160`.

The actual scenario for (C): a query `[500, 600]` and a record `[200,
300]`. The record starts before `qstart` but also ends before `qstart`.
After `lower_bound_seqname`, `i` starts at the first record with the right
chromosome — there may be many records that started before the query window
and whose ends also fall before it. Condition (C) skips them individually
while condition (B) will eventually cut off the scan once `start > qend`.

**Result array: pointers, not copies.** `hits` is an array of `GtfRecord *`,
each pointing directly into `t->recs`. No records are copied. This is cheap
and avoids the need for a separate free loop over the records. The consequence
is that the pointers become invalid the moment the table is sorted again or
freed — the caller must use the results before either operation.

**Returning `NULL` for zero results.** When `*count == 0`, the function
frees the hits array and returns `NULL`. An empty array (non-NULL pointer,
`*count == 0`) would also work, but `NULL` is idiomatic in C for "nothing
found" and saves the caller from allocating and freeing an empty buffer.
The caller checks `count > 0` before iterating, so the `NULL` case is
handled naturally.

### `cmd_overlap` — region string parsing

```c
static int cmd_overlap(FILE *fp, const char *region)
{
    const char *colon = strrchr(region, ':');
    // ...
    char seqname[GTF_FIELD_MAX];
    memcpy(seqname, region, seqlen);
    seqname[seqlen] = '\0';

    long qstart, qend;
    if (sscanf(colon + 1, "%ld-%ld", &qstart, &qend) != 2 || qstart > qend) { ... }

    GtfTable *t = gtf_table_load(fp);
    gtf_table_sort(t);

    size_t count;
    GtfRecord **hits = gtf_table_query(t, seqname, qstart, qend, &count);

    for (size_t i = 0; i < count; i++)
        gtf_print(hits[i]);

    free(hits);
    gtf_table_free(t);
}
```

**`strrchr` vs `strchr`.** `strchr(s, c)` finds the *first* occurrence of
`c` in `s`; `strrchr` finds the *last*. A seqname like `chr1:alt_loci`
contains a colon, so `strchr` would split it at the wrong place. `strrchr`
finds the colon that separates seqname from coordinates — the last one —
and handles any number of colons in the seqname itself.

**`sscanf(colon + 1, "%ld-%ld", &qstart, &qend)`.** `colon` points to the
`:` character. `colon + 1` is a pointer to the character immediately after
it — the start of the coordinate string `"1000-5000"`. `sscanf` reads two
`long` integers separated by a literal `-`. The return value is the number
of items successfully parsed; if it is not exactly `2`, the format is wrong.

**`qstart > qend` validation.** An inverted range is a user error. The
check is placed alongside the `sscanf` check so both validation conditions
produce the same error-and-return path.

**`free(hits)` but not the records.** After printing, only the hits array
is freed — not the records it points to. Those belong to the table and will
be freed by `gtf_table_free(t)`. Freeing a pointer you do not own corrupts
the heap; freeing a pointer you do own exactly once is the rule.

---

## 17. Interval merge — a scan on sorted data

### The problem

Genomic annotation files often contain overlapping or adjacent features that
should be collapsed into non-redundant intervals. Three exon records on the
same chromosome with overlapping coordinates represent the same transcribed
region; the merge command reduces them to one interval spanning the combined
range and records how many source intervals were collapsed.

### The algorithm: scan on sorted data

The algorithm is a *running-window sweep* on a sorted sequence:

1. Sort all records by (seqname, start, end) — groups everything by
   chromosome and ensures left-to-right encounter order.
2. Walk through the records maintaining a "current interval" (`cur`).
3. If the next record is on the same chromosome AND starts at or before
   `cur.end + 1` (overlapping or adjacent), extend `cur.end` if needed.
4. Otherwise, emit `cur` and start a new interval from the current record.
5. After the loop, emit the final `cur`.

```c
for (size_t i = 0; i < t->n; i++) {
    GtfRecord *r = &t->recs[i];
    if (strcmp(r->feature, feature) != 0) continue;

    if (!in_merge) { cur = *r; cur_count = 1; in_merge = 1; continue; }

    if (strcmp(r->seqname, cur.seqname) == 0 && r->start <= cur.end + 1) {
        if (r->end > cur.end) cur.end = r->end;   /* extend */
        cur_count++;
    } else {
        /* emit and start over */
        snprintf(merge_attrs, sizeof(merge_attrs),
                 "merged_count \"%zu\";", cur_count);
        cur.attributes = merge_attrs;
        gtf_print(&cur);
        cur = *r; cur_count = 1;
    }
}
if (in_merge) { /* emit the last interval */ }
```

**Why sort first?** Without sorting you would need to compare every record
against every other — O(n²). Sorting costs O(n log n), then the sweep is
O(n). One sort unlocks a single pass, which is the fundamental pattern
behind most "efficient scan" algorithms in bioinformatics.

**The adjacency condition `r->start <= cur.end + 1`.** GTF uses 1-based
closed coordinates. The intervals [900, 1200] and [1201, 1500] share no
bases, but 1201 is immediately after 1200 — they are adjacent with no gap.
Adding 1 merges them. Whether you want this is domain-specific; removing
the `+ 1` gives overlap-only merging.

### `cur = *r` — shallow struct copy

`GtfRecord` contains pointers (`const char *seqname`, `char *attributes`).
`cur = *r` is a *shallow copy*: it copies the pointer values, not the data
they point to. `cur.seqname` now holds the same address as `r->seqname`.

This is intentional and safe here because:
- `seqname`, `source`, `feature` are pool-owned and remain valid until
  `gtf_table_free` — we only read them through `cur`.
- `attributes` is owned by the original table record. We overwrite
  `cur.attributes` with `merge_attrs` before printing, so we never
  read or modify the original through `cur`.

### The local buffer trick

The merged record's attributes need a freshly formatted string for each
emitted interval. Formatting into a `char merge_attrs[64]` on the stack
and temporarily redirecting `cur.attributes` at it avoids any `malloc`:

```c
char merge_attrs[64];   /* lives for the whole function */

/* inside the loop, before gtf_print: */
snprintf(merge_attrs, sizeof(merge_attrs),
         "merged_count \"%zu\";", cur_count);
cur.attributes = merge_attrs;   /* point cur at the stack buffer */
gtf_print(&cur);                /* reads cur.attributes — buffer is alive */
cur = *r;                       /* overwrites the pointer; buffer still fine */
```

`gtf_print` only reads through the pointer. The buffer is alive for the
entire duration of the function, so there is no dangling-pointer risk.
This pattern — temporarily redirecting a pointer to a stack buffer for a
read-only operation — is common in C when you want to avoid heap allocation
for a short-lived value.

### Single-pass streaming output

`cmd_merge` does not collect all merged intervals before printing — it
emits each one as soon as it is finalised. Peak memory is just the
`GtfTable` (needed for sorting) plus a single `GtfRecord cur` on the
stack. The output streams out of the program as it is produced.

---

## 18. String interning — one allocation per unique string

### The memory problem

The original `GtfRecord` had four fixed-size character arrays:

```c
char seqname[256];      /*  256 bytes */
char source[256];       /*  256 bytes */
char feature[256];      /*  256 bytes */
char attributes[8192];  /* 8192 bytes */
                        /* ≈ 9 KB per record */
```

Across 1.7 million records: **~15 GB** just for the record array — well
beyond the RAM of most workstations.

The observation: `seqname`, `source`, and `feature` repeat heavily. A real
annotation file might have fewer than 100 unique chromosome names, fewer
than 10 unique sources, and fewer than 20 unique feature types. Storing a
256-byte copy of `"chr1"` in every one of the 300 000 records that live on
chromosome 1 is pure waste. The `attributes` field, however, is unique per
record — it contains gene IDs, transcript IDs, etc.

**String interning** stores exactly one copy of each unique string in a
pool, and every user of that string holds a pointer to the pooled copy
rather than its own copy. Two records on `chr1` share one `"chr1"` heap
allocation.

### `GtfRawRecord` vs `GtfRecord` — why two structs

The challenge: `gtf_parse_line` writes into a `GtfRecord` using `strncpy`.
If `seqname` becomes a `const char *` pointer rather than a `char[256]`
array, there is nowhere to write into — the pointer has no backing storage.

The solution is to split the type:

```c
/* Scratch buffer — allocated on the stack for one line's parse. */
typedef struct {
    char seqname[GTF_FIELD_MAX];    /* writable arrays */
    char source[GTF_FIELD_MAX];
    char feature[GTF_FIELD_MAX];
    long start, end; float score; char strand; int frame;
    char attributes[GTF_ATTR_MAX];
} GtfRawRecord;

/* Stored in the table — thin pointers, heap-allocated attributes. */
typedef struct {
    const char *seqname;   /* interned — owned by the pool */
    const char *source;    /* interned — owned by the pool */
    const char *feature;   /* interned — owned by the pool */
    long start, end; float score; char strand; int frame;
    char *attributes;      /* heap-allocated — owned by this record */
} GtfRecord;
```

`gtf_parse_line` fills a `GtfRawRecord` on the stack — a temporary that
disappears when the function returns. `gtf_table_load` then converts it to
a `GtfRecord` by interning the three short fields and copying attributes:

```c
GtfRecord *r = &t->recs[t->n];
r->seqname    = str_pool_intern(&t->pool, raw.seqname);
r->source     = str_pool_intern(&t->pool, raw.source);
r->feature    = str_pool_intern(&t->pool, raw.feature);
r->attributes = strndup_safe(raw.attributes, strlen(raw.attributes));
```

`stats.c` and `compare.c` use `GtfRawRecord` directly (they parse
line-by-line without building a table) and needed a single change each:
`GtfRecord rec` → `GtfRawRecord rec`.

### `StrPool` — the intern pool

```c
typedef struct {
    char  **data;   /* heap array of heap-allocated strings */
    size_t  n;
    size_t  cap;
} StrPool;
```

The pool is a dynamic array of `char *`. `str_pool_intern` does a linear
scan; if the string already exists it returns the stored pointer, otherwise
it appends a new copy:

```c
static const char *str_pool_intern(StrPool *pool, const char *s)
{
    for (size_t i = 0; i < pool->n; i++)
        if (strcmp(pool->data[i], s) == 0) return pool->data[i];

    /* not found — grow the pool and append */
    /* ... realloc pool->data ... */
    char *dup = strndup_safe(s, strlen(s));
    pool->data[pool->n++] = dup;
    return dup;
}
```

Linear scan is O(unique strings) per intern call. With ~100 unique seqnames
and ~20 unique features, this is at most ~120 comparisons per record —
negligible against the parsing overhead.

### `const char *` vs `char *` — communicating ownership

The interned fields are `const char *`: the caller may read but not write
or free them. The pool is the sole owner. `attributes` is `char *` because
it is uniquely owned by each record.

This is not just style. If `seqname` were `char *`, nothing in the type
system would stop a careless `free(rec->seqname)` — which would corrupt
the pool and cause a double-free when the pool itself is cleaned up.
`const char *` makes the prohibition explicit.

### Updated `gtf_table_free`

The free cascade now has three layers:

```c
void gtf_table_free(GtfTable *t) {
    /* 1. Each record's uniquely-owned attributes string */
    for (size_t i = 0; i < t->n; i++)
        free(t->recs[i].attributes);
    free(t->recs);

    /* 2. The pool's interned strings */
    for (size_t i = 0; i < t->pool.n; i++)
        free(t->pool.data[i]);
    free(t->pool.data);

    /* 3. The hash index (see §19) */
    htable_free(t->index);
    free(t);
}
```

Order matters: always free the contents before the container that holds
them. Freeing `t->recs` first (before freeing attributes) would not cause
a bug here because the attribute pointers are copied into each `GtfRecord`
individually — but it would make the free loop impossible.

### Memory impact

After interning, each `GtfRecord` is ~40 bytes of fixed-size fields plus
a pointer. The `attributes` string averages ~100 bytes (the real text, not
a fixed 8 KB buffer). Per record: ~140 bytes. Across 1.7 M records:
**~240 MB** — a **60× reduction** from the original fixed-buffer layout.

### Pointer equality as string equality

Once strings are interned, the same string value has exactly one address.
This means:

```c
/* Before interning: O(string_length) */
strcmp(rec_a->seqname, rec_b->seqname) == 0

/* After interning: O(1) */
rec_a->seqname == rec_b->seqname
```

This is more than a micro-optimisation. It is the key that makes the hash
table index (§19) fast: detecting where one seqname block ends and the
next begins during index construction becomes a pointer comparison — no
string data touched at all.

---

## 19. Hash table — O(1) seqname lookup

### Why a hash table here

`gtf_table_query` needs to jump to the first record with a given seqname
in the sorted array. The previous implementation used a binary search:
O(log n). A hash table gives O(1). For 1.7 M records, log₂(1 700 000) ≈ 21
steps vs one hash lookup — the absolute speed difference is small, but the
design teaches every fundamental hash table concept: hash functions,
collision resolution, load factor, resize, and pointer-identity keys.

The index also resolves an ambiguity in the binary search: the lower-bound
search could land in the middle of a seqname block or just before a
seqname that does not exist. The hash table returns `NULL` immediately for
absent chromosomes, and the exact first index for chromosomes that do exist.

### Open addressing with linear probing

Two main collision strategies:
- **Chaining**: each slot holds a linked list of entries.
- **Open addressing**: all entries live in the same flat array; on
  collision, probe adjacent slots.

We use open addressing with *linear probing* — on collision, try `h+1`,
`h+2`, etc., wrapping around. This keeps all data in one contiguous
allocation (cache-friendly) and requires no `malloc` per entry.

```c
typedef struct { const char *key; size_t val; } HEntry;

typedef struct {
    HEntry *slots;
    size_t  cap;   /* always a power of 2 */
    size_t  n;     /* live entries */
} HTable;
```

An empty slot is `key == NULL`. `calloc` in `htable_new` zeros the slot
array, so all slots start empty without an explicit initialisation loop.

**Why power-of-2 capacity?** The slot index is computed as
`hash & (cap - 1)` instead of `hash % cap`. Bitwise AND is a single
instruction; modulo by an arbitrary number requires division. The trick
only works when `cap` is a power of 2 because `cap - 1` is then a bitmask
with all bits below the capacity bit set.

### FNV-1a — the hash function

```c
static size_t fnv1a(const char *s)
{
    size_t h = 14695981039346656037ULL;   /* FNV offset basis */
    for (; *s; s++) {
        h ^= (unsigned char)*s;
        h *= 1099511628211ULL;            /* FNV prime */
    }
    return h;
}
```

FNV-1a processes one byte at a time: XOR the byte into the hash, then
multiply by a chosen prime. The XOR-before-multiply order (FNV-1a) has
better *avalanche* — a one-bit change in the input changes roughly half
the output bits — than multiply-before-XOR (FNV-1). The constants are the
standard 64-bit FNV parameters, derived to avoid common patterns that
produce clustering. No division, no branching, no tables.

### Load factor and resize

As a hash table fills up, the average *probe length* (how many slots must
be checked before finding the target or an empty slot) grows. The
conventional threshold is **75%**: when `n * 4 > cap * 3`, double the
capacity and rehash.

```c
static int htable_insert(HTable *ht, const char *key, size_t val)
{
    if (ht->n * 4 > ht->cap * 3)
        htable_resize(ht, ht->cap * 2);

    size_t h = fnv1a(key) & (ht->cap - 1);
    while (ht->slots[h].key) {
        /* pointer equality: all stored keys are interned */
        if (ht->slots[h].key == key) { ht->slots[h].val = val; return 0; }
        h = (h + 1) & (ht->cap - 1);
    }
    ht->slots[h].key = key;
    ht->slots[h].val = val;
    ht->n++;
    return 0;
}
```

`htable_resize` allocates a new slot array and re-inserts every live entry
using the new capacity mask:

```c
static int htable_resize(HTable *ht, size_t new_cap)
{
    HEntry *new_slots = calloc(new_cap, sizeof(HEntry));
    for (size_t i = 0; i < ht->cap; i++) {
        if (!ht->slots[i].key) continue;
        size_t h = fnv1a(ht->slots[i].key) & (new_cap - 1);
        while (new_slots[h].key) h = (h + 1) & (new_cap - 1);
        new_slots[h] = ht->slots[i];
    }
    free(ht->slots);
    ht->slots = new_slots;
    ht->cap   = new_cap;
    return 0;
}
```

Each entry is re-hashed because the slot index `hash & (cap-1)` changes
when `cap` doubles. Old slot positions are invalid in the new array.

### Building the index in `gtf_table_sort`

The index maps each seqname to the index of its first record in the sorted
array. After sorting, records with the same seqname are contiguous — the
first occurrence of each new seqname marks a block boundary.

```c
void gtf_table_sort(GtfTable *t)
{
    qsort(t->recs, t->n, sizeof(GtfRecord), rec_cmp);

    htable_free(t->index);
    t->index = htable_new(64);

    for (size_t i = 0; i < t->n; i++) {
        /* pointer inequality: interned strings share one address */
        if (i == 0 || t->recs[i].seqname != t->recs[i - 1].seqname)
            htable_insert(t->index, t->recs[i].seqname, i);
    }
}
```

Because `seqname` is interned (§18), detecting block boundaries is a
pointer comparison — O(1), no string data accessed. This is the payoff of
combining interning and hashing: the two techniques reinforce each other.

### Using the index in `gtf_table_query`

```c
if (t->index) {
    size_t *idx = htable_get(t->index, seqname);
    if (!idx) return NULL;          /* chromosome not in this file */
    lo       = *idx;
    interned = t->recs[lo].seqname; /* the canonical pool pointer */
} else {
    lo = lower_bound_seqname(t, seqname);   /* fallback */
}
```

`htable_get` uses `strcmp` because the caller's `seqname` (e.g., a stack
buffer in `cmd_overlap`) is not guaranteed to be the interned pointer:

```c
static size_t *htable_get(const HTable *ht, const char *key)
{
    size_t h = fnv1a(key) & (ht->cap - 1);
    while (ht->slots[h].key) {
        if (strcmp(ht->slots[h].key, key) == 0) return &ht->slots[h].val;
        h = (h + 1) & (ht->cap - 1);
    }
    return NULL;
}
```

But once the slot is found, `t->recs[lo].seqname` IS the interned pointer.
The inner scan loop then uses pointer equality instead of `strcmp`:

```c
for (size_t i = lo; i < t->n; i++) {
    if (t->recs[i].seqname != interned) break;  /* O(1): pointer compare */
    if (t->recs[i].start > qend)        break;
    if (t->recs[i].end   < qstart)      continue;
    /* collect hit */
}
```

`strcmp` appears exactly once (in `htable_get`) and only for as many probe
steps as there are hash collisions — near zero for a well-loaded table.
Every other string comparison in the hot path is gone.

### Ownership: the table borrows the keys

`htable_free` frees the slot array and the `HTable` struct, but **not**
the strings the slots point to:

```c
void htable_free(HTable *ht)
{
    if (!ht) return;
    free(ht->slots);
    free(ht);
}
```

The keys are interned pointers owned by `t->pool`. The hash table *borrows*
them — it did not allocate them and must not free them. The pool is freed
separately in `gtf_table_free`. This is the same owner/borrower pattern
that governs `gtf_table_query`'s result array: the caller owns the array
of pointers but not the `GtfRecord` objects the pointers point to.

### Three layers working together

| Layer | What it does | Key concept |
|---|---|---|
| `StrPool` (§18) | One allocation per unique string | Pointer equality = string equality |
| `gtf_table_sort` | Sorts + builds seqname→index map | Pointer comparison to find block edges |
| `gtf_table_query` | Hash lookup + pointer scan | O(1) jump, O(1) boundary check |

The three ideas form a chain: interning enables cheap equality; cheap
equality enables the hash index; the hash index collapses O(log n) binary
search into O(1) lookup and turns the inner scan loop from `strcmp`-per-
record into a pointer comparison.

---

## §20 BED output — coordinate conversion and format interop

### Why BED?

GTF and GFF3 are the standard exchange formats for gene annotation, but most
downstream tools in a genome analysis pipeline — bedtools, IGV, UCSC Genome
Browser, deeptools — consume **BED** instead.  BED is simpler: no attributes
column, no format variants, and it uses the coordinate system (0-based,
half-open) that C arrays and most programming languages use naturally.

Implementing a GTF→BED converter is therefore not busywork: it is one of the
most-used single operations in annotation-adjacent pipelines.

### BED coordinate system

The fundamental difference between GTF and BED:

| Property | GTF/GFF3 | BED |
|---|---|---|
| Start is | 1-based | 0-based |
| End is | inclusive (closed) | exclusive (open, "half-open") |
| A single base at position 1 | `start=1, end=1` | `start=0, end=1` |
| Length | `end − start + 1` | `end − start` |

Both representations carry the same information; the conversion is:

```
BED start = GTF start − 1
BED end   = GTF end          (no change)
```

Example: GTF `chr1 . gene 1001 2000` → BED `chr1 1000 2000 . 0 .`

The length is `2000 − 1001 + 1 = 1000` in GTF and `2000 − 1000 = 1000` in
BED — identical, as expected.

The 0-based convention exists because C arrays are 0-indexed.  A 100-byte
array has valid indices `[0, 100)` — that is, `arr[0]` through `arr[99]`.
BED mirrors that: to address bases 1000–2000 of a chromosome stored as a
C string, you would write `seq + 1000` and copy `2000 - 1000 = 1000` bytes.
With GTF you would write `seq + start - 1` — the off-by-one is always there
when interfacing with C strings.

### BED6 columns

The `bed` command produces BED6:

| # | Name | Source |
|---|------|--------|
| 1 | chrom | `GtfRecord.seqname` |
| 2 | chromStart | `GtfRecord.start − 1` |
| 3 | chromEnd | `GtfRecord.end` |
| 4 | name | attribute value (via `--name`) or `.` |
| 5 | score | integer 0–1000 (GTF score clamped), or `0` when NaN |
| 6 | strand | `GtfRecord.strand` (`+`, `-`, or `.`) |

### Score field: NaN, clamping, and casting

GTF stores the score as a floating-point number, but the BED spec requires
an integer in `[0, 1000]`.  GTF also uses `.` to mean "no score", which
the parser stores as `NAN` (from `<math.h>`).  Three things need to happen:

```c
int score = 0;
if (!isnan(r->score)) {
    float s = r->score;
    score = (s < 0.0f) ? 0 : (s > 1000.0f) ? 1000 : (int)s;
}
```

- `isnan()` detects the "no score" sentinel; those records get `0`.
- The ternary clamps out-of-range values; the cast truncates toward zero.
- `(int)s` is a *narrowing conversion* — this is safe here because we have
  already established `0 ≤ s ≤ 1000`, so the result fits in an `int`.

Without the clamp, a GTF score of `9999.7` would truncate to `9999`, violating
the BED spec.  Without `isnan`, casting `NAN` to `int` is undefined behaviour
in C — the result is implementation-defined and may trap on some architectures.

### Name field: attribute extraction per record

When `--name` is given, the implementation calls `gtf_attrs_parse` (or
`gff3_attrs_parse` for GFF3) for every record to extract one attribute value:

```c
GtfAttrs *attrs = (t->format == GTF_FMT_GFF3)
                  ? gff3_attrs_parse(r->attributes)
                  : gtf_attrs_parse(r->attributes);
if (attrs) {
    const char *val = gtf_attrs_get(attrs, name_attr);
    if (val) name = val;
}
```

This is the same pattern used by `cmd_attrs`.  Parsing costs one `malloc`
per record, which is acceptable for an export command that is typically run
once.  If performance becomes a concern, a targeted single-attribute scanner
(like the private `attr_get` in `stats.c`) could replace it without changing
the calling code.

The `name` pointer is valid until `gtf_attrs_free(attrs)` is called (the
string is part of the `GtfAttrs` allocation).  `printf` runs before
`gtf_attrs_free`, so the pointer is live.

### Dispatch: the --name flag pattern

`bed` uses the same flag-scanning dispatch pattern as `print` (§5):

```c
for (int i = 2; i < argc; i++) {
    if (strcmp(argv[i], "--name") == 0 && i + 1 < argc) {
        name_attr = argv[++i];               // consume next token
    } else if (strncmp(argv[i], "--name=", 7) == 0) {
        name_attr = argv[i] + 7;             // pointer arithmetic into argv[i]
    } else if (argv[i][0] != '-') {
        path = argv[i];                      // first non-flag token is the file
    } else {
        fprintf(stderr, "error: unknown option '%s'\n", argv[i]);
        return 1;
    }
}
```

Both `--name gene_id` (two tokens) and `--name=gene_id` (one token) are
supported.  `argv[i] + 7` is pointer arithmetic: skip the first 7 characters
(`--name=`) and point directly to the value inside the existing string.
No allocation, no copy.

### Composing with other commands

Because `bed` reads from a `FILE *` and writes to stdout, it can be piped
with other `gtfparse` commands.  The `filter` command output is a valid GTF
stream, which `bed` can read from `/dev/stdin`:

```sh
# Only gene features, named by gene_id
./gtfparse filter gene file.gtf | ./gtfparse bed --name gene_id /dev/stdin

# Merge overlapping exons, then convert to BED
./gtfparse merge exon file.gtf | ./gtfparse bed /dev/stdin
```

This composability is the payoff of the design decision to always read from
`FILE *` rather than a file path: any command can feed any other.

---

## §21 Output redirection — argv compaction and freopen

### The problem

Every command writes its output to `stdout` using `printf` or `gtf_print`.
The user can always redirect that with the shell (`> file.gtf`), but that
requires the shell.  Inside a script or when the flag position must be
controlled programmatically, a built-in `-o`/`--output` flag is more
convenient and self-documenting.

The challenge is that the flag must work with *every* command without
modifying any `cmd_*` function.

### stdout is just a FILE *

In C, `stdin`, `stdout`, and `stderr` are global variables of type `FILE *`,
defined in `<stdio.h>`.  They are initialised at program start by the C
runtime, but they are just ordinary pointers — you can replace what they
point to.

`freopen` does exactly that:

```c
FILE *freopen(const char *path, const char *mode, FILE *stream);
```

It closes `stream`'s current file descriptor, opens `path` in `mode`, and
re-uses the `FILE` object at the same address.  After the call, any code
that writes to `stdout` — including code in libraries you didn't write —
automatically writes to `path` instead.  The pointer value of `stdout`
doesn't change; only the underlying descriptor does.

```c
if (!freopen(out_path, "w", stdout)) {
    perror(out_path);   // prints: "out.gtf: No such file or directory"
    return 1;
}
```

`freopen` returns `NULL` on failure, leaving `stdout` in an indeterminate
state.  We treat failure as fatal and exit immediately.

### argv compaction

The flag must be removed from `argv` before any dispatch code runs, because
the `print` and `bed` dispatch blocks reject unknown flags:

```c
} else {
    fprintf(stderr, "error: unknown option '%s'\n", argv[i]);
    return 1;
}
```

The solution is to compact `argv` in-place: when the flag is found, shift
all later elements left by the number of consumed tokens and decrement
`argc`.

```c
const char *out_path = NULL;
for (int i = 1; i < argc; ) {
    if ((strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0)
            && i + 1 < argc) {
        out_path = argv[i + 1];
        for (int j = i; j < argc - 2; j++)   /* consume 2 tokens */
            argv[j] = argv[j + 2];
        argc -= 2;
    } else if (strncmp(argv[i], "--output=", 9) == 0) {
        out_path = argv[i] + 9;
        for (int j = i; j < argc - 1; j++)   /* consume 1 token */
            argv[j] = argv[j + 1];
        argc--;
    } else {
        i++;
    }
}
```

A few details worth noting:

**The loop does not increment `i` when it removes elements.**  After
shifting, the element that moved into position `i` is new and hasn't been
examined yet.  Incrementing would skip it.  This is the standard pattern
for in-place filtering of an array: only advance when you *keep* the current
element.

**`argv` is a `char *[]` — an array of pointers.**  Shifting does not copy
any string data; it only moves the pointers.  The strings themselves live
in read-only memory provided by the OS and are not touched.

**`argc` is a local variable.**  The C standard says `main`'s `argc` is a
local copy; modifying it does not affect anything outside `main`.  All
dispatch code after the pre-pass uses the modified `argc`, which no longer
counts the consumed flag tokens.

**The `--output=value` form uses pointer arithmetic.**  `argv[i] + 9` skips
the first nine characters (`--output=`) and points directly into the middle
of the existing string.  No allocation, no `strncpy`.

### Why a pre-pass rather than per-command handling

An alternative design would add `--output` parsing to each command's own
flag loop.  That would work but requires touching every dispatch block now
and every new one in the future.  The pre-pass keeps the concern in one
place: the `cmd_*` functions never need to know that an output file exists.

This is an instance of a general principle: **cross-cutting concerns belong
at a single boundary**.  Authentication, logging, and output redirection are
all examples — they affect every operation but should be handled once, not
scattered through every operation's implementation.

### stderr stays on the terminal

Diagnostic output — loading progress, record counts, error messages — all
use `fprintf(stderr, ...)`.  `freopen` only replaces `stdout`; `stderr` is
unaffected.  This means:

```sh
./gtfparse filter gene annotation.gtf -o genes.gtf
# terminal shows:  loading: 100000 records...
#                  matched=14234  total=1873083
# genes.gtf gets:  the GTF records
```

The Unix convention of separating data (stdout) from diagnostics (stderr)
is what makes this work cleanly.  It is also why piping works: in
`cmd1 | cmd2`, the shell connects cmd1's stdout to cmd2's stdin while both
commands' stderrs remain on the terminal.

---

## §22 sort — exposing library code as a command

### The function

```c
static int cmd_sort(FILE *fp)
{
    GtfTable *t = gtf_table_load(fp);
    if (!t) { fprintf(stderr, "error: out of memory\n"); return 1; }

    fprintf(stderr, "  sorting %zu records...\n", t->n);
    gtf_table_sort(t);

    for (size_t i = 0; i < t->n; i++)
        gtf_print(&t->recs[i]);

    fprintf(stderr, "records=%zu  format=%s\n", t->n,
            t->format == GTF_FMT_GFF3 ? "gff3" : "gtf");
    gtf_table_free(t);
    return 0;
}
```

This is the shortest `cmd_*` function in the codebase.  It calls
`gtf_table_sort` (§15), which already exists and is already used internally
by `overlap` and `merge`.  The only new work is the loop that calls
`gtf_print` for every record after sorting.

### Library code vs command code

The codebase has two layers:

- **Library layer** (`gtf.c`): functions that manipulate data — `gtf_table_load`,
  `gtf_table_sort`, `gtf_table_query`, `gtf_attrs_parse`, etc.  These know
  nothing about command-line arguments, flags, or what the user typed.
- **Command layer** (`main.c`): functions that parse arguments, open files,
  call library functions, and format output for the user.

`gtf_table_sort` existed before `sort` was a user command because `overlap`
needed it internally.  Adding the `sort` command required zero changes to the
library layer.  This is the payoff of keeping the two layers separate: library
functions accumulate value over time and can be reused in ways not anticipated
when they were written.

The same principle applies in larger systems.  A database engine's sort
routine can be exposed as `ORDER BY` in SQL, as a sort utility in the CLI,
and called internally by join algorithms — all without modifying the sort
routine itself.

### Why sort is useful as a command

**Reproducibility.**  GTF files from different sources, or produced by
different runs of the same tool, may have records in different orders.
Sorting normalises the order so that `diff` and checksum-based comparisons
work correctly.

**Interoperability.**  Tools like `bedtools` require sorted input.
`./gtfparse sort file.gtf | bedtools ...` feeds pre-sorted data to bedtools
without a separate sort step.

**Inspection.**  It is easier to read a large annotation file when all records
for a chromosome appear together and in coordinate order.

**Prerequisite visibility.**  `overlap` sorts internally before querying.
The explicit `sort` command lets users perform that step once, save the result,
and run multiple queries against the sorted file without re-sorting each time.

### The sort key

`gtf_table_sort` calls `qsort` with `rec_cmp`, which compares:

1. `seqname` — lexicographic (`strcmp`)
2. `start` — numeric (ascending)
3. `end` — numeric (ascending)

`seqname` comparison is lexicographic, which means `chr10` sorts before
`chr2` (because `'1' < '2'`).  This is the same behaviour as standard Unix
`sort` without `-V`.  Natural sort (`chr1, chr2, ..., chr10`) would require
a custom comparator that splits the string into text and numeric segments —
a more complex function, and one we have not needed yet.

---

## §23 keys — dynamic string deduplication and goto

### What the command does

`keys` answers the question: *what attribute keys exist anywhere in this
file?*  It parses every record's attributes column, collects unique key
names into a dynamic array, sorts them, and prints one per line.

### The key array

The unique keys are stored in a local `char **keys` dynamic array — the same
doubling-growth pattern used by `GtfTable` and `StrPool`, but local to
`cmd_keys` rather than in a named struct:

```c
char  **keys = NULL;
size_t  nk = 0, capk = 0;
```

`keys` starts as `NULL`.  The first time a new key is found, `realloc` grows
it to 16 slots; subsequent growth doubles the capacity.  Each entry is a
heap-allocated copy of the key string.

**Why copy?**  The key string lives inside a `GtfAttrs` struct which is freed
at the end of each record's iteration.  A pointer into a freed allocation is a
dangling pointer.  Copying with `malloc + strcpy` gives the key array its own
stable storage that outlives the per-record loop.

### Deduplication by linear scan

Before inserting a new key, the code scans every existing entry:

```c
int found = 0;
for (size_t m = 0; m < nk; m++)
    if (strcmp(keys[m], k) == 0) { found = 1; break; }
if (found) continue;
```

This is O(u) per candidate key, where u is the number of unique keys already
seen.  For the common case — GTF files have roughly 5–25 distinct attribute
keys — this is negligible compared to the O(n) cost of loading 1.7 million
records.  Replacing the linear scan with a hash table would be an
over-engineering for this problem size.

The choice between linear scan and hash table is a question of expected input
size.  When u is bounded and small, a linear scan is simpler, has lower
constant factors, and uses less memory.  The hash table in `gtf_table_sort`
is justified because n (number of chromosomes) can reach thousands; u here
is bounded by the format specification itself.

### goto for error cleanup

When a `malloc` or `realloc` fails mid-loop, the code needs to free everything
already allocated and return an error.  Without `goto`, this requires either
deeply nested `if`s or duplicating the cleanup code on every error path.

The `goto done` pattern consolidates all cleanup in one place:

```c
    if (!tmp) { gtf_attrs_free(a); rc = 1; goto done; }
    ...
    if (!copy) { gtf_attrs_free(a); rc = 1; goto done; }

done:
    if (rc) fprintf(stderr, "error: out of memory\n");

    for (size_t i = 0; i < nk; i++) free(keys[i]);
    free(keys);
    gtf_table_free(t);
    return rc;
```

The label `done:` is placed just before the cleanup block.  Every error path
sets `rc = 1` and jumps there; the normal path falls through to `done:`
naturally after the loop ends.  The cleanup runs unconditionally in both cases.

`goto` has a bad reputation from its misuse in early languages, but in C it
is the idiomatic way to handle multi-step cleanup on error.  The Linux kernel
uses this pattern in nearly every function that allocates more than one
resource.  The rule is: `goto` only jumps *forward*, never backward, and the
label is always near the end of the function.

### Sorting char ** with qsort

After collecting all unique keys, they are sorted with:

```c
qsort(keys, nk, sizeof(char *), cmp_str_ptr);
```

The comparator receives two `const void *` pointers, each pointing to one
element of the `char **` array — so each is a `char **`:

```c
static int cmp_str_ptr(const void *a, const void *b)
{
    return strcmp(*(const char **)a, *(const char **)b);
}
```

The double dereference `*(const char **)a` casts `a` from `void *` to
`char **`, then dereferences to get the `char *` string pointer.  This is the
standard pattern for sorting any array of pointers with `qsort` — the
comparator always receives a pointer *to the element*, not the element itself.

---

## §24 validate — per-record checks, feature sets, and exit codes

### What the command checks

Each record is tested against four independent groups of rules:

| Group | Condition checked |
|-------|------------------|
| `coord` | `start >= 1`; `end >= start` |
| `strand` | value is `+`, `-`, or `.` |
| `frame` | CDS / start_codon / stop_codon must have frame 0, 1, or 2 |
| `attr` | presence of required attribute keys (format-dependent) |

Each violation is printed to stdout immediately when found, with the record
number and check group name.  After all records are scanned, a summary line
goes to stderr.

### NULL-terminated feature sets

Several checks need to test whether the current record's feature belongs to
a predefined set.  Rather than a chain of `||` comparisons, the code uses
`NULL`-terminated `const char *` arrays and a small helper:

```c
static const char *const cds_like[] = { "CDS", "start_codon", "stop_codon", NULL };

static int feat_in(const char *feat, const char *const *set)
{
    for (; *set; set++)
        if (strcmp(feat, *set) == 0) return 1;
    return 0;
}
```

`static const char *const cds_like[]` — there are two layers of `const`:
- `const char *` — the string each pointer points to is read-only
- `const *` at the array level — the pointers themselves are read-only

`static` on the array means it lives in static storage (once, for the
lifetime of the program) rather than being rebuilt on the stack each time
the function is called.  Because the array never changes, this is the right
storage class.

The sentinel `NULL` at the end lets `feat_in` walk the array without needing
to know its length — a classic C idiom used by `argv`, `execv`, and many
library APIs.

### Why check frame only for CDS-like features

The GTF/GFF3 frame field specifies how many bases at the start of the feature
need to be skipped to reach the first complete codon.  It is only meaningful
for features that encode protein: `CDS`, `start_codon`, `stop_codon`.  For
`gene`, `exon`, and other non-coding features the field is conventionally `.`
(stored as `-1`), so checking frame == -1 on those would produce false
positives.

### Exit code as a signal

The function returns `nerrors > 0 ? 1 : 0`:

```c
return nerrors > 0 ? 1 : 0;
```

Returning a non-zero exit code when errors are found is the Unix convention
for tools that make a pass/fail determination (`diff`, `test`, `grep`).  It
lets the validator be used directly in shell scripts without parsing stdout:

```sh
if ./gtfparse validate annotation.gtf > /dev/null 2>&1; then
    echo "clean — proceeding with pipeline"
fi
```

This is distinct from returning an error code to signal a *program* failure
(like `malloc` returning `NULL`).  The distinction matters: exit code 1 here
means "the file has problems", not "the program crashed".

### What this validator does not check

Per-record checks catch individual malformed records but cannot detect
cross-record inconsistencies.  A complete validator would also check:

- Every `transcript_id` referenced by an exon must have a corresponding
  transcript record.
- Exons within a transcript must not overlap each other.
- CDS intervals must be contained within exon intervals.
- `Parent=` values in GFF3 must resolve to a record with a matching `ID=`.

These require grouping records by ID and comparing intervals across the
group — a hash-table-of-record-lists problem.  The current implementation
deliberately stops before that complexity; those checks are the natural next
step and a good exercise in combining the `HTable` (§19) with dynamic arrays
(§11).
