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

### The struct

```c
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
} GtfRecord;
```

A `struct` groups related values into a single named unit. Without `typedef`
you would write `struct GtfRecord rec` everywhere. The `typedef` lets you
write just `GtfRecord rec`.

**Fixed-size char arrays** (`char seqname[256]`) are the simplest way to
store strings in C. The entire storage is embedded directly inside the struct
— no pointers, no heap allocation, no `malloc`. When you declare `GtfRecord rec`
as a local variable, the compiler reserves exactly `sizeof(GtfRecord)` bytes
on the stack.

The trade-off: any string longer than 255 characters (plus null terminator)
will be silently truncated. In practice, NCBI RefSeq GTF files fit within
these limits. A future step could replace these arrays with heap-allocated
`char *` pointers to remove the cap entirely.

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
int  gtf_parse_line(const char *line, GtfRecord *rec);
void gtf_print(const GtfRecord *rec);
```

These are **declarations** (also called prototypes), not definitions. They
tell the compiler what types the functions accept and return, so it can
type-check call sites in `main.c` before it has seen `gtf.c`. The actual
code lives in `gtf.c`.

`const char *line` — a pointer to char, with `const` on the data. The
function promises not to modify the string it receives.

`GtfRecord *rec` — a pointer to a `GtfRecord`. The function will write
results into the struct the caller provides. This is how C functions
"return" multiple values — pass a pointer to where you want the output.

### GtfTable and GtfAttrs

`gtf.h` also declares two higher-level types built on top of `GtfRecord`:

**`GtfTable`** holds a heap-allocated, dynamically-grown array of every
record in a file. It is the in-memory representation of the whole file,
as opposed to streaming one record at a time. See section 11.

**`GtfAttr` / `GtfAttrs`** represent the parsed attributes column: a
dynamic array of `{char *key; char *value;}` pairs allocated on the heap
on demand. Rather than scanning the raw attribute string each time you need
a value, you parse once and look up by key. See section 13.

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

#### Copying fields into the struct

```c
strncpy(rec->seqname, fields[0], GTF_FIELD_MAX - 1);
rec->seqname[GTF_FIELD_MAX - 1] = '\0';
```

`rec` is a pointer to the caller's `GtfRecord`. The `->` operator is
shorthand: `rec->seqname` means `(*rec).seqname` — dereference the pointer
to get the struct, then access the field. Writing to `rec->seqname` writes
directly into the caller's memory; no copy of the struct is made.

We copy the bytes that `fields[0]` points to into the struct's fixed array.
This is the crucial step that makes our data safe — once copied, the struct
fields no longer depend on `buf`, which disappears when the function returns.

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
    if (strcmp(cmd, "compare") == 0) {
        if (argc < 4) { ... }
        return cmd_compare(argv[2], argv[3]);
    }

    /* filter / attrs: one string argument then a file */
    if (strcmp(cmd, "filter") == 0 || strcmp(cmd, "attrs") == 0) {
        if (argc < 4) { ... }
        FILE *fp = fopen(argv[3], "r");
        int rc = (strcmp(cmd, "filter") == 0)
                 ? cmd_filter(fp, argv[2])
                 : cmd_attrs(fp, argv[2]);
        fclose(fp);
        return rc;
    }

    /* print / stats: one file argument */
    if (argc < 3) { usage(argv[0]); return 1; }
    FILE *fp = fopen(argv[2], "r");
    // ...
    if      (strcmp(cmd, "print") == 0) ret = cmd_print(fp);
    else if (strcmp(cmd, "stats") == 0) ret = cmd_stats(fp);
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

**Three argument shapes.** Commands fall into three groups based on how
many arguments they need, and the shapes are checked from most-specific to
least-specific:

| Shape | Commands | `argv` layout |
|-------|----------|---------------|
| Two files | `compare` | `argv[2]`, `argv[3]` are both files |
| String + file | `filter`, `attrs` | `argv[2]` is a string, `argv[3]` is a file |
| One file | `print`, `stats` | `argv[2]` is a file |

Checking the special cases first is a common C idiom: each block either
returns or falls through to the next, so the final `fopen(argv[2])` is only
reached by commands that actually take exactly one file argument. If `compare`
were not checked first, `fopen(argv[2])` would try to open `"a.gtf"` as the
file for a two-file command — wrong.

**Sharing a dispatch case.** `filter` and `attrs` have identical argument
shapes (`<string> <file>`), so they are handled by a single `if` block.
A ternary chooses which function to call:

```c
int rc = (strcmp(cmd, "filter") == 0)
         ? cmd_filter(fp, argv[2])
         : cmd_attrs(fp, argv[2]);
```

This is cleaner than duplicating the `fopen`/`fclose` pair. The pattern
scales: when a third command with the same shape is added, only the ternary
needs to grow.

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

### Buffer size caps in `GtfRecord` — still present

`seqname`, `source`, and `feature` are capped at 255 characters.
`attributes` is capped at 8191 characters. Any value longer than these
limits is silently truncated. In practice, GTF files from NCBI RefSeq
(tested on 1.7 million records) fit within these limits.

Removing the cap requires replacing `char field[N]` with `char *field`
and heap-allocating each string with `malloc(strlen(src) + 1)`. That change
would turn `GtfRecord` from a value type (copyable by assignment) into a
type that owns heap memory and must be explicitly freed — a larger API
change that is deferred to a future step.

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

### Sorting and interval queries — addressed in Steps 5 and 7

`gtf_table_sort` sorts the table in place by `(seqname, start, end)`.
`gtf_table_query` uses binary search to locate a chromosome then scans
forward, giving O(log n + k) per query. See sections 15 and 16.

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
} GtfTable;
```

- `recs` — pointer to the heap-allocated array of records.
- `n` — the number of valid records currently stored.
- `cap` — the allocated capacity (always `>= n`).
- `n_skipped` — count of comment and blank lines encountered while loading.
- `n_errors` — count of lines that failed to parse.

The last two fields carry diagnostic information without requiring the
caller to maintain separate counters. They are written to `stderr` by
commands that care (`cmd_print`), and ignored by commands that do not
(`cmd_filter`, `cmd_attrs`).

### `gtf_table_load`

```c
GtfTable *gtf_table_load(FILE *fp)
{
    GtfTable *t = malloc(sizeof(GtfTable));
    // ...
    t->cap  = GTF_TABLE_INIT_CAP;   /* 4096 */
    t->recs = malloc(t->cap * sizeof(GtfRecord));

    char line[65536];
    GtfRecord rec;

    while (fgets(line, sizeof(line), fp)) {
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
```

**Three-way return from `gtf_parse_line`.** The function returns `1` for
comments/blanks (increment `n_skipped`, continue), `-1` for parse errors
(increment `n_errors`, continue), and `0` for a valid record. Only `0`
causes a record to be added to the table.

**Grow before insert.** The capacity check comes *before* storing the new
record. The invariant is: `t->n < t->cap` at the point of assignment.
Checking after would risk writing one past the end of the array.

**`realloc` via a temporary.** As explained in section 6, assigning
`realloc`'s return value directly to `t->recs` would leak the old allocation
if `realloc` returned `NULL`. The temporary `tmp` avoids this.

**`t->recs[t->n++] = rec` — struct assignment.** C copies a struct by value
when you assign it with `=`. Every field of `rec` is copied into the next
slot of the array. After this line, the table's copy is independent of the
local `rec` — modifying `rec` on the next iteration does not affect what
was stored.

**Failure path calls `gtf_table_free`.** If `realloc` fails mid-load, we
call `gtf_table_free(t)` before returning `NULL`. This frees the partially
filled table rather than leaking it. `gtf_table_free` handles a partially
initialised `GtfTable` safely because it only frees `t->recs` and `t` — both
of which are valid at that point.

### `gtf_table_free`

```c
void gtf_table_free(GtfTable *t)
{
    if (!t) return;
    free(t->recs);
    free(t);
}
```

This is simpler than `cmp_free_transcripts` (section 7) because `GtfRecord`
is a value type — it contains no pointers to separately allocated memory, so
there is nothing to free inside each record. A single `free(t->recs)` releases
the entire array at once.

Compare with `Transcript`, which has an `IntronSpan *introns` pointer in
each element: that design requires a loop to free each intron array before
freeing the outer array. `GtfRecord`'s fixed-size embedded arrays avoid
this complexity, at the cost of the size caps discussed in section 10.

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

**`strcmp` on a struct field.** `t->recs[i].feature` is a `char` array
embedded in the `GtfRecord`. `strcmp` takes two `const char *` arguments —
passing an array name decays to a pointer to its first element automatically.
This is the array-to-pointer decay rule: in most expression contexts, an
array name becomes a pointer to its first element.

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
consistency: `gtf_table_query` uses the same `strcmp` to binary-search, so
the sort order and search order always agree.

### `gtf_table_sort`

```c
void gtf_table_sort(GtfTable *t)
{
    qsort(t->recs, t->n, sizeof(GtfRecord), rec_cmp);
}
```

A one-line wrapper. `qsort` sorts the array **in place** — the `GtfRecord`
values physically move within `t->recs`. This has an important consequence:
any pointer you took into the array before sorting (`GtfRecord *p =
&t->recs[42]`) may now point to a completely different record. This is why
`gtf_table_query` is always called *after* `gtf_table_sort`, and why
`cmd_overlap` does not store any record pointers before sorting.

`qsort` is not guaranteed to be stable (equal elements may appear in any
order relative to each other). For our purposes this does not matter —
records at identical `(seqname, start, end)` coordinates are duplicates
from the annotation perspective.

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

### The lower-bound binary search

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

This is the **lower bound** pattern: find the first index where
`seqname >= target`. It is distinct from `bsearch`, which finds *an*
element equal to the key (and returns `NULL` if none exists). Lower bound
always returns a valid index, even when the target is not present — it
returns the position where the target *would* be inserted to keep the
array sorted.

The invariant the loop maintains is: every index `< lo` has
`seqname < target`, and every index `>= hi` has `seqname >= target`.
When `lo == hi`, that single position is the answer.

**`lo + (hi - lo) / 2` vs `(lo + hi) / 2`.** Both compute the midpoint,
but `(lo + hi)` can overflow when `lo` and `hi` are large `size_t` values
near `SIZE_MAX`. Subtracting first and adding `lo` afterward keeps the
intermediate value within bounds. This is a classic binary search bug;
the safe form is always preferred.

### `gtf_table_query`

```c
GtfRecord **gtf_table_query(const GtfTable *t, const char *seqname,
                             long qstart, long qend, size_t *count)
{
    *count = 0;
    size_t lo = lower_bound_seqname(t, seqname);

    size_t cap = 16;
    GtfRecord **hits = malloc(cap * sizeof(GtfRecord *));

    for (size_t i = lo; i < t->n; i++) {
        int cmp = strcmp(t->recs[i].seqname, seqname);
        if (cmp > 0) break;                      /* (A) past target seqname */
        if (t->recs[i].start > qend) break;      /* (B) start past query end */
        if (t->recs[i].end < qstart) continue;   /* (C) ends before query start */

        hits[(*count)++] = &t->recs[i];
        // grow hits if needed...
    }

    if (*count == 0) { free(hits); return NULL; }
    return hits;
}
```

The loop starts at `lo` — the first record that could possibly be on the
target chromosome — and uses three conditions to process the rest:

**(A) `cmp > 0` — break.** The records are sorted by `(seqname, start)`.
Once `seqname` is lexicographically greater than the target, every
subsequent record is also on a later chromosome. There can be no more
matches; stop immediately.

**(B) `start > qend` — break.** Within the target chromosome, records are
sorted by `start`. Once a record starts after the query ends, no subsequent
record can overlap (their starts are even later). Stop immediately.

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
