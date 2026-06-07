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
6. [Pointers — a complete picture](#6-pointers--a-complete-picture)
7. [Memory lifetimes and dangling pointers](#7-memory-lifetimes-and-dangling-pointers)
8. [Known limitations and what Step 2 fixes](#8-known-limitations-and-what-step-2-fixes)

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

### Why multiple files?

Splitting code across files is not just organisation — it controls what each
part of the program can see and depend on.

- `gtf.h` is the **public interface**: the struct definition and the two
  function signatures. Any file that includes it knows what a `GtfRecord` is
  and what `gtf_parse_line` expects.
- `gtf.c` is the **implementation**: the actual code. It is the only file
  that needs to know how parsing works internally.
- `main.c` is the **consumer**: it uses the interface without caring about
  implementation details.

This separation means you can change how parsing works inside `gtf.c` without
touching `main.c`, as long as the function signatures in `gtf.h` stay the same.

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
gtfparse: main.o gtf.o
    $(CC) $(CFLAGS) -o gtfparse main.o gtf.o $(LDFLAGS)

main.o: main.c gtf.h
    $(CC) $(CFLAGS) -c main.c

gtf.o: gtf.c gtf.h
    $(CC) $(CFLAGS) -c gtf.c
```

**Stage 1 — compile** (`-c` flag): each `.c` file is translated into machine
code and written to a `.o` (object) file. At this stage, calls to functions
in other files are left as unresolved placeholders.

**Stage 2 — link**: the linker combines the `.o` files, resolves every
placeholder, and writes the final executable.

The benefit of this split: `make` tracks which files have changed. If you
edit only `gtf.c`, only `gtf.o` is recompiled. `main.o` is reused as-is.
On large projects with hundreds of files this matters enormously.

You can inspect what symbols a `.o` file exports and imports:

```sh
nm gtf.o
```

Output:
```
T gtf_parse_line    # T = defined (exported) by this file
T gtf_print
U printf            # U = undefined (imported from elsewhere)
U strtok
U atof
```

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
will be silently truncated. Step 2 replaces these arrays with `char *` heap
pointers to remove that limit.

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
exists. See section 7 for a detailed discussion of this.

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
separation means that `./gtfparse file.gtf > out.gtf` sends records to the
file while counts still appear in the terminal.

---

## 5. `main.c` — the entry point

### Command-line arguments

```c
int main(int argc, char *argv[])
```

`argc` (argument count) is the number of strings on the command line,
including the program name itself. `argv` (argument vector) is an array of
pointers to those strings:

```
./gtfparse file.gtf

argv[0] → "./gtfparse\0"
argv[1] → "file.gtf\0"
argv[2] → NULL           (argv is always NULL-terminated)
argc    == 2
```

### File I/O

```c
FILE *fp = fopen(argv[1], "r");
if (!fp) {
    perror(argv[1]);
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
./gtfparse bad.gtf > /dev/null || echo "parse failed"
```

---

## 6. Pointers — a complete picture

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

---

## 7. Memory lifetimes and dangling pointers

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

### Heap memory (preview of Step 2)

`malloc(n)` allocates `n` bytes on the **heap** — memory that persists until
explicitly freed with `free()`, regardless of which function allocated it.
This removes the lifetime problem but introduces a new one: you must call
`free` exactly once for every `malloc`, or you get a memory leak. Step 2
will work through this in detail.

---

## 8. Known limitations and what Step 2 fixes

### Buffer size caps

`seqname`, `source`, and `feature` are capped at 255 characters.
`attributes` is capped at 8191 characters. Any value longer than these
limits is silently truncated. In practice, GTF files from NCBI RefSeq fit
within these limits, but the design is fragile.

**Fix in Step 2:** replace `char field[N]` with `char *field` and allocate
exactly as much heap memory as each string needs (`strlen(src) + 1` bytes).

### Line length cap

`buf[65536]` limits lines to 65535 characters. GTF attributes columns can
be long on heavily annotated genomes, but 64 KB is generous in practice.

**Fix in Step 2:** use `getline()`, which allocates a buffer large enough
for each line dynamically.

### `strtok` is not thread-safe

`strtok` uses a hidden static pointer. Calling it from two threads
simultaneously corrupts both calls. `strtok_r` is the thread-safe
alternative.

### No attribute parsing

The `attributes` field is stored as a raw string. Accessing individual
attributes (e.g. `gene_id`) requires a second pass of string parsing.

**Fix in Step 4:** parse `key "value";` pairs into an array of
`{ char *key; char *value; }` structs.
