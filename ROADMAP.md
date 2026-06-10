# GTFParty — Planned Features

Features identified as commonly needed in real annotation pipelines, ordered
roughly from simplest to most complex.  Each entry notes the primary C concept
it would introduce.

---

## Commands to add

| # | Command | Synopsis | C concepts introduced |
|---|---------|----------|-----------------------|
| 1 | ~~`sort`~~ ✓ | Sort records by (seqname, start, end) and write back out | Expose existing `gtf_table_sort`; teach that library code can be surfaced as a command with no new logic |
| 2 | `count` | Count records per feature type | Hash table of `feature → size_t`; formatted table output |
| 3 | ~~`keys`~~ ✓ | List all unique attribute keys in a file | Nested string scan; second hash table use-case (key deduplication) |
| 4 | ~~`validate`~~ ✓ | Check coordinate sanity, required attributes, strand values | State machine over records; structured error reporting to stderr |
| 5 | `intersect` | Records in file A that overlap any record in file B | Two-pointer merge on sorted arrays; classic sweep-line algorithm |
| 6 | `subtract` | Records in file A with no overlap in file B | Same sweep-line, inverted predicate |
| 7 | `introns` | Derive intron intervals from exon records per transcript | Per-transcript grouping; interval gap arithmetic |
| 8 | `utrs` | Derive UTR intervals from exon and CDS records | Multi-feature grouping; set difference of intervals |
| 9 | `extract` | Extract FASTA sequences for each feature given a genome | FASTA parsing; file seeking (`fseek`/`ftell`); reverse-complement |
| 10 | `hierarchy` | Print gene → transcript → exon trees (GFF3 Parent= chain) | Tree / linked-list data structure; recursive or stack-based traversal |

---

## Notes

**`sort`** is nearly free — `gtf_table_sort` already exists.  The only work
is adding a dispatch block and a `gtf_print`-per-record loop.

**`count` and `keys`** are good next steps: they reuse the `HTable` from §19
for a different purpose (counting instead of indexing), reinforcing the idea
that a data structure is a tool usable in many contexts.

**`intersect` / `subtract`** are the most-requested missing features in
annotation tooling (after BED output).  The two-pointer merge algorithm is a
fundamental technique: sort both inputs, advance two pointers in tandem, and
emit or suppress records based on overlap.  No extra memory proportional to
file size is needed.

**`introns` and `utrs`** teach interval arithmetic: introns are the gaps
between consecutive exons within a transcript (`prev.end + 1` to
`next.start - 1`); UTRs are the parts of exons outside the CDS interval.
Both require grouping records by transcript, which is a hash-table-of-lists
problem.

**`extract`** introduces a new I/O pattern: random access.  GTF coordinates
tell you exactly where in the genome sequence to look, but a FASTA file is
text with line-wrapping.  A common approach is a first pass to build a
byte-offset index (chromosome name → file offset), then `fseek` for each
feature.

**`hierarchy`** is the most structurally complex: GFF3's `Parent=` attribute
forms a directed acyclic graph.  Building it requires a hash map from ID to
node, a second pass to wire up parent pointers, and a traversal that respects
the tree order.
