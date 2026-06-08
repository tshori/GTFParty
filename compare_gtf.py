#!/usr/bin/env python3
"""
compare_gtf.py — compare two GTF files at gene and transcript level.

Comparison is by ID (gene_id / transcript_id strings), not coordinates.
This is the right approach when comparing two versions of the same annotation
(e.g. Ensembl v110 vs v111). Coordinate-based comparison is a different problem.

Dependencies (only for --venn):
    pip install matplotlib matplotlib-venn
"""
import argparse
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _attr(attributes, key):
    """Extract a GTF attribute value by key (quoted or unquoted)."""
    m = re.search(rf'\b{re.escape(key)} "([^"]+)"', attributes)
    if m:
        return m.group(1)
    m = re.search(rf'\b{re.escape(key)} (\S+?)(?:;|$)', attributes)
    if m:
        return m.group(1)
    return None


def parse_gtf(path):
    """
    Parse a GTF file.

    Returns
    -------
    records        : list of (raw_line, gene_id, transcript_id)
    gene_ids       : set of all gene_id strings seen
    transcript_ids : set of all transcript_id strings seen
    """
    records       = []
    gene_ids      = set()
    transcript_ids = set()

    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                print(f"  warning: {path}:{lineno}: fewer than 9 fields — skipped",
                      file=sys.stderr)
                continue
            gid = _attr(parts[8], 'gene_id')
            tid = _attr(parts[8], 'transcript_id')
            if gid:
                gene_ids.add(gid)
            if tid:
                transcript_ids.add(tid)
            records.append((line, gid, tid))

    return records, gene_ids, transcript_ids


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def _pct(n, total):
    return f"{100 * n / total:.1f}%" if total else "N/A"


def print_summary(label1, label2, genes1, genes2, txs1, txs2):
    shared_g = genes1 & genes2
    only_g1  = genes1 - genes2
    only_g2  = genes2 - genes1

    shared_t = txs1 & txs2
    only_t1  = txs1 - txs2
    only_t2  = txs2 - txs1

    # Summary goes to stderr so stdout stays clean for GTF output when piping.
    out = sys.stderr

    print(f"  GTF1 = {label1}", file=out)
    print(f"  GTF2 = {label2}", file=out)
    print(file=out)

    W = 22
    print("=== Gene-level comparison ===", file=out)
    print(f"  {'GTF1':<{W}}: {len(genes1):>10,}", file=out)
    print(f"  {'GTF2':<{W}}: {len(genes2):>10,}", file=out)
    print(f"  {'Shared':<{W}}: {len(shared_g):>10,}  "
          f"({_pct(len(shared_g), len(genes1))} of GTF1, "
          f"{_pct(len(shared_g), len(genes2))} of GTF2)", file=out)
    print(f"  {'Only in GTF1':<{W}}: {len(only_g1):>10,}  "
          f"({_pct(len(only_g1), len(genes1))} of GTF1)", file=out)
    print(f"  {'Only in GTF2':<{W}}: {len(only_g2):>10,}  "
          f"({_pct(len(only_g2), len(genes2))} of GTF2)", file=out)
    print(file=out)
    print("=== Transcript-level comparison ===", file=out)
    print(f"  {'GTF1':<{W}}: {len(txs1):>10,}", file=out)
    print(f"  {'GTF2':<{W}}: {len(txs2):>10,}", file=out)
    print(f"  {'Shared':<{W}}: {len(shared_t):>10,}  "
          f"({_pct(len(shared_t), len(txs1))} of GTF1, "
          f"{_pct(len(shared_t), len(txs2))} of GTF2)", file=out)
    print(f"  {'Only in GTF1':<{W}}: {len(only_t1):>10,}  "
          f"({_pct(len(only_t1), len(txs1))} of GTF1)", file=out)
    print(f"  {'Only in GTF2':<{W}}: {len(only_t2):>10,}  "
          f"({_pct(len(only_t2), len(txs2))} of GTF2)", file=out)


# ---------------------------------------------------------------------------
# GTF output
# ---------------------------------------------------------------------------

def write_output(recs1, recs2, genes1, genes2, mode, out_path):
    """
    Write a filtered or merged GTF.

    merge     — all GTF1 records, plus GTF2 records whose gene_id is absent
                from GTF1. Shared genes keep their GTF1 annotation.
    only-in-1 — GTF1 records whose gene_id does not appear in GTF2.
    only-in-2 — GTF2 records whose gene_id does not appear in GTF1.
    """
    out = open(out_path, 'w') if out_path else sys.stdout
    n = 0
    try:
        if mode == 'merge':
            for line, gid, _ in recs1:
                print(line, file=out)
                n += 1
            for line, gid, _ in recs2:
                if gid not in genes1:
                    print(line, file=out)
                    n += 1
        elif mode == 'only-in-1':
            for line, gid, _ in recs1:
                if gid not in genes2:
                    print(line, file=out)
                    n += 1
        elif mode == 'only-in-2':
            for line, gid, _ in recs2:
                if gid not in genes1:
                    print(line, file=out)
                    n += 1
    finally:
        if out_path:
            out.close()

    dest = out_path if out_path else "stdout"
    print(f"Wrote {n:,} records [{mode}] → {dest}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Venn diagram
# ---------------------------------------------------------------------------

def make_venn(genes1, genes2, txs1, txs2, out_path, label1, label2):
    try:
        import matplotlib.pyplot as plt
        from matplotlib_venn import venn2
    except ImportError:
        print(
            "Warning: matplotlib-venn not installed — Venn diagram skipped.\n"
            "Install with:  pip install matplotlib matplotlib-venn",
            file=sys.stderr,
        )
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"{label1}  vs  {label2}", fontsize=12)

    venn2([genes1, genes2], set_labels=(label1, label2), ax=axes[0])
    axes[0].set_title("Genes", fontsize=11)

    venn2([txs1, txs2], set_labels=(label1, label2), ax=axes[1])
    axes[1].set_title("Transcripts", fontsize=11)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Venn diagram saved → {out_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("gtf1", help="First GTF file")
    ap.add_argument("gtf2", help="Second GTF file")
    ap.add_argument(
        "--output",
        choices=["merge", "only-in-1", "only-in-2"],
        metavar="{merge,only-in-1,only-in-2}",
        help=(
            "Write a GTF to --out-file (or stdout). "
            "'merge' = all of GTF1 plus genes unique to GTF2; "
            "'only-in-1' = genes absent from GTF2; "
            "'only-in-2' = genes absent from GTF1."
        ),
    )
    ap.add_argument("--out-file", metavar="FILE",
                    help="Output GTF path (default: stdout)")
    ap.add_argument(
        "--venn", metavar="FILE",
        help="Save a two-panel Venn diagram (genes + transcripts) to FILE. "
             "Extension sets format: .png, .pdf, .svg",
    )
    ap.add_argument("--label1", metavar="NAME",
                    help="Display name for GTF1 (default: filename)")
    ap.add_argument("--label2", metavar="NAME",
                    help="Display name for GTF2 (default: filename)")
    args = ap.parse_args()

    label1 = args.label1 or Path(args.gtf1).name
    label2 = args.label2 or Path(args.gtf2).name

    print(f"Parsing {args.gtf1} ...", file=sys.stderr)
    recs1, genes1, txs1 = parse_gtf(args.gtf1)
    print(f"Parsing {args.gtf2} ...", file=sys.stderr)
    recs2, genes2, txs2 = parse_gtf(args.gtf2)
    print(file=sys.stderr)

    print_summary(label1, label2, genes1, genes2, txs1, txs2)

    if args.output:
        print(file=sys.stderr)
        write_output(recs1, recs2, genes1, genes2, args.output, args.out_file)

    if args.venn:
        make_venn(genes1, genes2, txs1, txs2, args.venn, label1, label2)


if __name__ == "__main__":
    main()
