#!/usr/bin/env python3
"""read_snippets -- a compact, genome-browser-style sample of the aligned reads
for the HTML report (zoom in on a coverage plot to see them).

Where: the whole reference when it is at most --whole-ref-max bp (plasmids,
phages, PCR products); otherwise windows of +/- --flank bp around every event
of the --events tables (variants, indels, breakpoints, low / high coverage
edges, cassette hotspots / sites), merged, at most --max-windows of them.

Which reads: primary alignments (no secondary / supplementary), sorted by
start and packed into at most --max-rows rows (a read that finds no free row
is dropped: like a browser's downsampling, every position keeps up to
--max-rows reads). The rows bound the output; --max-reads (default 0 = off)
is an extra cap per contig.

What is kept of a read: start, end (1-based, inclusive), flags (1 = reverse
strand, 2 = MAPQ < --min-mapq), soft clips on each side, and only what differs
from the reference: mismatches (base), insertions (length), deletions (length).

Output JSON: {contig: {"windows": [[s, e], ...], "rows": n,
  "ref": [[start, "SEQ"], ...]            reference bases of each window
  "reads": [[row, start, end, flags, lclip, rclip, ops...], ...]}}
ops = flat triples (type, pos, value): 0 = mismatch (value 0..4 = A C G T N),
1 = insertion after pos (value = length), 2 = deletion starting at pos (length).

Needs pysam.
"""
import argparse
import csv
import json
import sys

BASE = {"A": 0, "C": 1, "G": 2, "T": 3}


def read_events(paths, lengths):
    """(contig, start, end) 0-based half-open from the event tables: a CHROM/POS
    column pair (variants) or contig/start/end as the first columns."""
    ev = []
    for path in paths:
        try:
            with open(path, newline="") as fh:
                rows = [r for r in csv.reader(fh, delimiter="\t") if r and r[0].strip()]
        except OSError:
            continue
        if not rows:
            continue
        head = [h.lstrip("#") for h in rows[0]]
        body = rows[1:] if (rows[0][0].startswith("#") or "CHROM" in head or "contig" in head) else rows
        for r in body:
            if r[0].startswith("#"):
                continue
            try:
                if "CHROM" in head and "POS" in head:
                    c, p = r[head.index("CHROM")], int(r[head.index("POS")]) - 1
                    ev.append((c, p, p + 1))
                else:
                    c, s, e = r[0], int(r[1]), int(r[2])
                    if e - s > 2000:              # a long region: its two edges
                        ev += [(c, s, s + 1), (c, e - 1, e)]
                    else:
                        ev.append((c, s, max(e, s + 1)))
            except (ValueError, IndexError):
                continue
    return [x for x in ev if x[0] in lengths]


def windows(events, lengths, flank, max_windows):
    by = {}
    for c, s, e in sorted(events):
        w = [max(0, s - flank), min(lengths[c], e + flank)]
        lst = by.setdefault(c, [])
        if lst and w[0] <= lst[-1][1]:
            lst[-1][1] = max(lst[-1][1], w[1])
        else:
            lst.append(w)
    n = sum(len(v) for v in by.values())
    if n > max_windows:
        print(f"WARNING: {n} windows, keeping the first {max_windows} (--max-windows)", file=sys.stderr)
        kept, out = 0, {}
        for c, lst in by.items():
            out[c] = lst[:max(0, max_windows - kept)]
            kept += len(out[c])
        by = {c: v for c, v in out.items() if v}
    return by


def read_ops(r, ref, min_bq):
    """Differences of one alignment from the reference: flat (type, pos1, value)."""
    ops, q, p = [], r.query_sequence or "", r.reference_start
    quals = r.query_qualities
    qi = 0
    for op, n in r.cigartuples:
        if op in (0, 7, 8):                       # M = X
            for k in range(n):
                rb, qb = ref[p + k] if p + k < len(ref) else "N", q[qi + k] if qi + k < len(q) else "N"
                if qb != rb and rb != "N" and (quals is None or quals[qi + k] >= min_bq):
                    ops += [0, p + k + 1, BASE.get(qb, 4)]
            p += n
            qi += n
        elif op == 1:                             # I
            ops += [1, p, n]
            qi += n
        elif op in (2, 3):                        # D N
            ops += [2, p + 1, n]
            p += n
        elif op == 4:                             # S
            qi += n
    return ops


def run(args):
    import pysam
    bam = pysam.AlignmentFile(args.bam)
    fa = pysam.FastaFile(args.fasta)
    lengths = dict(zip(bam.references, bam.lengths))
    bad = [c for c in lengths if c in (".", "*")]
    if bad:
        sys.exit(f"read_snippets: contig(s) named {bad}: rerun from a fresh out_dir (see stage_ref)")
    total = sum(lengths.values())
    if total <= args.whole_ref_max:
        wins = {c: [[0, L]] for c, L in lengths.items()}
        print(f"reference {total:,} bp <= {args.whole_ref_max:,}: whole reference")
    else:
        wins = windows(read_events(args.events, lengths), lengths, args.flank, args.max_windows)
        print(f"reference {total:,} bp: {sum(len(v) for v in wins.values())} window(s) "
              f"(+/- {args.flank} bp around events)")
    out = {}
    for c, ws in wins.items():
        ref = fa.fetch(c).upper()
        # Reads are streamed in start order (windows are sorted and do not
        # overlap) and packed as they come: nothing per contig is held in memory
        # but the rows' ends and the reads kept.
        row_end, packed, n_seen, prev_end, full = [], [], 0, 0, False
        for s, e in ws:
            for r in bam.fetch(c, s, e):
                if r.flag & 0x904 or not r.cigartuples:
                    continue
                if r.reference_start < prev_end:     # already seen in the previous window
                    continue
                n_seen += 1
                ct = r.cigartuples
                lclip = ct[0][1] if ct[0][0] == 4 else 0
                rclip = ct[-1][1] if ct[-1][0] == 4 else 0
                s0 = r.reference_start - lclip               # room for the clipped tails
                for i, e0 in enumerate(row_end):
                    if s0 > e0 + 1:
                        row = i
                        break
                else:
                    if len(row_end) >= args.max_rows:
                        continue
                    row_end.append(0)
                    row = len(row_end) - 1
                row_end[row] = r.reference_end + rclip
                flags = (1 if r.is_reverse else 0) | (2 if r.mapping_quality < args.min_mapq else 0)
                packed.append([row, r.reference_start + 1, r.reference_end, flags, lclip, rclip]
                              + read_ops(r, ref, args.min_bq))
                if args.max_reads and len(packed) >= args.max_reads:
                    print(f"WARNING: {c}: --max-reads {args.max_reads} reached at {r.reference_end:,}: "
                          "no reads stored after it", file=sys.stderr)
                    full = True
                    break
            prev_end = e
            if full:
                break
        out[c] = {"windows": [[s + 1, e] for s, e in ws], "rows": len(row_end),
                  "ref": [[s + 1, ref[s:e]] for s, e in ws], "reads": packed}
        print(f"{c}: {n_seen} read(s) in the windows, {len(packed)} kept in {len(row_end)} row(s)")
    with open(args.out, "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(f"-> {args.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True, help="reference FASTA (indexed)")
    ap.add_argument("--events", nargs="*", default=[], help="TSV tables of events (see above)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--whole-ref-max", type=int, default=200000)
    ap.add_argument("--flank", type=int, default=500)
    ap.add_argument("--max-windows", type=int, default=300)
    ap.add_argument("--max-rows", type=int, default=40)
    ap.add_argument("--max-reads", type=int, default=0,
                    help="stop after this many reads per contig (0 = no limit: --max-rows bounds the output)")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-bq", type=int, default=0, help="mismatches below this base quality are not kept")
    run(ap.parse_args())


if __name__ == "__main__":
    main()
