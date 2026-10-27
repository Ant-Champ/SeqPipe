#!/usr/bin/env python3
"""sv_scan -- large indels and structural breakpoints from a sorted, indexed BAM.

Ported from alignqc.py (step_indels). Two evidence types (short reads cannot
span big indels, so both are needed):
  1. CIGAR indels >= min-indel bp, clustered by position   -> --out-indels
  2. Soft-clip pile-ups (>= min-clip bp) = breakpoints from large deletions /
     insertions / rearrangements; n_split counts reads that also carry an SA
     (supplementary alignment) tag                          -> --out-breakpoints

Coordinates are 0-based half-open (BED-like). Empty BAMs / no events give
header-only tables.

--out-indels columns ('#' header line):
  #contig  start  end  type(INS|DEL)  n_reads  median_len  min_len  max_len
  depth  frac_reads
--out-breakpoints columns ('#' header line):
  #contig  start  end  side(L|R)  n_reads  n_split  depth  frac_reads

Python >= 3.8; needs pysam.
"""
import argparse
import statistics
import sys
from collections import defaultdict
from pathlib import Path


class Progress:
    """Minimal progress bar on stderr (no extra dependency). Redraws in place
    on a terminal; in logs (not a tty) prints one line per 10%."""
    def __init__(self, label, total):
        self.label, self.total, self.last = label, max(total, 1), -1
        self.tty = sys.stderr.isatty()

    def set(self, n):
        pct = min(100, 100 * n // self.total)
        shown = pct if self.tty else pct // 10 * 10
        if shown != self.last:
            self.last = shown
            print(f"{self.label} [{'#' * (pct // 4):<25}] {pct:3d}%  ({n:,})",
                  end="\r" if self.tty else "\n", file=sys.stderr, flush=True)

    def close(self, n):
        self.set(self.total if n >= self.total else n)
        if self.tty:
            print(file=sys.stderr)


def _cluster(events, tol):
    """Group position-sorted tuples (pos, ...) whose positions are within tol."""
    if not events:                                  # guard: nothing to group
        return []
    events.sort()
    clusters, cur = [], [events[0]]
    for e in events[1:]:
        if e[0] - cur[-1][0] <= tol:
            cur.append(e)
        else:
            clusters.append(cur)
            cur = [e]
    clusters.append(cur)
    return clusters


def run(args):
    import pysam
    for p in (args.out_indels, args.out_breakpoints):
        Path(p).parent.mkdir(parents=True, exist_ok=True)

    indel_events = defaultdict(list)            # (contig, 'INS'|'DEL') -> [(pos, len)]
    clip_events = defaultdict(lambda: [0, 0])   # (contig, pos, side) -> [n_reads, n_split]

    bam = pysam.AlignmentFile(args.bam, "rb")
    try:
        total = bam.mapped                      # from the .bai, instant
    except (ValueError, AttributeError):
        total = 0
    prog = Progress("scanning reads", total)
    n = 0
    for n, r in enumerate(bam.fetch(), 1):
        if n % 200_000 == 0:
            prog.set(n)
        # unmapped|secondary|qcfail|duplicate|supplementary in one test
        if r.flag & 0xF04 or r.mapping_quality < args.min_mapq:
            continue
        ct = r.cigartuples
        if not ct or len(ct) == 1:              # plain e.g. 75M: no indel, no clip
            continue
        contig, ref = r.reference_name, r.reference_start
        for op, ln in ct:
            if op in (0, 7, 8, 3):              # M, =, X, N consume reference
                ref += ln
            elif op == 1 and ln >= args.min_indel:              # insertion
                indel_events[(contig, "INS")].append((ref, ln))
            elif op == 2:                                       # deletion
                if ln >= args.min_indel:
                    indel_events[(contig, "DEL")].append((ref, ln))
                ref += ln
        split = 1 if r.has_tag("SA") else 0
        if ct[0][0] == 4 and ct[0][1] >= args.min_clip:         # left soft clip
            e = clip_events[(contig, r.reference_start, "L")]
            e[0] += 1; e[1] += split
        if ct[-1][0] == 4 and ct[-1][1] >= args.min_clip:       # right soft clip
            e = clip_events[(contig, r.reference_end - 1, "R")]
            e[0] += 1; e[1] += split
    prog.close(n)
    print("writing indel / breakpoint tables")

    def depth_at(contig, pos):
        return bam.count(contig, max(pos, 0), max(pos, 0) + 1, read_callback="all")

    # ---- indels
    with open(args.out_indels, "w") as out:
        out.write("#contig\tstart\tend\ttype\tn_reads\tmedian_len\tmin_len\tmax_len\t"
                  "depth\tfrac_reads\n")
        for (contig, kind), ev in sorted(indel_events.items()):
            for c in _cluster(ev, args.cluster_tol):
                if len(c) < args.min_support:
                    continue
                start = int(statistics.median(x[0] for x in c))
                lens = [x[1] for x in c]
                med = int(statistics.median(lens))
                end = start + (med if kind == "DEL" else 1)
                d = depth_at(contig, start)
                out.write(f"{contig}\t{start}\t{end}\t{kind}\t{len(c)}\t{med}\t"
                          f"{min(lens)}\t{max(lens)}\t{d}\t{len(c)/max(d,1):.2f}\n")

    # ---- breakpoints (merge jittered neighbouring positions on the same side)
    by_side = defaultdict(list)
    for (contig, pos, side), (nr, ns) in clip_events.items():
        by_side[(contig, side)].append((pos, nr, ns))
    contig_len = dict(zip(bam.references, bam.lengths))
    with open(args.out_breakpoints, "w") as out:
        out.write("#contig\tstart\tend\tside\tn_reads\tn_split\tdepth\tfrac_reads\n")
        for (contig, side), ev in sorted(by_side.items()):
            for c in _cluster(ev, args.cluster_tol):
                nr = sum(x[1] for x in c)
                if nr < args.min_support:
                    continue
                pos = max(c, key=lambda x: x[1])[0]             # strongest position
                # clips at the contig ends = reads crossing the origin of a circular replicon
                if pos < args.edge or pos >= contig_len[contig] - args.edge:
                    continue
                d = depth_at(contig, pos)
                out.write(f"{contig}\t{pos}\t{pos+1}\t{side}\t{nr}\t{sum(x[2] for x in c)}\t"
                          f"{d}\t{nr/max(d,1):.2f}\n")
    bam.close()
    print(f"done -> {args.out_indels}, {args.out_breakpoints}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="sorted, indexed BAM")
    ap.add_argument("--out-indels", required=True, help="indels TSV")
    ap.add_argument("--out-breakpoints", required=True, help="breakpoints TSV")
    ap.add_argument("--min-mapq", type=int, default=20, help="minimum mapping quality [20]")
    ap.add_argument("--min-indel", type=int, default=10, help="minimum CIGAR indel size (bp) [10]")
    ap.add_argument("--min-clip", type=int, default=20,
                    help="minimum soft-clip length to call a breakpoint (bp) [20]")
    ap.add_argument("--min-support", type=int, default=3, help="minimum reads supporting an event [3]")
    ap.add_argument("--cluster-tol", type=int, default=5,
                    help="position tolerance when clustering events (bp) [5]")
    ap.add_argument("--edge", type=int, default=10,
                    help="ignore breakpoints this close to a contig end (origin of circular replicons); 0 = keep all [10]")
    args = ap.parse_args()
    try:
        run(args)
    except ImportError as e:
        sys.exit(f"sv_scan: missing dependency: {e}")
    except (OSError, ValueError) as e:
        sys.exit(f"sv_scan: error: {e}")


if __name__ == "__main__":
    main()
