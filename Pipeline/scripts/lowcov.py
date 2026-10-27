#!/usr/bin/env python3
"""lowcov -- per-contig coverage summary + merged low-coverage regions.

Ported from alignqc.py (step_lowcov). Depth comes from `samtools depth -aa -Q
<min-mapq>` (every position, including zeros and reference contigs with no
reads) and is held in memory at 4 bytes/base: fine for bacteria/plasmids/phage.
A base is "low" if depth < max(low_depth, low_frac * median depth); low bases
closer than --low-merge-gap are merged and regions shorter than --low-min-len
are dropped.

Coordinates are 0-based half-open (BED-like).

--out-summary columns ('#' header line):
  #contig  length  mean_depth  median_depth  pct_zero  pct_low
--out-low columns ('#' header line):
  #contig  start  end  length  mean_depth  min_depth

Python >= 3.8; needs samtools on PATH.
"""
import argparse
import subprocess
import sys
from array import array
from collections import Counter
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


def _median_from_hist(hist, n):
    seen, half = 0, n / 2
    for k in sorted(hist):
        seen += hist[k]
        if seen >= half:
            return k


def low_regions(arr, thr, merge_gap, min_len):
    """Merged low-depth regions of one contig -> [(start, end, mean, min)]."""
    regions = []                                # [start, end, sum, n, min]
    for i, d in enumerate(arr):
        if d >= thr:
            continue
        if regions and i - regions[-1][1] <= merge_gap:
            r = regions[-1]
            r[1], r[2], r[3], r[4] = i + 1, r[2] + d, r[3] + 1, min(r[4], d)
        else:
            regions.append([i, i + 1, d, 1, d])
    return [(s, e, tot / n, mn) for s, e, tot, n, mn in regions if e - s >= min_len]


def run(args):
    with open(args.fai) as fai:
        genome_len = sum(int(l.split("\t")[1]) for l in fai if l.strip())

    # ---- pass 1: read depths, contig -> array (index = pos-1)
    # -aa (not -a): also report contigs with no reads at all, e.g. a lost plasmid
    depth = {}
    prog = Progress("reading depth", genome_len)
    try:
        proc = subprocess.Popen(["samtools", "depth", "-aa", "-Q", str(args.min_mapq), args.bam],
                                stdout=subprocess.PIPE, text=True)
    except OSError as e:
        sys.exit(f"lowcov: cannot run samtools: {e}")
    n = 0
    for n, line in enumerate(proc.stdout, 1):
        if n % 500_000 == 0:
            prog.set(n)
        contig, _pos, d = line.rstrip("\n").split("\t")
        depth.setdefault(contig, array("I")).append(int(d))
    prog.close(n)
    if proc.wait() != 0 or not depth:
        sys.exit("lowcov: samtools depth failed or returned nothing")

    hist = Counter()
    for arr in depth.values():
        hist.update(arr)
    median = _median_from_hist(hist, sum(hist.values()))
    # A base is "low" if depth < max(low_depth, low_frac * median depth)
    thr = max(args.low_depth, args.low_frac * median)
    print(f"median depth={median}  low-coverage threshold=<{thr:.1f}")

    for p in (args.out_summary, args.out_low):
        Path(p).parent.mkdir(parents=True, exist_ok=True)

    # ---- summary table
    with open(args.out_summary, "w") as out:
        out.write("#contig\tlength\tmean_depth\tmedian_depth\tpct_zero\tpct_low\n")
        for contig, arr in depth.items():
            n, h = len(arr), Counter(arr)
            low = sum(v for k, v in h.items() if k < thr)
            out.write(f"{contig}\t{n}\t{sum(arr)/n:.1f}\t{_median_from_hist(h, n)}\t"
                      f"{100*h[0]/n:.2f}\t{100*low/n:.2f}\n")

    # ---- pass 2: low regions, merging gaps <= low_merge_gap
    with open(args.out_low, "w") as out:
        out.write("#contig\tstart\tend\tlength\tmean_depth\tmin_depth\n")
        for contig, arr in depth.items():
            for start, end, mean, mn in low_regions(arr, thr, args.low_merge_gap, args.low_min_len):
                out.write(f"{contig}\t{start}\t{end}\t{end-start}\t{mean:.1f}\t{mn}\n")
    print(f"-> {args.out_low}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="sorted, indexed BAM")
    ap.add_argument("--fai", required=True, help="reference .fai (genome length)")
    ap.add_argument("--out-summary", required=True, help="coverage summary TSV")
    ap.add_argument("--out-low", required=True, help="low-coverage regions TSV")
    ap.add_argument("--min-mapq", type=int, default=20, help="minimum mapping quality [20]")
    ap.add_argument("--low-depth", type=int, default=10, help="absolute depth floor for 'low' [10]")
    ap.add_argument("--low-frac", type=float, default=0.25,
                    help="fraction of median depth below which a base is 'low' [0.25]")
    ap.add_argument("--low-min-len", type=int, default=20,
                    help="minimum low-coverage region length (bp) [20]")
    ap.add_argument("--low-merge-gap", type=int, default=10,
                    help="merge low regions closer than this (bp) [10]")
    args = ap.parse_args()
    try:
        run(args)
    except (OSError, ValueError) as e:
        sys.exit(f"lowcov: error: {e}")


if __name__ == "__main__":
    main()
