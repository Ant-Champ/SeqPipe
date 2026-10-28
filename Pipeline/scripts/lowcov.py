#!/usr/bin/env python3
"""lowcov -- per-contig coverage summary + merged low- and high-coverage regions.

Ported from alignqc.py (step_lowcov). Depth comes from `samtools depth -aa -Q
<min-mapq>` (every position, including zeros and reference contigs with no
reads) and is held in memory at 4 bytes/base: fine for bacteria/plasmids/phage.
A base is "low" if its depth, averaged over --low-smooth bp around it, is below
max(low_depth, low_frac * median depth), capped at --low-max-thr-frac * median
(so an absolute floor close to the median cannot turn ordinary Poisson dips of
a low-depth sample into low regions). Low bases closer than --low-merge-gap
(not counting masked bases) are merged; regions shorter than --low-min-len, or
whose mean depth is above --low-region-max-frac * median, are dropped (a real
deletion is ~0x over its whole length; a cluster of dips is not).

--mask (BED, e.g. the unmappable regions of mappability.py): masked bases are
never low and are left out of the median; a low region can run across them.

High coverage (duplications, amplifications): the replicon is cut into
--high-window bp windows; a window is "high" if its mean depth is above
--high-frac x the median depth OF THAT REPLICON (so a multicopy plasmid is not
flagged as a whole: its copy number is in the summary). High windows closer
than --high-merge-gap are merged and regions shorter than --high-min-len are
dropped. fold = region mean / replicon median (~2 for one extra copy).

Coordinates are 0-based half-open (BED-like).

--out-summary columns ('#' header line):
  #contig  length  mean_depth  median_depth  pct_zero  pct_low  pct_masked
--out-low columns ('#' header line):
  #contig  start  end  length  mean_depth  min_depth  masked_bp  rel_depth
  (mean_depth / min_depth over the unmasked bases; rel_depth = mean / median)
--out-high columns ('#' header line):
  #contig  start  end  length  mean_depth  replicon_median  fold

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


def smooth(arr, window):
    """Centered rolling mean of arr over `window` bp (window <= 1: arr itself)."""
    if window <= 1:
        return arr
    n, h = len(arr), window // 2
    pre = array("q", [0]) * (n + 1)
    run = 0
    for i, d in enumerate(arr):
        run += d
        pre[i + 1] = run
    out = array("f", bytes(4 * n))
    for i in range(n):
        a, b = max(0, i - h), min(n, i + h + 1)
        out[i] = (pre[b] - pre[a]) / (b - a)
    return out


def low_regions(arr, sm, mask, thr, merge_gap, min_len):
    """Merged low-depth regions of one contig -> [(start, end, mean, min, masked_bp)].
    sm = smoothed depth (decides "low"), arr = raw depth (mean/min reported),
    mask = bytearray (1 = masked) or None. Gaps are counted in unmasked bases."""
    regions = []                                # [start, end, sum, n, min]
    gap = None                                  # unmasked non-low bases since the last low base
    for i in range(len(arr)):
        if mask is not None and mask[i]:
            continue
        if sm[i] >= thr:
            if gap is not None:
                gap += 1
            continue
        d = arr[i]
        if regions and gap is not None and gap <= merge_gap:
            r = regions[-1]
            r[1], r[2], r[3], r[4] = i + 1, r[2] + d, r[3] + 1, min(r[4], d)
        else:
            regions.append([i, i + 1, d, 1, d])
        gap = 0
    out = []
    for s, e, tot, n, mn in regions:
        if e - s < min_len:
            continue
        masked = sum(mask[s:e]) if mask is not None else 0
        # mean / min over the unmasked bases of the whole region (gaps included)
        vals = [arr[j] for j in range(s, e) if mask is None or not mask[j]]
        out.append((s, e, sum(vals) / len(vals), min(vals), masked))
    return out


def read_mask(path, lengths):
    """BED -> {contig: bytearray (1 = masked)}."""
    mask = {}
    if not path:
        return mask
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c, s, e = line.split("\t")[:3]
            if c not in lengths:
                continue
            m = mask.setdefault(c, bytearray(lengths[c]))
            s, e = max(0, int(s)), min(lengths[c], int(e))
            m[s:e] = b"\x01" * (e - s)
    return mask


def high_regions(arr, median, frac, window, merge_gap, min_len):
    """Merged high-depth windows of one contig -> [(start, end, mean)]."""
    if median <= 0:
        return []
    thr, regions = frac * median, []            # [start, end]
    for ws in range(0, len(arr), window):
        w = arr[ws:ws + window]
        if sum(w) / len(w) <= thr:
            continue
        if regions and ws - regions[-1][1] <= merge_gap:
            regions[-1][1] = ws + len(w)
        else:
            regions.append([ws, ws + len(w)])
    return [(s, e, sum(arr[s:e]) / (e - s)) for s, e in regions if e - s >= min_len]


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

    mask = read_mask(args.mask, {c: len(a) for c, a in depth.items()})
    if args.mask:
        n_masked = sum(sum(m) for m in mask.values())
        print(f"mask {args.mask}: {n_masked:,} bp ignored ({100 * n_masked / max(genome_len, 1):.2f}%)")

    hist = Counter()
    for c, arr in depth.items():
        m = mask.get(c)
        hist.update(arr if m is None else (d for d, x in zip(arr, m) if not x))
    median = _median_from_hist(hist, sum(hist.values())) or 0
    # A base is "low" if its (smoothed) depth < max(low_depth, low_frac * median),
    # never above low_max_thr_frac * median
    thr = max(args.low_depth, args.low_frac * median)
    cap = args.low_max_thr_frac * median
    if thr > cap:
        print(f"WARNING: low-coverage threshold {thr:.1f} is above {args.low_max_thr_frac} x median "
              f"({median}): capped at {cap:.1f}")
        thr = cap
    region_max = args.low_region_max_frac * median
    print(f"median depth={median}  low-coverage threshold=<{thr:.1f} (smoothed over {args.low_smooth} bp)  "
          f"region mean must be <={region_max:.1f}")

    for p in (args.out_summary, args.out_low, args.out_high):
        Path(p).parent.mkdir(parents=True, exist_ok=True)

    # ---- summary table
    contig_median = {}
    with open(args.out_summary, "w") as out:
        out.write("#contig\tlength\tmean_depth\tmedian_depth\tpct_zero\tpct_low\tpct_masked\n")
        for contig, arr in depth.items():
            n, h = len(arr), Counter(arr)
            low = sum(v for k, v in h.items() if k < thr)
            contig_median[contig] = _median_from_hist(h, n)
            masked = sum(mask[contig]) if contig in mask else 0
            out.write(f"{contig}\t{n}\t{sum(arr)/n:.1f}\t{contig_median[contig]}\t"
                      f"{100*h[0]/n:.2f}\t{100*low/n:.2f}\t{100*masked/n:.2f}\n")

    # ---- pass 2: low regions, merging gaps <= low_merge_gap
    with open(args.out_low, "w") as out:
        out.write("#contig\tstart\tend\tlength\tmean_depth\tmin_depth\tmasked_bp\trel_depth\n")
        dropped = 0
        for contig, arr in depth.items():
            sm = smooth(arr, args.low_smooth)
            for start, end, mean, mn, masked in low_regions(arr, sm, mask.get(contig), thr,
                                                            args.low_merge_gap, args.low_min_len):
                if mean > region_max:
                    dropped += 1
                    continue
                rel = mean / median if median else 0
                out.write(f"{contig}\t{start}\t{end}\t{end-start}\t{mean:.1f}\t{mn}\t{masked}\t{rel:.2f}\n")
        print(f"{dropped} candidate low region(s) dropped: mean depth above {region_max:.1f}")
    print(f"-> {args.out_low}")

    # ---- pass 3: high regions, against each replicon's own median
    with open(args.out_high, "w") as out:
        out.write("#contig\tstart\tend\tlength\tmean_depth\treplicon_median\tfold\n")
        for contig, arr in depth.items():
            med = contig_median[contig]
            print(f"{contig}: median depth={med}  high-coverage threshold=>{args.high_frac * med:.1f}")
            for start, end, mean in high_regions(arr, med, args.high_frac, args.high_window,
                                                 args.high_merge_gap, args.high_min_len):
                out.write(f"{contig}\t{start}\t{end}\t{end-start}\t{mean:.1f}\t{med}\t{mean/med:.2f}\n")
    print(f"-> {args.out_high}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="sorted, indexed BAM")
    ap.add_argument("--fai", required=True, help="reference .fai (genome length)")
    ap.add_argument("--out-summary", required=True, help="coverage summary TSV")
    ap.add_argument("--out-low", required=True, help="low-coverage regions TSV")
    ap.add_argument("--out-high", required=True, help="high-coverage regions TSV")
    ap.add_argument("--min-mapq", type=int, default=20, help="minimum mapping quality [20]")
    ap.add_argument("--low-depth", type=int, default=10, help="absolute depth floor for 'low' [10]")
    ap.add_argument("--low-frac", type=float, default=0.25,
                    help="fraction of median depth below which a base is 'low' [0.25]")
    ap.add_argument("--low-min-len", type=int, default=20,
                    help="minimum low-coverage region length (bp) [20]")
    ap.add_argument("--low-merge-gap", type=int, default=10,
                    help="merge low regions closer than this (bp) [10]")
    ap.add_argument("--low-smooth", type=int, default=50,
                    help="average the depth over this many bp before calling 'low' (1 = off) [50]")
    ap.add_argument("--low-max-thr-frac", type=float, default=0.5,
                    help="the low threshold never exceeds this x the median depth [0.5]")
    ap.add_argument("--low-region-max-frac", type=float, default=0.25,
                    help="drop a low region whose mean depth is above this x the median [0.25]")
    ap.add_argument("--mask", default=None,
                    help="BED of bases to ignore (e.g. unmappable regions) [none]")
    ap.add_argument("--high-frac", type=float, default=1.7,
                    help="window is 'high' above this x the replicon median depth [1.7]")
    ap.add_argument("--high-window", type=int, default=100,
                    help="window size (bp) for the high-coverage scan [100]")
    ap.add_argument("--high-min-len", type=int, default=1000,
                    help="minimum high-coverage region length (bp) [1000]")
    ap.add_argument("--high-merge-gap", type=int, default=300,
                    help="merge high regions closer than this (bp) [300]")
    args = ap.parse_args()
    try:
        run(args)
    except (OSError, ValueError) as e:
        sys.exit(f"lowcov: error: {e}")


if __name__ == "__main__":
    main()
