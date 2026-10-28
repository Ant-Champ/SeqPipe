#!/usr/bin/env python3
"""mappability -- regions of a reference where short paired reads cannot be
placed uniquely (repeats: rRNA operons, IS elements, tRNA clusters...).

Every --step bp, a perfect read pair is cut from the reference (R1 forward at
i, R2 reverse-complement ending at i + --insert; reads wrap around the origin
of circular records, as real reads do), aligned back with bwa-mem2 exactly as
the sample reads are, and counted per base like `samtools depth -Q`
(--min-mapq; unmapped / secondary alignments skipped). A perfect pair gives
2 x read_len / step reads per base; a base is "unmappable" when it gets less
than --min-frac of that. Real reads behave the same there, so the low-coverage
scan (lowcov.py --mask) ignores these bases.

--out-bed : contig  start  end  (0-based half-open, merged, '#' header)
Needs bwa-mem2 on PATH. Python >= 3.8, biopython.
"""
import argparse
import gzip
import os
import re
import subprocess
import sys
import tempfile
from array import array

from Bio import SeqIO

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def simulate(records, read_len, insert, step, r1_path, r2_path):
    qual = "I" * read_len
    n = 0
    with gzip.open(r1_path, "wt", compresslevel=1) as f1, gzip.open(r2_path, "wt", compresslevel=1) as f2:
        for rec in records:
            seq = str(rec.seq).upper()
            L = len(seq)
            circular = rec.annotations.get("topology", "").lower() == "circular"
            if circular:
                seq = seq + seq[:insert]
                last = L                       # every start position
            else:
                last = L - insert + 1          # pairs must fit in the record
            if insert > L:
                print(f"WARNING: {rec.id} ({L} bp) is shorter than the insert: skipped", file=sys.stderr)
                continue
            for i in range(0, max(last, 0), step):
                r1 = seq[i:i + read_len]
                r2 = seq[i + insert - read_len:i + insert].translate(_COMP)[::-1]
                name = f"{rec.id}_{i}"
                f1.write(f"@{name}/1\n{r1}\n+\n{qual}\n")
                f2.write(f"@{name}/2\n{r2}\n+\n{qual}\n")
                n += 1
    return n


def count_depth(sam, lengths, min_mapq):
    """Per-base depth of the alignments of a SAM stream (difference arrays)."""
    diff = {c: array("i", bytes(4 * (L + 1))) for c, L in lengths.items()}
    for line in sam:
        if line.startswith("@"):
            continue
        f = line.split("\t", 6)
        flag, contig, pos, mapq, cigar = int(f[1]), f[2], int(f[3]) - 1, int(f[4]), f[5]
        if flag & 0x704 or mapq < min_mapq or contig not in diff:   # as samtools depth: unmapped, secondary, qcfail, dup
            continue
        span = sum(int(n) for n, op in _CIGAR.findall(cigar) if op in "MDN=X")
        d, L = diff[contig], lengths[contig]
        end = pos + span
        if end <= L:
            d[pos] += 1
            d[end] -= 1
        else:                                  # alignment over the origin of a circular record
            d[pos] += 1
            d[L] -= 1
            d[0] += 1
            d[end - L] -= 1
    depth = {}
    for c, d in diff.items():
        run, out = 0, array("i", bytes(4 * lengths[c]))
        for i in range(lengths[c]):
            run += d[i]
            out[i] = run
        depth[c] = out
    return depth


def unmappable_regions(arr, thr):
    regions, start = [], None
    for i, v in enumerate(arr):
        if v < thr:
            if start is None:
                start = i
        elif start is not None:
            regions.append((start, i))
            start = None
    if start is not None:
        regions.append((start, len(arr)))
    return regions


def run(args):
    records = list(SeqIO.parse(args.genbank, "genbank"))
    if not records:
        sys.exit(f"mappability: no GenBank record in {args.genbank}")
    lengths = {r.id: len(r.seq) for r in records}
    os.makedirs(os.path.dirname(os.path.abspath(args.out_bed)), exist_ok=True)

    with tempfile.TemporaryDirectory(dir=args.tmpdir) as tmp:
        ref = os.path.join(tmp, "ref.fasta")
        SeqIO.write(records, ref, "fasta")
        r1, r2 = os.path.join(tmp, "r1.fq.gz"), os.path.join(tmp, "r2.fq.gz")
        n = simulate(records, args.read_len, args.insert, args.step, r1, r2)
        print(f"simulated {n:,} pairs ({args.read_len} bp, insert {args.insert}, step {args.step})", flush=True)

        subprocess.run(["bwa-mem2", "index", ref], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        proc = subprocess.Popen(["bwa-mem2", "mem", "-t", str(args.threads), ref, r1, r2],
                                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        depth = count_depth(proc.stdout, lengths, args.min_mapq)
        if proc.wait() != 0:
            sys.exit("mappability: bwa-mem2 mem failed")

    expected = 2 * args.read_len / args.step
    thr = args.min_frac * expected
    total_bad = 0
    with open(args.out_bed, "w") as out:
        out.write(f"#contig\tstart\tend\t# PE{args.read_len} insert {args.insert} MAPQ>={args.min_mapq}: "
                  f"depth < {thr:.1f} of {expected:.0f} expected\n")
        for rec in records:
            regs = unmappable_regions(depth[rec.id], thr)
            bad = sum(e - s for s, e in regs)
            total_bad += bad
            print(f"{rec.id}: {len(regs)} unmappable region(s), {bad:,} bp "
                  f"({100 * bad / max(lengths[rec.id], 1):.2f}%)")
            for s, e in regs:
                out.write(f"{rec.id}\t{s}\t{e}\n")
    print(f"total unmappable: {total_bad:,} bp -> {args.out_bed}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genbank", required=True, help="reference GenBank file (all records)")
    ap.add_argument("--out-bed", required=True, help="unmappable regions (BED)")
    ap.add_argument("--read-len", type=int, default=75, help="read length (bp) [75]")
    ap.add_argument("--insert", type=int, default=300, help="fragment length (bp) [300]")
    ap.add_argument("--step", type=int, default=2, help="one pair every STEP bp [2]")
    ap.add_argument("--min-mapq", type=int, default=20, help="as in the low-coverage scan [20]")
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="unmappable below this fraction of the expected depth [0.5]")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--tmpdir", default=os.environ.get("TMPDIR"), help="temporary folder [$TMPDIR]")
    args = ap.parse_args()
    if args.insert < args.read_len:
        sys.exit("mappability: --insert must be >= --read-len")
    run(args)


if __name__ == "__main__":
    main()
