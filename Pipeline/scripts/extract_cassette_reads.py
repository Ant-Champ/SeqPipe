#!/usr/bin/env python3
"""
Extract the original FASTQ reads that align to the cassette.

Read names are collected from a BAM already restricted to the cassette
reference (primary/supplementary, mapped alignments; secondary ones are
skipped so a read is never duplicated). The original FASTQ(s) are then
filtered by those names.

Long reads : --reads R.fastq.gz            --out O.fastq.gz
Paired     : --reads R1.fastq.gz R2.fastq.gz --out O_R1.fastq.gz O_R2.fastq.gz
             Both mates of every read name that hits the cassette are kept
             (a pair is kept when either mate aligned), so the two output
             files stay in sync.

A trailing /1 or /2 on read names is ignored when matching, because
aligners strip it from BAM query names.
"""
import argparse
import gzip

import pysam


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def norm_name(raw):
    """First whitespace-delimited token, without a leading '@' or trailing /1 /2."""
    name = raw.strip().split()[0] if raw.strip() else ""
    if name.startswith("@"):
        name = name[1:]
    if len(name) > 2 and name[-2] == "/" and name[-1] in "12":
        name = name[:-2]
    return name


def bam_names(bam_path):
    """Query names of mapped, non-secondary alignments in the BAM."""
    names = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary:
                continue
            names.add(norm_name(read.query_name))
    return names


def filter_fastq(path_in, path_out, names):
    """Write FASTQ records of path_in whose name is in names. Returns (n_in, n_kept)."""
    n = kept = 0
    with open_text(path_in) as fh, gzip.open(path_out, "wt") as out:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                raise RuntimeError(f"Malformed FASTQ (incomplete record): {path_in}")
            n += 1
            if norm_name(h) in names:
                out.write(h + seq + plus + qual)
                kept += 1
    return n, kept


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True, help="BAM restricted to the cassette reference")
    p.add_argument("--reads", required=True, nargs="+", metavar="FASTQ",
                   help="input FASTQ (one file, or R1 R2)")
    p.add_argument("--out", required=True, nargs="+", metavar="FASTQ",
                   help="output FASTQ.gz (same number of files as --reads)")
    args = p.parse_args()

    if len(args.reads) not in (1, 2) or len(args.reads) != len(args.out):
        p.error("--reads and --out must each take the same number of files (1 or 2)")

    names = bam_names(args.bam)
    print(f"Cassette-aligned read names: {len(names)}")

    for r_in, r_out in zip(args.reads, args.out):
        n, kept = filter_fastq(r_in, r_out, names)
        print(f"{r_in}: {n} reads, extracted {kept} -> {r_out}")


if __name__ == "__main__":
    main()
