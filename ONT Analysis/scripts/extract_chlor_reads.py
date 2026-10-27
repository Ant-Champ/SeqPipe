#!/usr/bin/env python3
import argparse
import gzip
import pysam


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required=True)
    p.add_argument("--reads", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    names = set()

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary:
                continue
            names.add(read.query_name)

    n = 0
    kept = 0

    with open_text(args.reads) as fh, gzip.open(args.out, "wt") as out:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()

            if not qual:
                raise RuntimeError("Malformed FASTQ: incomplete record.")

            n += 1
            name = h[1:].strip().split()[0] if h.startswith("@") else h.strip().split()[0]

            if name in names:
                out.write(h)
                out.write(seq)
                out.write(plus)
                out.write(qual)
                kept += 1

    print(f"Panel-aligned read names: {len(names)}")
    print(f"Input FASTQ reads: {n}")
    print(f"Extracted unique reads: {kept}")


if __name__ == "__main__":
    main()
