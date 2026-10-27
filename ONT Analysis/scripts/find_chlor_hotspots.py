#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
import pysam


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required=True)
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-bed", required=True)
    p.add_argument("--max-gap", type=int, default=1000)
    p.add_argument("--flank", type=int, default=500)
    p.add_argument("--min-reads", type=int, default=2)
    args = p.parse_args()

    # Collect one genomic interval per alignment record.
    # Supplementary alignments are retained; secondary alignments are
    # ignored because the panel alignment was run with --secondary=no.
    by_chrom = defaultdict(list)

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_id < 0:
                continue
            if read.is_secondary:
                continue
            chrom = bam.get_reference_name(read.reference_id)
            start = read.reference_start
            end = read.reference_end
            if end is None or end <= start:
                continue
            by_chrom[chrom].append((start, end, read.query_name))

    hotspots = []

    for chrom, intervals in by_chrom.items():
        intervals.sort(key=lambda x: (x[0], x[1]))

        cluster = []
        cluster_end = None

        def flush():
            if not cluster:
                return
            starts = [x[0] for x in cluster]
            ends = [x[1] for x in cluster]
            names = {x[2] for x in cluster}
            hotspots.append(
                {
                    "chrom": chrom,
                    "start": min(starts),
                    "end": max(ends),
                    "n_alignments": len(cluster),
                    "n_reads": len(names),
                }
            )

        for interval in intervals:
            start, end, name = interval

            if not cluster:
                cluster = [interval]
                cluster_end = end
                continue

            if start <= cluster_end + args.max_gap:
                cluster.append(interval)
                cluster_end = max(cluster_end, end)
            else:
                flush()
                cluster = [interval]
                cluster_end = end

        flush()

    # Keep only clusters supported by at least min_reads distinct reads.
    hotspots = [h for h in hotspots if h["n_reads"] >= args.min_reads]

    # Sort deterministically by chromosome and coordinate.
    hotspots.sort(key=lambda h: (h["chrom"], h["start"], h["end"]))

    # Assign IDs and expanded coordinates.
    for i, h in enumerate(hotspots, start=1):
        h["id"] = f"hotspot_{i:03d}"
        h["expanded_start"] = max(0, h["start"] - args.flank)
        h["expanded_end"] = h["end"] + args.flank

    with open(args.out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "hotspot",
            "chrom",
            "start",
            "end",
            "n_alignments",
            "n_reads",
            "expanded_start",
            "expanded_end",
        ])
        for h in hotspots:
            writer.writerow([
                h["id"],
                h["chrom"],
                h["start"],
                h["end"],
                h["n_alignments"],
                h["n_reads"],
                h["expanded_start"],
                h["expanded_end"],
            ])

    # BED is 0-based, half-open.
    with open(args.out_bed, "w") as fh:
        for h in hotspots:
            fh.write(
                f'{h["chrom"]}\t{h["expanded_start"]}\t{h["expanded_end"]}'
                f'\t{h["id"]}\t{h["n_reads"]}\n'
            )

    print(f"Found {len(hotspots)} Chlor_cassette hotspot(s).")
    for h in hotspots:
        print(
            f'{h["id"]}: {h["chrom"]}:{h["start"] + 1}-{h["end"]} '
            f'({h["n_reads"]} reads; '
            f'IGV {h["expanded_start"] + 1}-{h["expanded_end"]})'
        )


if __name__ == "__main__":
    main()
