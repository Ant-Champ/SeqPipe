#!/usr/bin/env python3
"""
aggregate_coverage.1.py -- concatenate per-sample `samtools coverage`
tables into one run-level long-format table, with a sample column
prepended. Shared by two branches:
  * plasmid-panel screen : one row per sample x candidate plasmid
  * reference coverage   : one row per sample x replicon in that
                            sample's own ref.fasta (chromosome and/or
                            plasmid records)

Both branches use the same -name_col (default "replicon") and pass their
own -source tag, so their outputs share one schema and can be
concatenated downstream with no further reconciling (see
merge_coverage_summaries in the Snakefile).

Input   : one `samtools coverage` TSV per sample (-in), plus the matching
          sample name for each, in the same order (-samples)
Output  : one TSV (-out), one row per sample x reference sequence,
          columns: sample, [source,] <name_col>, startpos, endpos,
          numreads, covbases, coverage, meandepth, meanbaseq, meanmapq
"""

import argparse
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Aggregate per-sample `samtools coverage` tables into one run-level table.")
    p.add_argument("-in", dest="infiles", nargs="*", default=[],
                   help="Per-sample coverage TSVs from `samtools coverage`")
    p.add_argument("-samples", nargs="*", default=[],
                   help="Sample name for each -in file, same order/count")
    p.add_argument("-out", required=True, help="Output aggregated TSV")
    p.add_argument("-name_col", default="replicon",
                   help="Column name for the reference-sequence identity "
                        "(samtools coverage's #rname) [replicon]")
    p.add_argument("-source", default=None,
                   help="If set, adds a constant 'source' column with this "
                        "value (e.g. 'reference' or 'plasmid_panel'), so rows "
                        "stay distinguishable after merging with another "
                        "aggregated table [unset]")
    return p.parse_args()


def main():
    args = parse_args()
    if len(args.infiles) != len(args.samples):
        sys.exit("ERROR: -in and -samples must have the same number of entries "
                  f"({len(args.infiles)} vs {len(args.samples)})")

    if not args.infiles:
        # every sample was skipped (e.g. all low_diversity): an empty table, not an error
        cols = ["sample"] + (["source"] if args.source else []) + [
            args.name_col, "startpos", "endpos", "numreads", "covbases",
            "coverage", "meandepth", "meanbaseq", "meanmapq"]
        pd.DataFrame(columns=cols).to_csv(args.out, sep="\t", index=False)
        sys.stderr.write(f"no samples to aggregate: wrote an empty {args.out}\n")
        return

    frames = []
    for path, sample in zip(args.infiles, args.samples):
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception as err:
            sys.stderr.write(f"WARNING: cannot read {path} ({err}); skipping\n")
            continue
        df = df.rename(columns={"#rname": args.name_col})
        df.insert(0, "sample", sample)
        if args.source:
            df.insert(1, "source", args.source)
        frames.append(df)

    if not frames:
        sys.exit("ERROR: no readable coverage tables, nothing to aggregate")

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(args.out, sep="\t", index=False)
    sys.stderr.write(f"wrote {args.out} ({len(combined)} sample x {args.name_col} rows)\n")


if __name__ == "__main__":
    main()
