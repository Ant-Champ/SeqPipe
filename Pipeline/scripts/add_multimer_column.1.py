#!/usr/bin/env python3
"""
add_multimer_column.1.py -- append the sample-level multimer verdict as a
column to a per-variant TSV.

The multimer call describes the *molecule*, not any individual variant: every
variant gets the verdict of its own replicon (CHROM column; the sample-level
call when the table has no CHROM column or the replicon is not in the summary).
Its purpose is to make the final table self-contained: when tables from
several samples are concatenated, each variant still knows whether it came
from a multimeric plasmid.

Inputs
  -in       per-variant TSV (e.g. <sample>.annotated.gatc.simple.tsv)
  -summary  TSV written by multimer.2.py -summary (one row per replicon)
Output
  -out      same table with the extra column(s) appended

If the variant table has no data rows the header still gains the column, so
the schema stays identical across samples.
"""

import argparse
import csv
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description="Append the multimer yes/no verdict to a variant TSV.")
    p.add_argument("-in", dest="infile", required=True,
                   help="Input per-variant TSV")
    p.add_argument("-summary", required=True,
                   help="Summary TSV from multimer.2.py -summary")
    p.add_argument("-out", required=True, help="Output TSV")
    p.add_argument("-col", default="multimer",
                   help="Name of the added column [multimer]")
    p.add_argument("-extra", action="store_true",
                   help="Also append multimer_copies and multimer_pct")
    p.add_argument("-no_header", action="store_true",
                   help="Input table has no header line; the value is simply "
                        "appended to every row")
    p.add_argument("-default", default="NA",
                   help="Value used if the summary is missing or unreadable [NA]")
    return p.parse_args()


def read_summary(path, default):
    """Return the summary rows (one per replicon) as dicts, or [] on any problem."""
    try:
        with open(path, newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    except (OSError, csv.Error) as err:
        sys.stderr.write("WARNING: cannot read %s (%s); using '%s'\n"
                         % (path, err, default))
        return []


def _pct(r):
    try:
        return float(r.get("multimer_pct", 0))
    except ValueError:
        return 0.0


def sample_row(rows):
    """The row behind the sample-level call: a multimeric replicon if any,
    else one that could be judged, else the one with the most reads."""
    for call in ("yes", "no"):
        hits = [r for r in rows if r.get("multimer") == call]
        if hits:
            return max(hits, key=_pct)
    return max(rows, key=lambda r: int(r.get("reads_mapped") or 0)) if rows else {}


def values(row, args):
    vals = [row.get("multimer", args.default) or args.default]
    if args.extra:
        vals += [row.get("multimer_copies", args.default),
                 row.get("multimer_pct", args.default)]
    return vals


def main():
    args = parse_args()
    rows = read_summary(args.summary, args.default)
    by_replicon = {r.get("reference"): r for r in rows}
    fallback = sample_row(rows)

    names = [args.col] + (["multimer_copies", "multimer_pct"] if args.extra else [])

    n_rows = 0
    with open(args.infile, newline="") as fin, \
         open(args.out, "w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")

        first, chrom_col = True, None
        for row in reader:
            if not row:                     # preserve nothing, skip blank lines
                continue
            if first and not args.no_header:
                writer.writerow(row + names)
                # each variant gets the verdict of its own replicon
                chrom_col = row.index("CHROM") if "CHROM" in row else None
                first = False
                continue
            first = False
            src = fallback
            if chrom_col is not None and chrom_col < len(row):
                src = by_replicon.get(row[chrom_col], fallback)
            writer.writerow(row + values(src, args))
            n_rows += 1

        if first:                           # completely empty input file
            writer.writerow(names)

    sys.stderr.write("added %s (per replicon: %s) to %d row(s) -> %s\n"
                     % (args.col, ", ".join("%s=%s" % (r.get("reference"), r.get("multimer"))
                                            for r in rows) or "none", n_rows, args.out))


if __name__ == "__main__":
    main()
