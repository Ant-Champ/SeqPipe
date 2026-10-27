#!/usr/bin/env python3
"""
add_multimer_column.1.py -- append the sample-level multimer verdict as a
column to a per-variant TSV.

The multimer call describes the *molecule*, not any individual variant, so
the added column carries the same value on every row of a given sample.
Its purpose is to make the final table self-contained: when tables from
several samples are concatenated, each variant still knows whether it came
from a multimeric plasmid.

Inputs
  -in       per-variant TSV (e.g. <sample>.annotated.gatc.simple.tsv)
  -summary  one-line TSV written by multimer.1.py -summary
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
                   help="Summary TSV from multimer.1.py -summary")
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
    """Return the single summary row as a dict, or {} on any problem."""
    try:
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        return rows[0] if rows else {}
    except (OSError, csv.Error) as err:
        sys.stderr.write("WARNING: cannot read %s (%s); using '%s'\n"
                         % (path, err, default))
        return {}


def main():
    args = parse_args()
    summary = read_summary(args.summary, args.default)

    added = [(args.col, summary.get("multimer", args.default) or args.default)]
    if args.extra:
        added += [
            ("multimer_copies", summary.get("multimer_copies", args.default)),
            ("multimer_pct",    summary.get("multimer_pct", args.default)),
        ]

    n_rows = 0
    with open(args.infile, newline="") as fin, \
         open(args.out, "w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")

        first = True
        for row in reader:
            if not row:                     # preserve nothing, skip blank lines
                continue
            if first and not args.no_header:
                writer.writerow(row + [name for name, _ in added])
                first = False
                continue
            first = False
            writer.writerow(row + [val for _, val in added])
            n_rows += 1

        if first:                           # completely empty input file
            writer.writerow([name for name, _ in added])

    sys.stderr.write("added %s=%s to %d row(s) -> %s\n"
                     % (args.col, added[0][1], n_rows, args.out))


if __name__ == "__main__":
    main()
