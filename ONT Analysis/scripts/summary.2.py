#!/usr/bin/env python3
"""
summary.2.py -- per-run summary table, one row per sample.

Changes from summary.1.py
  * added multimer / multimer_copies / multimer_pct columns, read from
    <sample>.multimer.summary.tsv
  * fixed single_cluster being written to the wrong row (used `idx`, the
    leftover loop variable from the read-counting loop, instead of `i`)
  * output path is an explicit argument so Snakemake controls it
  * bare `except:` replaced by `except Exception` with the reason reported
  * metrics accumulated as records instead of assigning into columns
    pre-filled with '' (that pattern raises TypeError on pandas >= 3.0)
  * total_bp counted with a plain 4-line FASTQ reader instead of SeqIO
    (same result, much faster on full ONT runs)
"""

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd
from Bio import SeqIO


# ----------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Generate ONT run summary")
    p.add_argument("-out_dir", required=True, help="Pipeline output directory")
    p.add_argument("-sample_sheet", required=True, help="Path to the sample sheet")
    p.add_argument("-fastq_dir", required=True, help="Raw fastq directory")
    p.add_argument("-out", default=None,
                   help="Summary output path [<out_dir>/summary.tsv]")
    p.add_argument("-fastq_base_name", default="",
                   help="Prefix before the barcode in raw fastq filenames "
                        "(e.g. 'SQK-RBK114-96_')")
    return p.parse_args()


def warn(msg):
    sys.stderr.write("WARNING: %s\n" % msg)


# ----------------------------------------------------------------------
def fastq_total_bp(path):
    """Total sequenced bases in a 4-line-per-record gzipped FASTQ."""
    total = 0
    with gzip.open(path, "rt") as fh:
        for n, line in enumerate(fh):
            if n % 4 == 1:                       # sequence line
                total += len(line.rstrip("\n"))
    return total


def mean_depth(path):
    """Mean depth from `samtools depth -a` output; 0 if the file is empty."""
    depths = []
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                depths.append(int(fields[2]))
    return round(float(np.mean(depths)), 2) if depths else 0


def read_multimer_summary(path):
    """One-row TSV from multimer.1.py -summary, as a dict ({} on failure)."""
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
        return df.iloc[0].to_dict() if len(df) else {}
    except Exception as err:
        warn("cannot read %s (%s)" % (path, err))
        return {}


def count_single_cluster(out_dir, sample):
    """True when trycycler produced exactly one cluster."""
    path = os.path.join(out_dir, sample, "trycycler")
    if not os.path.isdir(path):
        return False
    subdirs = [d for d in os.listdir(path)
               if os.path.isdir(os.path.join(path, d))]
    return len(subdirs) == 1


# ----------------------------------------------------------------------
def main():
    args = parse_args()
    outdir = os.path.join(args.out_dir, "")
    fastqdir = os.path.join(args.fastq_dir, "")
    out_path = args.out or os.path.join(outdir, "summary.tsv")

    sample_df = pd.read_csv(args.sample_sheet, sep="\t")

    # Metrics are accumulated as records and joined at the end.  Do NOT
    # pre-create columns with `sample_df[col] = ""` and assign into them:
    # that fixes the column dtype to str, and pandas >= 3.0 then rejects
    # every numeric assignment with a TypeError.
    records = []

    # ------------------------------------------------------------------
    # raw yield per barcode
    # ------------------------------------------------------------------
    total_bp = []
    for _, sample in sample_df.iterrows():
        fq = os.path.join(fastqdir,
                          args.fastq_base_name + str(sample["barcode"]) + ".fastq.gz")
        try:
            total_bp.append(fastq_total_bp(fq))
        except Exception as err:
            warn("reads issue with sample %s (%s): %s"
                 % (sample["sample_name"], fq, err))
            total_bp.append(0)

    # ------------------------------------------------------------------
    # per-sample metrics
    # ------------------------------------------------------------------
    for (_, row), bp in zip(sample_df.iterrows(), total_bp):
        sample = row["sample_name"]
        base = os.path.join(outdir, sample, sample)
        rec = {"total_bp": bp}

        # reference length ------------------------------------------------
        try:
            ref_len = len(list(SeqIO.parse(base + ".ref.fasta", "fasta"))[0].seq)
            rec["REF_length"] = ref_len
            rec["coverage_approx"] = round(bp / ref_len, 2) if ref_len else "ERR"
        except Exception as err:
            warn("no reference for %s (%s)" % (sample, err))
            rec["REF_length"] = "ERR"
            rec["coverage_approx"] = "ERR"

        # observed coverage -----------------------------------------------
        try:
            rec["coverage"] = mean_depth(base + ".cov.txt")
        except Exception as err:
            warn("no coverage for %s (%s)" % (sample, err))
            rec["coverage"] = "ERR"

        # variants ----------------------------------------------------------
        try:
            rec["detected_mutations"] = len(pd.read_csv(base + ".simple.tsv", sep="\t"))
        except Exception as err:
            warn("no variant table for %s (%s)" % (sample, err))
            rec["detected_mutations"] = "ERR"

        # coverage anomalies -------------------------------------------------
        try:
            with open(base + ".cov.log") as fh:
                rec["detected_cov_anomalies"] = len(fh.readlines())
        except Exception as err:
            warn("no coverage log for %s (%s)" % (sample, err))
            rec["detected_cov_anomalies"] = "ERR"

        # multimer ------------------------------------------------------------
        multi = read_multimer_summary(base + ".multimer.summary.tsv")
        rec["multimer"] = multi.get("multimer", "NA")
        rec["multimer_copies"] = multi.get("multimer_copies", "")
        rec["multimer_pct"] = multi.get("multimer_pct", "")

        # trycycler clusters ---------------------------------------------------
        rec["single_cluster"] = count_single_cluster(outdir, sample)

        records.append(rec)

    metrics = pd.DataFrame.from_records(records, index=sample_df.index)
    sample_df = pd.concat([sample_df, metrics], axis=1)

    sample_df.to_csv(out_path, sep="\t", index=False)
    sys.stderr.write("wrote %s (%d samples)\n" % (out_path, len(sample_df)))


if __name__ == "__main__":
    main()
