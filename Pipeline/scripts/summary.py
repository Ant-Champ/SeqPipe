#!/usr/bin/env python3
"""
summary.py -- run summary, one row per sample (sheet order).  Used through
Snakemake's `script:` directive (rule run_summary, rules/summary.smk).

snakemake.params : samples (list of dicts), out_dir, toggles (dict)
snakemake.output : summary.tsv

Only small per-sample result files are read (seqkit stats table, flagstat,
samtools coverage, result tables); never the FASTQs.

Cell values: 'NA' = analysis off / not applicable, 'ERR' = an expected file
could not be parsed (the exception is logged to stderr).
Pure python (csv + re); no pandas needed.
"""

import csv
import os
import re
import sys

NA, ERR = "NA", "ERR"

COLUMNS = [
    "sample_name", "platform", "reference",
    "raw_reads", "raw_bases", "filt_reads", "filt_bases", "filt_N50",
    "ref_length", "ref_records", "coverage_approx",
    "pct_mapped", "mean_depth", "min_replicon_breadth",
    "detected_mutations", "low_cov_regions", "large_indels", "breakpoints",
    "multimer", "multimer_copies", "multimer_pct",
    "assembly_contigs", "assembly_circular", "cassette_insertions",
]


def warn(sample, what, err):
    sys.stderr.write("WARNING: sample %s: cannot read %s: %s: %s\n"
                     % (sample, what, type(err).__name__, err))


# ----------------------------------------------------------------------
# parsers
# ----------------------------------------------------------------------
def to_int(x):
    return int(str(x).replace(",", "").strip())


def read_stats(path, raw_names, filt_names, is_ont):
    """Parse `seqkit stats -a -T` (raw files first, then filtered files).

    Returns dict raw_reads raw_bases filt_reads filt_bases filt_N50.
    Rows are matched by file basename; mates are summed. N50 for ONT only.
    """
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")]
    by_name = {os.path.basename(r["file"]): r for r in rows}

    def pick(names):
        out = []
        for n in names:
            if n not in by_name:
                raise KeyError("no row for %s in %s" % (n, path))
            out.append(by_name[n])
        return out

    raw, filt = pick(raw_names), pick(filt_names)
    res = {
        "raw_reads":  sum(to_int(r["num_seqs"]) for r in raw),
        "raw_bases":  sum(to_int(r["sum_len"]) for r in raw),
        "filt_reads": sum(to_int(r["num_seqs"]) for r in filt),
        "filt_bases": sum(to_int(r["sum_len"]) for r in filt),
        "filt_N50":   NA,
    }
    if is_ont:
        res["filt_N50"] = to_int(filt[0]["N50"])
    return res


def fasta_stats(path):
    """(total length of all records, number of records) of a FASTA."""
    total, n = 0, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
            else:
                total += len(line.strip())
    return total, n


_MAPPED_RE = re.compile(r"^\d+\s*\+\s*\d+\s+(primary mapped|mapped)\s+\(([\d.]+)%")


def flagstat_pct_mapped(path):
    """'primary mapped (xx%' if present, else 'mapped (xx%'; as a string."""
    found = {}
    with open(path) as fh:
        for line in fh:
            m = _MAPPED_RE.match(line.strip())
            if m:
                found[m.group(1)] = m.group(2)
    if "primary mapped" in found:
        return found["primary mapped"]
    if "mapped" in found:
        return found["mapped"]
    return NA          # e.g. 0 reads -> 'mapped (N/A : N/A)'


def coverage_stats(path):
    """samtools coverage -> (length-weighted mean depth, min breadth %)."""
    wsum, length, min_cov = 0.0, 0, None
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        raise ValueError("empty file")
    header = [h.lstrip("#") for h in rows[0]]
    idx = {h: i for i, h in enumerate(header)}
    for r in rows[1:]:
        if not r or not r[0].strip():
            continue
        ln = int(r[idx["endpos"]]) - int(r[idx["startpos"]]) + 1
        wsum += float(r[idx["meandepth"]]) * ln
        length += ln
        c = float(r[idx["coverage"]])
        min_cov = c if min_cov is None else min(min_cov, c)
    if length == 0:
        raise ValueError("no replicons")
    return round(wsum / length, 2), round(min_cov, 2)


def count_rows(path):
    """Data rows of a TSV.  Blank and '#' lines are skipped; the first
    remaining line is a column header unless a '#' header line was seen."""
    n, saw_hash = 0, False
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                saw_hash = True
                continue
            n += 1
    if not saw_hash and n > 0:
        n -= 1
    return n


def multimer_summary(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError("no data row")
    r = rows[0]
    return r.get("multimer", ERR), r.get("multimer_copies", ERR), r.get("multimer_pct", ERR)


def assembly_stats(path):
    """Flye assembly_info.txt -> (contigs, circular contigs) or 'skipped'."""
    contigs, circ = 0, 0
    with open(path) as fh:
        first = fh.readline()
        if first.lstrip().startswith("# SKIPPED"):
            return "skipped", "skipped"
        fh.seek(0)
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            contigs += 1
            if len(f) > 3 and f[3].strip().upper() in ("Y", "+", "YES", "TRUE"):
                circ += 1
    return contigs, circ


# ----------------------------------------------------------------------
def try_(sample, what, fn, default=ERR):
    try:
        return fn()
    except Exception as err:
        warn(sample, what, err)
        return default


def summarize_sample(s, out_dir, tg):
    name, platform = s["sample_name"], s["platform"]
    ont = platform == "ont"
    d = os.path.join(out_dir, name)
    base = os.path.join(d, name)
    rec = {c: NA for c in COLUMNS}
    rec.update(sample_name=name, platform=platform, reference=s.get("reference", NA))

    # read stats ---------------------------------------------------------
    if ont:
        raw_n, filt_n = ["%s.fastq.gz" % name], ["%s.filt.fastq.gz" % name]
    else:
        raw_n = ["%s_R1.fastq.gz" % name, "%s_R2.fastq.gz" % name]
        filt_n = ["%s.filt_R1.fastq.gz" % name, "%s.filt_R2.fastq.gz" % name]
    st = try_(name, "read_stats", lambda: read_stats(
        os.path.join(d, "qc", name + ".read_stats.tsv"), raw_n, filt_n, ont), default=None)
    if st is None:
        for k in ("raw_reads", "raw_bases", "filt_reads", "filt_bases", "filt_N50"):
            rec[k] = ERR
    else:
        rec.update(st)

    # reference ----------------------------------------------------------
    fa = try_(name, "ref.fasta", lambda: fasta_stats(base + ".ref.fasta"), default=None)
    if fa is None:
        rec["ref_length"] = rec["ref_records"] = rec["coverage_approx"] = ERR
    else:
        rec["ref_length"], rec["ref_records"] = fa
        if st is None or fa[0] == 0:
            rec["coverage_approx"] = ERR
        else:
            rec["coverage_approx"] = "%.1f" % (st["filt_bases"] / fa[0])

    # alignment ----------------------------------------------------------
    rec["pct_mapped"] = try_(name, "flagstat", lambda: flagstat_pct_mapped(base + ".flagstat.txt"))
    cov = try_(name, "ref_coverage", lambda: coverage_stats(base + ".ref_coverage.tsv"), default=None)
    if cov is None:
        rec["mean_depth"] = rec["min_replicon_breadth"] = ERR
    else:
        rec["mean_depth"], rec["min_replicon_breadth"] = cov

    # optional analyses --------------------------------------------------
    if tg.get("variants"):
        rec["detected_mutations"] = try_(name, "simple.tsv", lambda: count_rows(base + ".simple.tsv"))
    if tg.get("coverage"):
        rec["low_cov_regions"] = try_(name, "low_coverage table", lambda: count_rows(
            os.path.join(d, "annotated", name + ".low_coverage.annotated.tsv")))
    if tg.get("sv"):
        rec["large_indels"] = try_(name, "indels table", lambda: count_rows(
            os.path.join(d, "annotated", name + ".indels.annotated.tsv")))
        rec["breakpoints"] = try_(name, "breakpoints table", lambda: count_rows(
            os.path.join(d, "annotated", name + ".breakpoints.annotated.tsv")))
    if tg.get("multimer") and ont:
        m = try_(name, "multimer summary", lambda: multimer_summary(
            base + ".multimer.summary.tsv"), default=(ERR, ERR, ERR))
        rec["multimer"], rec["multimer_copies"], rec["multimer_pct"] = m
    if tg.get("denovo"):
        a = try_(name, "assembly_info", lambda: assembly_stats(base + ".assembly_info.txt"),
                 default=(ERR, ERR))
        rec["assembly_contigs"], rec["assembly_circular"] = a
    if tg.get("cassettes"):
        rec["cassette_insertions"] = try_(name, "cassette_insertions", lambda: cassette_sites(
            os.path.join(d, "annotated", name + ".cassette_insertions.annotated.tsv")))
    return rec


def cassette_sites(path):
    """'cassette@contig:start-end(genes)' per insertion site, ';'-joined; 'none' if no site.
    Columns: contig start end length status cassette orientation ... features."""
    sites = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            genes = sorted({f.split(":", 1)[1].split("(")[0] for f in c[-1].split(";") if ":" in f})
            where = f"{c[0]}:{c[1]}-{c[2]}" if c[4] == "paired" else f"{c[0]}:{c[1]}({c[4]})"
            sites.append(f"{c[5]}@{where}" + (f"[{','.join(genes)}]" if genes else ""))
    return ";".join(sites) if sites else "none"


def run(samples, out_dir, toggles, out_path):
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for s in samples:
            rec = summarize_sample(s, out_dir, toggles)
            w.writerow([rec[c] for c in COLUMNS])
    sys.stderr.write("wrote %s (%d samples)\n" % (out_path, len(samples)))


if "snakemake" in globals():
    _sm = globals()["snakemake"]
    run(list(_sm.params.samples), str(_sm.params.out_dir), dict(_sm.params.toggles),
        str(_sm.output.summary))
