#!/usr/bin/env python3
"""annotate_regions -- overlap a region table with GenBank features.

Ported from alignqc.py (annotate_table). Appends a 'features' column to any TSV
whose first 3 columns are contig/start/end (0-based half-open, '#' header), e.g.
the low_coverage / indels / breakpoints tables. Matching is contig-aware: the
contig must equal the GenBank record.id. A region overlaps a feature if
fs < end and fe > start. Each hit is "type:name(product)" (name = gene, else
locus_tag, else label, else '?'); hits are joined by ';', no hit gives '.'.
Header-only input gives header-only output. Feature types in --skip-types
are ignored.

Output columns: the input columns + features

Python >= 3.8; needs biopython.
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path


def load_features(gb, skip):
    from Bio import SeqIO
    feats = defaultdict(list)                       # contig -> [(start, end, label)]
    for rec in SeqIO.parse(gb, "genbank"):
        for f in rec.features:
            if f.type in skip:
                continue
            q = f.qualifiers
            nm = (q.get("gene") or q.get("locus_tag") or q.get("label") or ["?"])[0]
            prod = (q.get("product") or q.get("note") or [""])[0]
            label = f"{f.type}:{nm}" + (f"({prod})" if prod else "")
            # ';' separates features in the output column, tabs / newlines would break the TSV
            label = " ".join(label.replace(";", ",").split())
            for part in f.location.parts:       # join() / across the origin: each part on its own
                feats[rec.id].append((int(part.start), int(part.end), label))
    return feats


def annotate_table(gb, tsv_in, tsv_out, skip_types="source,gene"):
    skip = set(skip_types.split(","))
    feats = load_features(gb, skip)
    Path(tsv_out).parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_in) as fin, open(tsv_out, "w") as fout:
        for line in fin:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                fout.write(line + "\tfeatures\n")
                continue
            c = line.split("\t")
            contig, s, e = c[0], int(c[1]), int(c[2])
            hits = list(dict.fromkeys(lab for fs, fe, lab in feats.get(contig, []) if fs < e and fe > s))
            fout.write(line + "\t" + (";".join(hits) if hits else ".") + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genbank", required=True, help="GenBank file (record IDs = contig names)")
    ap.add_argument("--in", dest="tsv_in", required=True, help="input TSV (contig/start/end first)")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--skip-types", default="source,gene",
                    help="comma-separated GenBank feature types to ignore [source,gene]")
    args = ap.parse_args()
    try:
        annotate_table(args.genbank, args.tsv_in, args.out, args.skip_types)
    except ImportError as e:
        sys.exit(f"annotate_regions: missing dependency: {e}")
    except (OSError, ValueError, IndexError) as e:
        sys.exit(f"annotate_regions: error: {e}")


if __name__ == "__main__":
    main()
