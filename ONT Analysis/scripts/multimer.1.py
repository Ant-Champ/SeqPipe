#!/usr/bin/env python3
"""
multimer.1.py -- detect head-to-tail plasmid multimers (concatemers) from
Oxford Nanopore reads mapped against a *monomeric* reference.

Rationale
---------
A circular monomer of length L cannot yield a read longer than L.  Any read
clearly longer than L is therefore either

  * a genuine multimer (dimer, trimer, ...), i.e. n tandem head-to-tail copies
    of the same sequence, or
  * a library artefact (ligation chimera, inverted duplex, RCA product).
  * a contamination from host dna

A head-to-tail multimer carries no unique junction sequence -- it is the
monomer repeated -- so it is indistinguishable from a monomer at the level of
consensus/variants/coverage.  The only usable signals are the *length* of the
molecule and the *topology* of its alignment blocks:

  true multimer  : all blocks on the same strand, and every junction is
                   continuous around the circle, i.e.
                       (t_start(i+1) - t_end(i))  ==  query_gap   (mod L)
  chimera        : strand flip, or a jump to an unrelated reference position.

Input  : PAF produced by `minimap2 -x map-ont -c --secondary=no ref.fa reads.fq`
Output : per-read TSV (-out) + human readable summary on stdout.
"""

import argparse
import sys
from collections import defaultdict, Counter


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Detect head-to-tail plasmid multimers from ONT reads.")
    p.add_argument("-paf", required=True,
                   help="PAF of reads vs monomeric reference")
    p.add_argument("-ref", required=True,
                   help="Monomeric reference FASTA (used for plasmid length L)")
    p.add_argument("-out", required=True,
                   help="Output per-read classification TSV")
    p.add_argument("-min_frac", type=float, default=0.10,
                   help="Min fraction of full-length reads that must be "
                        "multimeric to call the sample multimeric [0.10]")
    p.add_argument("-min_reads", type=int, default=10,
                   help="Min number of supporting multimer reads [10]")
    p.add_argument("-tol", type=float, default=0.02,
                   help="Junction tolerance, as a fraction of L (min 50 bp) [0.02]")
    p.add_argument("-min_qcov", type=float, default=0.85,
                   help="Min fraction of the read that must align to the "
                        "plasmid; below this the read carries foreign DNA [0.85]")
    p.add_argument("-cut", type=float, default=1.8,
                   help="read_length / L above which a read is concatemeric [1.8]")
    p.add_argument("-min_mapq", type=int, default=0,
                   help="Discard alignment blocks below this MAPQ [0]")
    return p.parse_args()


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def read_fasta_lengths(path):
    """Return {seq_name: length} from a FASTA file."""
    lengths, name, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line)
    if name is not None:
        lengths[name] = n
    return lengths


def read_paf(path, min_mapq):
    """Group PAF records by read name.  Secondary alignments are dropped."""
    per_read = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            tp = "P"
            for tag in f[12:]:
                if tag.startswith("tp:A:"):
                    tp = tag.split(":")[-1]
            if tp == "S":                      # secondary
                continue
            mapq = int(f[11])
            if mapq < min_mapq:
                continue
            per_read[f[0]].append({
                "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]),
                "strand": f[4], "t": f[5], "tlen": int(f[6]),
                "ts": int(f[7]), "te": int(f[8]),
                "nmatch": int(f[9]), "alen": int(f[10]), "mapq": mapq,
            })
    return per_read


def merged_length(intervals):
    """Total length covered by a set of (start, end) intervals."""
    if not intervals:
        return 0
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None:
            cur_s, cur_e = s, e
        elif s > cur_e:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    total += cur_e - cur_s
    return total


def circular_delta(a, b, L):
    """Smallest circular distance between two offsets."""
    return min((a - b) % L, (b - a) % L)


# ----------------------------------------------------------------------
# per-read classification
# ----------------------------------------------------------------------
def classify(blocks, L, tol_bp, min_qcov, cut):
    blocks = sorted(blocks, key=lambda b: b["qs"])
    qlen = blocks[0]["qlen"]

    qcov = merged_length([(b["qs"], b["qe"]) for b in blocks]) / float(qlen)
    ref_span = sum(b["te"] - b["ts"] for b in blocks)
    len_ratio = qlen / float(L)
    span_ratio = ref_span / float(L)
    strands = {b["strand"] for b in blocks}
    same_strand = len(strands) == 1

    # --- junction continuity around the circle ------------------------
    seamless, worst_junction = True, 0
    for a, b in zip(blocks, blocks[1:]):
        gap_q = b["qs"] - a["qe"]                 # may be slightly negative
        if a["strand"] == "+" and b["strand"] == "+":
            delta = (b["ts"] - a["te"]) % L
        elif a["strand"] == "-" and b["strand"] == "-":
            delta = (a["ts"] - b["te"]) % L
        else:
            seamless = False
            worst_junction = L
            break
        d = circular_delta(delta, gap_q, L)
        worst_junction = max(worst_junction, d)
        if d > tol_bp:
            seamless = False

    # --- label --------------------------------------------------------
    if qcov < min_qcov:
        label, copies = "foreign", 0            # read is not purely plasmid
    elif len_ratio < 0.9:
        label, copies = "partial", 0            # fragment, uninformative
    elif len_ratio < cut:
        label, copies = "monomer", 1
    elif same_strand and seamless:
        label, copies = "multimer", int(round(len_ratio))
    else:
        label, copies = "chimera", 0

    return {
        "qlen": qlen, "len_ratio": len_ratio, "span_ratio": span_ratio,
        "n_blocks": len(blocks), "strand": "/".join(sorted(strands)),
        "q_cov": qcov, "seamless": seamless, "max_junction_bp": worst_junction,
        "label": label, "copies": copies,
    }


# ----------------------------------------------------------------------
def histogram(ratios, width=50, bin_size=0.1, max_ratio=6.0):
    """Text histogram of read_length / L."""
    counts = Counter(min(int(r / bin_size), int(max_ratio / bin_size))
                     for r in ratios)
    if not counts:
        return []
    top = max(counts.values())
    lines = []
    for i in range(0, int(max_ratio / bin_size) + 1):
        c = counts.get(i, 0)
        if c == 0:
            continue
        bar = "#" * max(1, int(width * c / top))
        lines.append("  %4.1f - %4.1f x L | %-*s %d"
                     % (i * bin_size, (i + 1) * bin_size, width, bar, c))
    return lines


# ----------------------------------------------------------------------
def main():
    args = parse_args()

    ref_lengths = read_fasta_lengths(args.ref)
    if not ref_lengths:
        sys.exit("ERROR: no sequence found in %s" % args.ref)

    per_read = read_paf(args.paf, args.min_mapq)

    # dominant reference contig (plasmids: normally a single sequence)
    aligned_per_target = Counter()
    for blocks in per_read.values():
        for b in blocks:
            aligned_per_target[b["t"]] += b["alen"]
    if aligned_per_target:
        target = aligned_per_target.most_common(1)[0][0]
    else:
        target = max(ref_lengths, key=ref_lengths.get)

    if target not in ref_lengths:
        sys.exit("ERROR: reference %s absent from %s" % (target, args.ref))
    L = ref_lengths[target]
    tol_bp = max(50.0, args.tol * L)

    results = {}
    for qname, blocks in per_read.items():
        blocks = [b for b in blocks if b["t"] == target]
        if blocks:
            results[qname] = classify(blocks, L, tol_bp, args.min_qcov, args.cut)

    # ---------------- per-read TSV ------------------------------------
    with open(args.out, "w") as fh:
        fh.write("read\tread_len\tlen_ratio\tref_span_ratio\tn_blocks\t"
                 "strand\tquery_cov\tseamless\tmax_junction_bp\t"
                 "class\tcopies\n")
        for qname in sorted(results, key=lambda k: -results[k]["len_ratio"]):
            r = results[qname]
            fh.write("%s\t%d\t%.2f\t%.2f\t%d\t%s\t%.3f\t%s\t%d\t%s\t%d\n"
                     % (qname, r["qlen"], r["len_ratio"], r["span_ratio"],
                        r["n_blocks"], r["strand"], r["q_cov"],
                        "yes" if r["seamless"] else "no",
                        int(r["max_junction_bp"]), r["label"], r["copies"]))

    # ---------------- summary ----------------------------------------
    labels = Counter(r["label"] for r in results.values())
    n_mono = labels["monomer"]
    n_multi = labels["multimer"]
    n_chim = labels["chimera"]
    n_full = n_mono + n_multi                      # informative denominator
    frac = (n_multi / float(n_full)) if n_full else 0.0

    bases_mono = sum(r["qlen"] for r in results.values() if r["label"] == "monomer")
    bases_multi = sum(r["qlen"] for r in results.values() if r["label"] == "multimer")
    frac_bases = (bases_multi / float(bases_mono + bases_multi)
                  if (bases_mono + bases_multi) else 0.0)

    copies = Counter(r["copies"] for r in results.values() if r["label"] == "multimer")

    print("=" * 68)
    print("PLASMID MULTIMER DETECTION")
    print("=" * 68)
    print("reference           : %s" % target)
    print("monomer length (L)  : %d bp" % L)
    print("junction tolerance  : %d bp" % tol_bp)
    print("reads mapped        : %d" % len(results))
    print("")
    print("read classification")
    print("  partial  (<0.9 L)                 : %d" % labels["partial"])
    print("  monomer  (0.9 - %.1f L)            : %d" % (args.cut, n_mono))
    print("  multimer (>=%.1f L, seamless h2t)   : %d" % (args.cut, n_multi))
    print("  chimera  (>=%.1f L, inconsistent)   : %d" % (args.cut, n_chim))
    print("  foreign  (<%.0f%% of read on plasmid) : %d"
          % (args.min_qcov * 100, labels["foreign"]))
    print("")
    if copies:
        print("multimer order (copies per molecule)")
        for n in sorted(copies):
            print("  %d-mer : %d reads" % (n, copies[n]))
        print("")
    print("multimer fraction of full-length reads : %.1f%% (%d / %d)"
          % (frac * 100, n_multi, n_full))
    print("multimer fraction of full-length bases : %.1f%%" % (frac_bases * 100))
    if n_full:
        print("chimera / multimer ratio               : %.2f"
              % (n_chim / float(max(n_multi, 1))))
    print("")
    print("read length distribution (units of L)")
    for line in histogram([r["len_ratio"] for r in results.values()]):
        print(line)
    print("")

    # ---------------- verdict ----------------------------------------
    if n_full == 0:
        verdict = "UNDETERMINED : no full-length read mapped"
    elif n_multi < args.min_reads or frac < args.min_frac:
        verdict = ("MONOMER : no significant multimer signal "
                   "(%.1f%% < %.1f%% or %d < %d supporting reads)"
                   % (frac * 100, args.min_frac * 100, n_multi, args.min_reads))
    else:
        dominant = copies.most_common(1)[0][0] if copies else 2
        state = "PREDOMINANTLY MULTIMERIC" if frac >= 0.80 else "MIXED monomer / multimer"
        verdict = ("MULTIMER : %s, dominant form = %d-mer (%d x %d bp = %d bp), "
                   "%.1f%% of full-length reads"
                   % (state, dominant, dominant, L, dominant * L, frac * 100))

    print("VERDICT: %s" % verdict)

    if n_chim > n_multi and n_chim > 5:
        print("WARNING: chimeric long reads outnumber seamless concatemers; "
              "suspect ligation artefacts rather than true multimers.")
    print("=" * 68)


if __name__ == "__main__":
    main()
