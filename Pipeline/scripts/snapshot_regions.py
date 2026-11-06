#!/usr/bin/env python3
"""snapshot_regions -- genome-browser-style PNG of every region of a table, drawn
with matplotlib (no Java, no X server: works on any compute node). The
pipeline's tool for read-level pictures: reuse it for new ones.

For each region (--regions: '#' header, contig start end length status cassette
..., 0-based half-open; e.g. the cassette hotspots table) +/- --flank bp:
  * genes / features of the GenBank reference (arrows, names)
  * depth of the sample reads vs its genome (--all-bam)
  * the reads of the genome + cassettes alignment (--cassette-bam) over the
    region, stacked: grey = ordinary read, light = MAPQ < --min-mapq, red ends =
    soft clips, green = read whose mate is on a cassette, purple = read split
    to a cassette (SA tag); at most --max-rows rows
  * dashed lines at the junction positions of the table's 'junctions' column
    (pos:reads), if there is one.
The view is at most --max-view bp: a longer region is shown around its most
supported junction (else its middle), and the picture says so. At most
--max-pictures regions are drawn (the ones with most reads first).
One PNG per region in --outdir: NN_<status>_<cassette>_<contig>_<start>-<end>.png
(1-based coordinates in the name). No region -> NO_CASSETTE_INSERTIONS.txt.

Debugging: pictures_debug.tsv in --outdir always has one row per region (window
drawn, contig found in each BAM, reads fetched / drawn by kind, rows, depth
min / median / max, features shown, seconds, png or the error). A region that
fails is logged with its traceback and the others are still drawn. --debug
also writes <png>.reads.tsv: every read fetched (name, flag, MAPQ, position,
CIGAR, mate, SA, kind, drawn or not).

Needs pysam, matplotlib, biopython.
"""
import argparse
import re
import sys
import time
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.patches import FancyArrow, Rectangle  # noqa: E402

COL = {"read": "#9aa4ae", "lowq": "#d6dbe0", "clip": "#d64034", "mate": "#2e9e5b",
       "split": "#8a4fc2", "depth": "#4a7fb5", "gene": "#5b6670", "junction": "#e08a1e"}


def read_regions(path):
    regions, cols = [], None
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                cols = [c.lstrip("#") for c in f]
                continue
            if len(f) < 6:
                raise SystemExit(f"regions table needs >= 6 columns: {line.strip()}")
            junctions = []
            if cols and "junctions" in cols and cols.index("junctions") < len(f):
                for x in f[cols.index("junctions")].split(","):
                    if ":" in x:
                        p, n = x.split(":")
                        junctions.append((int(p), int(n)))
            regions.append(dict(contig=f[0], start=int(f[1]), end=int(f[2]), status=f[4],
                                cassette=f[5], junctions=junctions, row=dict(zip(cols or [], f))))
    return regions


def features(genbank):
    from Bio import SeqIO
    out = {}
    for rec in SeqIO.parse(genbank, "genbank"):
        fs = []
        for f in rec.features:
            if f.type in ("source", "gene"):
                continue
            q = f.qualifiers
            name = (q.get("gene") or q.get("locus_tag") or q.get("label") or [f.type])[0]
            for part in f.location.parts:       # join() / across the origin: each part on its own
                fs.append((int(part.start), int(part.end), part.strand or f.location.strand or 1, name))
        out[rec.id] = fs
    return out


def stack(reads, max_rows):
    """Greedy row packing: [(row, read)]; reads beyond max_rows are dropped."""
    ends, placed = [], []
    for r in sorted(reads, key=lambda r: r["start"]):
        for i, e in enumerate(ends):
            if r["start"] > e + 2:
                ends[i] = r["end"]
                placed.append((i, r))
                break
        else:
            if len(ends) < max_rows:
                ends.append(r["end"])
                placed.append((len(ends) - 1, r))
    return placed, len(ends)


def draw(reg, a, b, args, feats, all_bam, cas_bam, cass_set, out_png, sample, reads_tsv=None):
    """Draw one region; returns the debug statistics of what was found / drawn."""
    contig = reg["contig"]
    st = {"in_all_bam": contig in all_bam.references, "in_cassette_bam": contig in cas_bam.references}
    # depth (all reads vs genome), like samtools depth
    depth = [0] * (b - a)
    if contig in all_bam.references:
        cov = all_bam.count_coverage(contig, a, b, quality_threshold=0,
                                     read_callback=lambda r: not (r.flag & 0x704))
        depth = [sum(x) for x in zip(*cov)]
    # reads (genome + cassettes alignment)
    reads = []
    if contig in cas_bam.references:
        for r in cas_bam.fetch(contig, a, b):
            if r.is_unmapped or r.flag & 0x900 or not r.cigartuples:
                continue
            ct = r.cigartuples
            kind = "read"
            if r.has_tag("SA") and any(x.split(",")[0] in cass_set for x in r.get_tag("SA").split(";") if x):
                kind = "split"
            elif r.is_paired and not r.mate_is_unmapped and r.next_reference_name in cass_set:
                kind = "mate"
            elif r.mapping_quality < args.min_mapq:
                kind = "lowq"
            reads.append(dict(start=r.reference_start, end=r.reference_end, kind=kind,
                              lclip=ct[0][1] if ct[0][0] == 4 else 0,
                              rclip=ct[-1][1] if ct[-1][0] == 4 else 0,
                              rev=r.is_reverse,
                              dbg=(r.query_name, r.flag, r.mapping_quality, r.cigarstring,
                                   r.next_reference_name if r.is_paired else ".",
                                   r.get_tag("SA") if r.has_tag("SA") else ".")
                              if reads_tsv else None))
    placed, nrows = stack(reads, args.max_rows)
    drawn = {id(r) for _, r in placed}
    sd = sorted(depth)
    st.update(reads_fetched=len(reads), reads_drawn=len(placed), rows=nrows,
              depth_min=sd[0] if sd else 0, depth_median=sd[len(sd) // 2] if sd else 0,
              depth_max=sd[-1] if sd else 0,
              features=sum(1 for s_, e_, _, _ in feats.get(contig, []) if e_ >= a and s_ <= b),
              **{f"n_{k}": sum(1 for r in reads if r["kind"] == k) for k in ("read", "lowq", "mate", "split")},
              n_clipped=sum(1 for r in reads if r["lclip"] or r["rclip"]))
    if reads_tsv:
        with open(reads_tsv, "w") as fh:
            fh.write("name\tflag\tmapq\tstart\tend\tcigar\tmate_contig\tSA\tkind\tdrawn\n")
            for r in reads:
                n, fl, mq, cg, mc, sa = r["dbg"]
                fh.write(f"{n}\t{fl}\t{mq}\t{r['start']}\t{r['end']}\t{cg}\t{mc}\t{sa}\t{r['kind']}\t"
                         f"{'yes' if id(r) in drawn else 'no'}\n")

    h_reads = max(1.2, 0.07 * nrows)
    fig = plt.figure(figsize=(14, 2.4 + h_reads))
    gs = fig.add_gridspec(3, 1, height_ratios=[0.6, 1.2, h_reads], hspace=0.08, top=0.93)
    ax_g, ax_d, ax_r = (fig.add_subplot(gs[i]) for i in range(3))
    for ax in (ax_g, ax_d, ax_r):
        ax.set_xlim(a, b)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    # features
    ax_g.set_ylim(-1.2, 1.4)
    ax_g.set_yticks([])
    ax_g.set_xticks([])
    ax_g.spines["left"].set_visible(False)
    ax_g.spines["bottom"].set_visible(False)
    for s, e, strand, name in feats.get(contig, []):
        if e < a or s > b:
            continue
        s2, e2 = max(s, a), min(e, b)
        y = 0.3 if strand >= 0 else -0.3
        hl = min((b - a) * 0.01, e2 - s2)
        x0, dx = (s2, e2 - s2) if strand >= 0 else (e2, -(e2 - s2))
        ax_g.add_patch(FancyArrow(x0, y, dx, 0, width=0.35, head_width=0.55, head_length=hl,
                                  length_includes_head=True, color=COL["gene"], lw=0))
        ax_g.text((s2 + e2) / 2, y + (0.45 if strand >= 0 else -0.75), name, ha="center",
                  fontsize=8, clip_on=True)

    # depth
    ax_d.fill_between(range(a, b), depth, step="post", color=COL["depth"], lw=0)
    ax_d.set_ylabel("depth", fontsize=8)
    ax_d.set_xticks([])
    ax_d.tick_params(labelsize=7)

    # reads
    for row, r in placed:
        y = nrows - row
        c = COL[r["kind"]]
        ax_r.add_patch(Rectangle((r["start"], y - 0.4), r["end"] - r["start"], 0.8, color=c, lw=0))
        if r["lclip"]:
            ax_r.add_patch(Rectangle((r["start"] - r["lclip"], y - 0.4), r["lclip"], 0.8,
                                     color=COL["clip"], lw=0))
        if r["rclip"]:
            ax_r.add_patch(Rectangle((r["end"], y - 0.4), r["rclip"], 0.8, color=COL["clip"], lw=0))
    ax_r.set_ylim(0, nrows + 1)
    ax_r.set_yticks([])
    ax_r.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax_r.tick_params(labelsize=8)
    ax_r.set_xlabel(f"{contig} (bp)", fontsize=8)

    top = ax_d.get_ylim()[1]
    for k, (p, n) in enumerate(reg["junctions"]):
        for ax in (ax_d, ax_r):
            ax.axvline(p, color=COL["junction"], ls="--", lw=1)
        ax_d.text(p, top * (1 - 0.18 * k), f" junction {p:,} ({n} reads)", color=COL["junction"],
                  fontsize=7, va="top")

    counts = "  ".join(f"{k}={reg['row'][k]}" for k in ("n_reads", "n_mates", "n_genome_clip",
                                                       "n_cassette_clip", "n_split") if k in reg["row"])
    view = "" if (a <= reg["start"] and b >= reg["end"]) else \
        f"   (view {a + 1:,}-{b:,} of a {reg['end'] - reg['start']:,} bp hotspot)"
    fig.suptitle(f"{sample}  {reg['cassette']}  {reg['status']}  {contig}:{reg['start'] + 1:,}-{reg['end']:,}"
                 f"   {counts}{view}", fontsize=10, x=0.01, y=0.99, ha="left")
    handles = [Rectangle((0, 0), 1, 1, color=COL[k]) for k in ("read", "lowq", "clip", "mate", "split")]
    fig.legend(handles, ["read", f"MAPQ < {args.min_mapq}", "soft clip", "mate on a cassette",
                         "split to a cassette"], loc="lower right", fontsize=7, ncol=5, frameon=False,
               bbox_to_anchor=(0.99, -0.02))
    if len(placed) < len(reads):
        ax_r.text(a, 0, f" {len(reads) - len(placed)} more reads not shown", fontsize=7, va="bottom")
    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    return st


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--all-bam", required=True, help="sample reads vs its genome")
    ap.add_argument("--cassette-bam", required=True, help="sample reads vs genome + cassettes")
    ap.add_argument("--genbank", required=True)
    ap.add_argument("--cassettes", nargs="+", required=True, help="cassette contig names")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample", default="")
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--max-rows", type=int, default=150)
    ap.add_argument("--max-view", type=int, default=10000, help="widest view drawn (bp) [10000]")
    ap.add_argument("--max-pictures", type=int, default=50, help="most regions drawn [50]")
    ap.add_argument("--dpi", type=int, default=110)
    ap.add_argument("--debug", action="store_true", help="also write <png>.reads.tsv (every read fetched)")
    args = ap.parse_args()

    import pysam
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    regions = read_regions(args.regions)
    if not regions:
        (out / "NO_CASSETTE_INSERTIONS.txt").write_text("No cassette hotspot / insertion site was found.\n")
        return
    import matplotlib as _mpl
    print(f"snapshot_regions: {len(regions)} region(s) in {args.regions}; pysam {pysam.__version__}, "
          f"matplotlib {_mpl.__version__}, python {sys.version.split()[0]}", flush=True)
    feats = features(args.genbank)
    all_bam = pysam.AlignmentFile(args.all_bam)
    cas_bam = pysam.AlignmentFile(args.cassette_bam)
    lens = dict(zip(all_bam.references, all_bam.lengths))
    bad = [c for c in set(all_bam.references) | set(cas_bam.references) if c in (".", "*")]
    if bad:
        sys.exit(f"snapshot_regions: contig(s) named {bad} in the BAMs: htslib reads them as 'the whole "
                 f"file', so every region would load all reads. Rerun the pipeline from the start "
                 f"(the staged reference now renames them).")
    print(f"all-bam contigs: {', '.join(all_bam.references)}", flush=True)
    print(f"cassette-bam contigs: {', '.join(cas_bam.references)}", flush=True)
    print(f"genbank records with features: {', '.join(f'{k} ({len(v)})' for k, v in feats.items())}", flush=True)
    if len(regions) > args.max_pictures:
        keep = sorted(regions, key=lambda r: -int(r["row"].get("n_reads", 0) or 0))[:args.max_pictures]
        print(f"WARNING: {len(regions)} regions, drawing the {args.max_pictures} with most reads "
              f"(--max-pictures)")
        regions = [r for r in regions if r in keep]
    debug_rows = []
    for i, reg in enumerate(regions, 1):
        a = max(0, reg["start"] - args.flank)
        b = min(lens.get(reg["contig"], reg["end"] + args.flank), reg["end"] + args.flank)
        if b - a > args.max_view:                   # long hotspot: around its best junction
            c = reg["junctions"][0][0] if reg["junctions"] else (reg["start"] + reg["end"]) // 2
            a = max(0, c - args.max_view // 2)
            b = min(lens.get(reg["contig"], a + args.max_view), a + args.max_view)
        label = re.sub(r"[^A-Za-z0-9._-]+", "_", f"{reg['status']}_{reg['cassette']}").strip("_")
        tag = re.sub(r"[^A-Za-z0-9._-]+", "_", reg["contig"])
        png = out / f"{i:02d}_{label}_{tag}_{reg['start'] + 1}-{reg['end']}.png"
        print(f"[{i}/{len(regions)}] {reg['contig']}:{reg['start']}-{reg['end']} ({reg['status']} "
              f"{reg['cassette']}): drawing {reg['contig']}:{a}-{b} ({b - a:,} bp)", flush=True)
        t0 = time.time()
        row = dict(n=i, contig=reg["contig"], start=reg["start"], end=reg["end"], status=reg["status"],
                   cassette=reg["cassette"], view_start=a, view_end=b, view_bp=b - a,
                   junctions=",".join(f"{p}:{n}" for p, n in reg["junctions"]) or ".")
        try:
            row.update(draw(reg, a, b, args, feats, all_bam, cas_bam, set(args.cassettes), png, args.sample,
                            reads_tsv=str(png) + ".reads.tsv" if args.debug else None))
            row["png"], row["error"] = png.name, "."
            print(f"    {row['reads_fetched']} reads fetched, {row['reads_drawn']} drawn in {row['rows']} rows, "
                  f"depth {row['depth_min']}-{row['depth_max']} (median {row['depth_median']}); "
                  f"{time.time() - t0:.1f} s -> {png}", flush=True)
            if not row["in_all_bam"] or not row["in_cassette_bam"]:
                print(f"    WARNING: contig {reg['contig']} missing from "
                      f"{'the all-reads BAM' if not row['in_all_bam'] else 'the genome+cassettes BAM'}",
                      flush=True)
        except Exception as err:                    # one bad region must not stop the others
            row["png"], row["error"] = ".", f"{type(err).__name__}: {err}"
            print(f"    ERROR: {row['error']}", flush=True)
            traceback.print_exc(file=sys.stdout)
            plt.close("all")
        row["seconds"] = f"{time.time() - t0:.1f}"
        debug_rows.append(row)
    cols = ["n", "contig", "start", "end", "status", "cassette", "view_start", "view_end", "view_bp",
            "junctions", "in_all_bam", "in_cassette_bam", "reads_fetched", "reads_drawn", "rows", "n_read",
            "n_lowq", "n_mate", "n_split", "n_clipped", "depth_min", "depth_median", "depth_max",
            "features", "seconds", "png", "error"]
    with open(out / "pictures_debug.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in debug_rows:
            fh.write("\t".join(str(r.get(c, ".")) for c in cols) + "\n")
    print(f"-> {out / 'pictures_debug.tsv'}")


if __name__ == "__main__":
    main()
