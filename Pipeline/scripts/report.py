#!/usr/bin/env python3
"""
report.py -- one self-contained HTML report of the run.  Used through
Snakemake's `script:` directive (rule html_report, rules/report.smk).

snakemake.params : samples (list of dicts), out_dir, toggles (dict), embed_igv
snakemake.output : html
snakemake.config : shown in the report ("Run settings")

Python gathers every small result table into one JSON object (tables are
kept as text: columns + rows), the binned depth (depth_profile), the genes
of each GenBank reference and the IGV snapshots (base64 PNG). The page
itself (HTML/CSS/JS below) builds both views from that object:
  * By sample   : one page per sample with every analysis
  * By analysis : one page per analysis with every sample
Nothing is loaded from the internet: the file opens offline.
"""

import base64
import csv
import datetime
import glob
import html
import json
import os
import sys

MAX_ROWS = 5000          # rows kept per table (the TSV on disk is complete)

# (key, title, toggle or None = always, description shown at the top of the section)
SECTIONS = [
    ("overview",  "Overview",  None,
     "One row per sample (summary.tsv). Click a sample name to open its page."),
    ("qc",        "Read QC",   None,
     "seqkit stats of the raw and filtered reads; fastp report for short reads."),
    ("coverage",  "Alignment & coverage", None,
     "Depth along every replicon (drag on a plot to zoom, double-click to reset; "
     "red = low-coverage region, teal = high-coverage region, green = cassette insertion site, light green = cassette hotspot, orange = breakpoint, purple = indel, "
     "black ticks = variants). samtools coverage per replicon and flagstat."),
    ("variants",  "Variants",  "variants",
     "Small variants against the reference (medaka for ONT, bcftools for short reads). "
     "Click a row to see it on the coverage plot."),
    ("lowcov",    "Low / high coverage", "coverage",
     "Regions below the low-coverage threshold (deletions) and above the high-coverage threshold "
     "of their replicon (duplications, amplifications; fold ~2 = one extra copy), with the genes they hit."),
    ("sv",        "Large indels & breakpoints", "sv",
     "Large indels from CIGARs and clipped-read breakpoint clusters."),
    ("cassettes", "Cassette insertions", "cassettes",
     "Where each cassette of cassette_dir is inserted (junction reads), the depth along each cassette "
     "(all reads, including those of the ends cassettes share) and the IGV snapshots."),
    ("plasmids",  "Plasmid panel", "plasmid_panel",
     "Breadth (% covered) and mean depth of every panel plasmid in every sample."),
    ("multimer",  "Plasmid multimers", "multimer",
     "Head-to-tail multimers from full-length ONT reads (ONT samples only)."),
    ("assembly",  "De novo assembly", "denovo",
     "Contigs of the Flye (ONT) or Unicycler (short reads) assembly."),
]

READ_STATS_COLS = ["file", "num_seqs", "sum_len", "min_len", "avg_len", "max_len",
                   "N50", "Q20(%)", "Q30(%)", "AvgQual", "GC(%)"]


def log(msg):
    sys.stderr.write(msg + "\n")


# ----------------------------------------------------------------------
# readers
# ----------------------------------------------------------------------
def read_table(path, keep=None):
    """TSV -> {"columns", "rows", "total", "path"}; header = first line (a
    leading '#' is dropped). None if the file is missing."""
    if not path or not os.path.isfile(path):
        return None
    with open(path, newline="") as fh:
        lines = [l.rstrip("\n").rstrip("\r") for l in fh if l.strip()]
    if not lines:
        return {"columns": [], "rows": [], "total": 0}
    if lines[0].lstrip().startswith("# SKIPPED"):
        return {"columns": [], "rows": [], "total": 0, "note": " ".join(l.lstrip("# ") for l in lines)}
    cols = lines[0].lstrip("#").split("\t")
    rows = [l.split("\t") for l in lines[1:] if not l.startswith("#")]
    if keep:
        idx = [cols.index(c) for c in keep if c in cols]
        cols = [cols[i] for i in idx]
        rows = [[r[i] if i < len(r) else "" for i in idx] for r in rows]
    t = {"columns": cols, "rows": rows[:MAX_ROWS], "total": len(rows)}
    return t


def read_text(path):
    if not path or not os.path.isfile(path):
        return None
    with open(path) as fh:
        return fh.read()


def read_depth(path):
    """depth_bins.tsv -> {contig: {"len", "bin", "v": [mean depth per bin]}} (file order)."""
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c, s, e, d = line.rstrip("\n").split("\t")
            rec = out.setdefault(c, {"len": 0, "bin": int(e) - int(s) + 1, "v": []})
            rec["len"] = int(e)
            rec["v"].append(round(float(d), 1))
    return out


def read_genes(gb_path, types=("CDS", "rRNA", "tRNA", "ncRNA")):
    """GenBank -> {contig id: [[start, end, strand, name, product], ...]} (features of
    these types; types=None: every feature but 'source')."""
    from Bio import SeqIO
    genes = {}
    for rec in SeqIO.parse(gb_path, "genbank"):
        feats = []
        for f in rec.features:
            if (f.type not in types) if types else f.type == "source":
                continue
            q = f.qualifiers
            name = (q.get("gene") or q.get("locus_tag") or q.get("label") or [f.type])[0]
            prod = (q.get("product") or q.get("note") or [""])[0]
            feats.append([int(f.location.start) + 1, int(f.location.end),
                          1 if (f.location.strand or 1) > 0 else -1, name, prod[:80]])
        genes[rec.id] = feats
    return genes


def png_data_uri(path):
    with open(path, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode("ascii")


def col(t, *names):
    """Index of the first matching column name in table t, else None."""
    if not t:
        return None
    for n in names:
        if n in t["columns"]:
            return t["columns"].index(n)
    return None


def events_from(t, kind, label_cols):
    """Plot overlays from a table that has contig/start(/end) columns."""
    ci, si, ei = col(t, "contig", "CHROM"), col(t, "start", "POS"), col(t, "end")
    if ci is None or si is None:
        return []
    li = [col(t, c) for c in label_cols]
    ev = []
    for r in t["rows"]:
        try:
            s = int(r[si])
            e = int(r[ei]) if ei is not None else s
        except (ValueError, IndexError):
            continue
        label = " ".join(r[i] for i in li if i is not None and i < len(r))
        ev.append([kind, r[ci], s, max(s, e), label])
    return ev


# ----------------------------------------------------------------------
# per sample
# ----------------------------------------------------------------------
def gather_sample(s, out_dir, tg, embed_igv):
    name, ont = s["sample_name"], s["platform"] == "ont"
    d = os.path.join(out_dir, name)
    base = os.path.join(d, name)
    ann = os.path.join(d, "annotated", name)
    blocks = {}            # section key -> list of blocks
    events = []

    def add(sec, title, table=None, text=None, images=None, links=None, note=None):
        blocks.setdefault(sec, []).append({k: v for k, v in dict(
            title=title, table=table, text=text, images=images, links=links, note=note).items()
            if v is not None})

    def rel(p):
        return os.path.relpath(p, out_dir).replace(os.sep, "/")

    # QC
    add("qc", "Read statistics (seqkit)", table=read_table(
        os.path.join(d, "qc", name + ".read_stats.tsv"), keep=READ_STATS_COLS))
    fastp = os.path.join(d, "qc", name + ".fastp.html")
    if os.path.isfile(fastp):
        add("qc", "fastp", links=[[rel(fastp), "fastp report (filtering, quality, adapters)"]])

    # depth check: a low-diversity sample has nothing more to show
    dc = read_table(os.path.join(d, "qc", name + ".depth_check.tsv"))
    low_div = None
    if dc and dc["rows"] and dc["rows"][0][-1] != "ok":
        r = dict(zip(dc["columns"], dc["rows"][0]))
        low_div = (f"Low diversity: {int(r['filt_bases']):,} filtered bases for a {int(r['ref_length']):,} bp "
                   f"reference = {r['theoretical_depth']}x theoretical depth, below min_theoretical_depth "
                   f"({r['min_theoretical_depth']}x). The sample was not aligned or analysed.")
        add("qc", "Depth check", table=dc)
        return {
            "platform": s["platform"], "reference": s["reference"],
            "reads": [x for x in (s.get("reads1"), s.get("reads2")) if x],
            "blocks": blocks, "events": [], "depth": {}, "low_diversity": low_div,
        }

    # alignment & coverage (the plot is drawn from depth / genes / events)
    add("coverage", "Coverage per replicon (samtools coverage)", table=read_table(base + ".ref_coverage.tsv"))
    add("coverage", "Alignment flags (samtools flagstat)", text=read_text(base + ".flagstat.txt"))

    if tg.get("variants"):
        t = read_table(base + ".annotated.gatc.simple.tsv") if tg.get("gatc") else None
        t = t or read_table(base + ".simple.tsv")
        add("variants", "Variants", table=t)
        events += events_from(t, "variant", ["REF", "ALT"])
    if tg.get("coverage"):
        t = read_table(ann + ".low_coverage.annotated.tsv")
        add("lowcov", "Low-coverage regions", table=t)
        events += events_from(t, "lowcov", ["mean_depth"])
        t = read_table(ann + ".high_coverage.annotated.tsv")
        add("lowcov", "High-coverage regions", table=t)
        events += events_from(t, "highcov", ["fold"])
    if tg.get("sv"):
        t = read_table(ann + ".indels.annotated.tsv")
        add("sv", "Large indels", table=t)
        events += events_from(t, "indel", ["type", "median_len"])
        t = read_table(ann + ".breakpoints.annotated.tsv")
        add("sv", "Breakpoints (clipped-read clusters)", table=t)
        events += events_from(t, "breakpoint", ["side", "n_reads"])
    if tg.get("cassettes"):
        t = read_table(ann + ".cassette_insertions.annotated.tsv")
        add("cassettes", "Insertion sites", table=t)
        events += events_from(t, "cassette", ["cassette", "orientation"])
        t = read_table(ann + ".cassette_hotspots.annotated.tsv")
        add("cassettes", "Hotspots (everywhere reads on a cassette lead back to the genome)", table=t)
        events += events_from(t, "hotspot", ["cassette", "n_reads"])
        add("cassettes", "Junctions", table=read_table(ann + ".cassette_junctions.annotated.tsv"))
        add("cassettes", "Reads per cassette", table=read_table(
            os.path.join(d, "cassettes", name + ".cassette_counts.tsv")))
        igv_dir = os.path.join(d, "cassettes", "igv_snapshots")
        pngs = sorted(glob.glob(os.path.join(igv_dir, "*.png")))
        if pngs:
            imgs = [[png_data_uri(p) if embed_igv else rel(p), os.path.basename(p)[:-4]] for p in pngs]
            add("cassettes", "IGV snapshots", images=imgs)
        elif os.path.isfile(os.path.join(igv_dir, "IGV_FAILED.txt")):
            add("cassettes", "IGV snapshots", note="IGV failed for this sample (see "
                + rel(os.path.join(igv_dir, "IGV_FAILED.txt")) + ").")
        elif os.path.isfile(os.path.join(igv_dir, "NO_CASSETTE_INSERTIONS.txt")):
            add("cassettes", "IGV snapshots", note="No insertion site, so no snapshot.")
        elif tg.get("igv"):
            add("cassettes", "IGV snapshots", note="No snapshot found in " + rel(igv_dir))
        else:
            add("cassettes", "IGV snapshots", note="IGV snapshots are off (set run_igv_snapshots: true).")
    if tg.get("plasmid_panel"):
        add("plasmids", "Plasmid panel coverage (samtools coverage)",
            table=read_table(base + ".plasmid_coverage.tsv"))
    if tg.get("multimer"):
        if ont:
            add("multimer", "Multimer summary", table=read_table(base + ".multimer.summary.tsv"))
        else:
            add("multimer", "Multimer summary", note="ONT samples only.")
    if tg.get("denovo"):
        t = read_table(base + ".assembly_info.txt")
        add("assembly", "Assembly contigs", table=t, note=t.get("note") if t else None,
            links=[[rel(base + ".assembly.fasta"), "assembly.fasta"],
                   [rel(base + ".assembly_graph.gfa"), "assembly_graph.gfa (open in Bandage)"]])

    # a table that is None (file missing) becomes a visible note
    for sec in blocks.values():
        for b in sec:
            if "table" in b and b["table"] is None:
                del b["table"]
                b.setdefault("note", "File not found (step failed or skipped).")
    return {
        "platform": s["platform"], "reference": s["reference"],
        "reads": [x for x in (s.get("reads1"), s.get("reads2")) if x],
        "blocks": blocks, "events": events,
        "depth": read_depth(os.path.join(d, "coverage", name + ".depth_bins.tsv")),
        "cdepth": read_depth(os.path.join(d, "cassettes", name + ".cassette_depth.tsv")) if tg.get("cassettes") else {},
    }


def cassette_features(cassettes):
    """{cassette name: [[start, end, strand, name, note], ...]} from the cassette
    GenBank files (FASTA cassettes have none)."""
    from Bio import SeqIO
    out = {}
    for c in cassettes:
        if c.get("format") != "genbank":
            continue
        try:
            rec = list(SeqIO.parse(c["path"], "genbank"))[c["index"]]
            feats = []
            for f in rec.features:
                if f.type == "source":
                    continue
                q = f.qualifiers
                name = (q.get("gene") or q.get("label") or q.get("locus_tag") or q.get("name")
                        or q.get("product") or [f.type])[0]
                note = (q.get("product") or q.get("note") or [f.type])[0]
                feats.append([int(f.location.start) + 1, int(f.location.end),
                              1 if (f.location.strand or 1) > 0 else -1, name, str(note)[:80]])
            out[c["name"]] = feats
        except Exception as err:           # annotations are a nice-to-have
            log("WARNING: features of cassette %s: %s" % (c["name"], err))
    return out


def main(samples, out_dir, toggles, embed_igv, config, out_html, cassettes=()):
    data = {
        "created": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "out_dir": os.path.abspath(out_dir),
        "toggles": toggles,
        "sections": [{"key": k, "title": t, "desc": desc}
                     for k, t, tog, desc in SECTIONS if tog is None or toggles.get(tog)],
        "summary": read_table(os.path.join(out_dir, "summary.tsv")),
        "coverage_summary": read_table(os.path.join(out_dir, "coverage_summary.tsv")),
        "config": {k: v for k, v in sorted(config.items()) if isinstance(v, (str, int, float, bool))},
        "order": [s["sample_name"] for s in samples],
        "samples": {}, "genes": {}, "ref_of": {},
        "cassette_genes": cassette_features(cassettes),
    }
    for s in samples:
        name = s["sample_name"]
        log("sample " + name)
        data["samples"][name] = gather_sample(s, out_dir, toggles, embed_igv)
        ref = s["reference"]
        data["ref_of"][name] = ref
        if ref not in data["genes"]:
            try:
                data["genes"][ref] = read_genes(os.path.join(out_dir, name, name + ".ref.gb"))
            except Exception as err:          # genes are a nice-to-have on the plot
                log("WARNING: genes of %s: %s" % (ref, err))
                data["genes"][ref] = {}

    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    # named after the run folder (the config's folder: snakemake runs there)
    page = PAGE.replace("__TITLE__", html.escape("Run report - " + os.path.basename(os.getcwd()))) \
               .replace("__DATA__", payload)
    with open(out_html, "w", encoding="utf-8") as fh:
        fh.write(page)
    log("wrote %s (%.1f MB)" % (out_html, os.path.getsize(out_html) / 1e6))


# ----------------------------------------------------------------------
# the page
# ----------------------------------------------------------------------
PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root{--bg:#f7f7f5;--panel:#fff;--ink:#1d1d1b;--muted:#6b6b66;--line:#e2e1dc;--accent:#2f5d8a;
 --depth:#4a7fb5;--low:rgba(214,64,52,.18);--lowedge:#d64034;--high:rgba(23,150,150,.2);--highedge:#179696;--cas:#2e9e5b;--bp:#e08a1e;--indel:#8a4fc2;--var:#1d1d1b;
 --hot:#8fcf9f;--gene:#9aa4ae;--good:#2e9e5b;--warn:#d9a21b;--bad:#d64034;--chip:#eef2f6}
@media (prefers-color-scheme: dark){:root{--bg:#161615;--panel:#1f1f1e;--ink:#ecebe6;--muted:#9c9b95;--line:#33332f;
 --accent:#8db6e0;--depth:#6f9fd1;--low:rgba(240,100,90,.22);--var:#ecebe6;--gene:#6c7680;--chip:#26303a}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
header{position:sticky;top:0;z-index:5;background:var(--panel);border-bottom:1px solid var(--line);
 display:flex;align-items:center;gap:16px;padding:10px 16px;flex-wrap:wrap}
header h1{font-size:16px;margin:0;font-weight:600}
header .meta{color:var(--muted);font-size:12px}
.seg{display:inline-flex;border:1px solid var(--line);border-radius:8px;overflow:hidden;margin-left:auto}
.seg button{border:0;background:transparent;color:var(--ink);padding:6px 14px;font:inherit;cursor:pointer}
.seg button.on{background:var(--accent);color:#fff}
.wrap{display:flex;min-height:calc(100vh - 50px)}
nav{width:230px;flex:none;border-right:1px solid var(--line);padding:10px 8px;position:sticky;top:50px;
 height:calc(100vh - 50px);overflow:auto}
nav a{display:block;padding:6px 10px;border-radius:6px;color:var(--ink);text-decoration:none;cursor:pointer}
nav a:hover{background:var(--chip)}
nav a.on{background:var(--accent);color:#fff}
nav a small{display:block;color:var(--muted);font-size:11px}
nav a.on small{color:#e6eef7}
main{flex:1;min-width:0;padding:16px 20px 60px}
h2{font-size:20px;margin:4px 0 4px}
h3{font-size:15px;margin:22px 0 6px;padding-bottom:4px;border-bottom:1px solid var(--line)}
h4{font-size:13px;margin:14px 0 6px;color:var(--muted);font-weight:600}
.desc{color:var(--muted);margin:0 0 10px}
.tiles{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:8px;margin:12px 0}
.tile{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:8px 10px}
.tile b{display:block;font-size:18px;font-weight:600}
.tile span{color:var(--muted);font-size:12px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:10px 12px;margin:8px 0}
.tbl{overflow:auto;max-height:520px;border:1px solid var(--line);border-radius:6px;background:var(--panel)}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{padding:4px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top;white-space:nowrap}
td.longtxt{white-space:normal}
td.longtxt .clip{width:420px;max-height:2.9em;overflow:hidden;cursor:pointer}
td.longtxt .clip.open{max-height:none}
td.longtxt .more{color:var(--accent);font-size:11px;cursor:pointer;white-space:nowrap}
th{position:sticky;top:0;background:var(--panel);cursor:pointer;user-select:none;font-weight:600}
th:after{content:" \2195";color:var(--muted);font-size:10px}
tr.go{cursor:pointer} tr.go:hover td{background:var(--chip)}
.tools{display:flex;gap:8px;align-items:center;margin:4px 0;color:var(--muted);font-size:12px}
.tools input{padding:3px 6px;border:1px solid var(--line);border-radius:5px;background:var(--panel);color:var(--ink)}
.note{color:var(--muted);font-style:italic}
.banner{background:var(--bad);color:#fff;border-radius:8px;padding:10px 14px;margin:10px 0;font-weight:600}
nav a .lowdiv{color:var(--bad);font-weight:600}
nav a.on .lowdiv{color:#fff}
pre{background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:8px;overflow:auto;font-size:12px}
.chip{display:inline-block;padding:1px 7px;border-radius:10px;background:var(--chip);font-size:11px;margin-right:4px}
.plot{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:6px 8px 4px;margin:8px 0}
.plot .bar{display:flex;gap:8px;align-items:center;font-size:12px;color:var(--muted)}
.plot .bar b{color:var(--ink)}
.plot .bar button{border:1px solid var(--line);background:var(--panel);color:var(--ink);border-radius:5px;padding:1px 8px;cursor:pointer}
.plot canvas{width:100%;display:block;cursor:crosshair}
.tip{position:fixed;pointer-events:none;background:var(--panel);border:1px solid var(--line);border-radius:6px;
 padding:6px 8px;font-size:12px;box-shadow:0 4px 14px rgba(0,0,0,.15);display:none;z-index:9;max-width:360px}
.gallery{display:grid;grid-template-columns:1fr;gap:12px}
.gallery figure{margin:0;background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:6px}
.gallery img{width:100%;display:block;cursor:zoom-in}
.gallery figcaption{font-size:12px;color:var(--muted);padding-top:4px}
.ok{color:var(--good)} .warn{color:var(--warn)} .bad{color:var(--bad)}
td.heat{text-align:center}
a{color:var(--accent)}
@media (max-width:760px){.wrap{display:block}nav{position:static;width:auto;height:auto;border-right:0;border-bottom:1px solid var(--line);
 display:flex;flex-wrap:wrap}main{padding:12px 16px}}
</style></head><body>
<header><h1 id="title"></h1><span class="meta" id="meta"></span>
<div class="seg"><button id="m-sample">By sample</button><button id="m-tool">By analysis</button></div></header>
<div class="wrap"><nav id="nav"></nav><main id="main"></main></div>
<div class="tip" id="tip"></div>
<script id="data" type="application/json">__DATA__</script>
<script>
"use strict";
const D = JSON.parse(document.getElementById("data").textContent);
const $ = (t, a = {}, ...kids) => { const e = document.createElement(t);
  for (const [k, v] of Object.entries(a)) { if (k === "class") e.className = v; else if (k === "html") e.innerHTML = v;
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v); else e.setAttribute(k, v); }
  for (const k of kids.flat()) if (k != null) e.append(k.nodeType ? k : document.createTextNode(k)); return e; };
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const fmtBp = x => x >= 1e6 ? (x / 1e6).toFixed(x >= 1e7 ? 1 : 2) + " Mb" : x >= 1e3 ? (x / 1e3).toFixed(x >= 1e4 ? 0 : 1) + " kb" : x + " bp";
const num = v => { const n = parseFloat(String(v).replace(/,/g, "")); return isNaN(n) ? null : n; };
const SEC = Object.fromEntries(D.sections.map(s => [s.key, s]));
let state = { mode: "sample", sample: D.order[0], tool: "overview" };

document.getElementById("title").textContent = document.title;
document.getElementById("meta").textContent = `${D.order.length} sample(s) · ${D.created} · ${D.out_dir}`;

// ---------------------------------------------------------------- tables
function table(t, opts = {}) {
  // t = {columns, rows, total}; opts.sample = sample for row->plot links, opts.sampleCol = index of a sample column
  const box = $("div");
  if (!t.rows.length) { box.append($("p", { class: "note" }, "No rows.")); return box; }
  let rows = t.rows.slice(), sortI = -1, asc = true, filter = "", limit = 200;
  const ci = ["contig", "CHROM", "rname"].map(c => t.columns.indexOf(c)).find(i => i >= 0);
  const si = ["start", "POS", "startpos"].map(c => t.columns.indexOf(c)).find(i => i >= 0);
  const ei = ["end", "endpos"].map(c => t.columns.indexOf(c)).find(i => i >= 0);
  const linkable = ci !== undefined && si !== undefined && (opts.sample || opts.sampleCol !== undefined);
  const tools = $("div", { class: "tools" });
  const info = $("span");
  if (t.rows.length > 15) tools.append($("input", { placeholder: "filter rows…", oninput: e => { filter = e.target.value.toLowerCase(); draw(); } }));
  tools.append(info);
  const wrapTbl = $("div", { class: "tbl" });
  box.append(tools, wrapTbl);
  const longCols = new Set(t.columns.map((c, i) => /features|annotation|description|product|verdict/i.test(c) ? i : -1));
  function draw() {
    let r = filter ? rows.filter(x => x.join("\t").toLowerCase().includes(filter)) : rows;
    const shown = r.slice(0, limit);
    const tb = $("tbody");
    for (const row of shown) {
      const tr = $("tr", linkable ? { class: "go", title: "show on the coverage plot" } : {});
      row.forEach((v, i) => {
        let td = longCols.has(i) ? longCell(v) : $("td", {}, v);
        if (opts.sampleCol === i && D.samples[v]) td = $("td", {}, $("a", { onclick: ev => { ev.stopPropagation(); go("sample", v); } }, v));
        if (opts.heat && opts.heat(i, v)) td.className = opts.heat(i, v);
        tr.append(td);
      });
      if (linkable) tr.addEventListener("click", () => {
        const s = opts.sample || row[opts.sampleCol], a = num(row[si]), b = ei !== undefined ? num(row[ei]) : a;
        goTo(s, row[ci], a, b ?? a);
      });
      tb.append(tr);
    }
    const thead = $("thead", {}, $("tr", {}, t.columns.map((c, i) => $("th", { onclick: () => sortBy(i) }, c))));
    wrapTbl.replaceChildren($("table", {}, thead, tb));
    info.replaceChildren(`${r.length} row(s)` + (t.total > t.rows.length ? ` (first ${t.rows.length} of ${t.total} kept; full table on disk)` : ""));
    if (r.length > limit) info.append(" · ", $("a", { onclick: () => { limit = Infinity; draw(); } }, `show all ${r.length}`));
  }
  function sortBy(i) {
    asc = sortI === i ? !asc : true; sortI = i;
    const numeric = rows.every(r => r[i] === "" || r[i] === "NA" || num(r[i]) !== null);
    rows.sort((x, y) => { const a = x[i], b = y[i]; const c = numeric ? (num(a) ?? -Infinity) - (num(b) ?? -Infinity) : String(a).localeCompare(String(b));
      return asc ? c : -c; });
    draw();
  }
  draw(); return box;
}

// long text (e.g. every gene a region hits): 2 lines, click to expand / collapse
function longCell(v) {
  const items = String(v).split(";").map(x => x.trim()).filter(Boolean);
  const clip = $("div", { class: "clip", title: "click to show all" }, items.join("; "));
  const td = $("td", { class: "longtxt" }, clip);
  if (items.length > 3) td.append($("span", { class: "more" }, ` ${items.length} items`));
  const toggle = ev => { ev.stopPropagation(); clip.classList.toggle("open"); };
  clip.addEventListener("click", toggle);
  td.querySelector(".more")?.addEventListener("click", toggle);
  return td;
}

function block(b, sample) {
  const out = [$("h4", {}, b.title)];
  if (b.note) out.push($("p", { class: "note" }, b.note));
  if (b.table) out.push(table(b.table, { sample }));
  if (b.text) out.push($("pre", {}, b.text));
  if (b.links) out.push($("p", {}, b.links.map(([href, label]) => $("a", { href, target: "_blank" }, label)).flatMap((a, i) => i ? [" · ", a] : [a])));
  if (b.images) out.push(gallery(b.images));
  return out;
}

function gallery(imgs) {
  return $("div", { class: "gallery" }, imgs.map(([src, cap]) =>
    $("figure", {}, $("img", { src, alt: cap, loading: "lazy", onclick: () => window.open().document.write(`<img src="${src}" style="max-width:none">`) }),
      $("figcaption", {}, cap))));
}

// combine one block title across samples (adds a "sample" column)
function combined(secKey, title) {
  const cols = [], rows = []; let total = 0; const notes = [];
  for (const s of D.order) {
    const b = (D.samples[s].blocks[secKey] || []).find(x => x.title === title);
    if (!b) continue;
    if (!b.table) { if (b.note) notes.push(`${s}: ${b.note}`); continue; }
    if (!cols.length) cols.push("sample", ...b.table.columns);
    for (const r of b.table.rows) rows.push([s, ...r]);
    total += b.table.total;
  }
  return { table: cols.length ? { columns: cols, rows, total: Math.max(total, rows.length) } : null, notes };
}

// ---------------------------------------------------------------- coverage plot
const PLOTS = {};
const TRACK = { depth: 120, ev: 18, gene: 30, axis: 18 };
const EVCOL = { lowcov: "--lowedge", highcov: "--highedge", cassette: "--cas", hotspot: "--hot", breakpoint: "--bp", indel: "--indel", variant: "--var" };

function plot(sample, contig, cassette = false) {
  const smp = D.samples[sample], dp = (cassette ? smp.cdepth : smp.depth)[contig];
  const genes = (cassette ? D.cassette_genes || {} : D.genes[D.ref_of[sample]] || {})[contig] || [];
  const ev = cassette ? [] : smp.events.filter(e => e[1] === contig);
  const sorted = dp.v.slice().sort((a, b) => a - b);
  const yMax = Math.max(1, sorted[Math.floor(sorted.length * 0.99)] * 1.15 || 1);
  const median = sorted[Math.floor(sorted.length / 2)] || 0;
  const st = { a: 1, b: dp.len };
  const cv = $("canvas");
  const pos = $("b"), reset = $("button", { onclick: () => { st.a = 1; st.b = dp.len; draw(); } }, "reset"),
        out = $("button", { onclick: () => zoom((st.a + st.b) / 2, (st.b - st.a) * 3) }, "zoom out");
  const box = $("div", { class: "plot" },
    $("div", { class: "bar" }, $("b", {}, contig), `${fmtBp(dp.len)} · median depth ${median}× · bin ${dp.bin} bp`,
      $("span", { style: "margin-left:auto" }, pos), out, reset), cv);
  const H = TRACK.depth + TRACK.ev + TRACK.gene + TRACK.axis;
  let W = 800, drag = null;
  function zoom(c, w) { w = Math.max(w, Math.min(200, dp.len)); st.a = Math.max(1, Math.round(c - w / 2)); st.b = Math.min(dp.len, Math.round(st.a + w)); draw(); }
  const X = p => (p - st.a) / (st.b - st.a + 1) * W;
  const P = x => st.a + x / W * (st.b - st.a + 1);
  function draw(sel) {
    W = cv.clientWidth || 800; const r = devicePixelRatio || 1;
    cv.width = W * r; cv.height = H * r; cv.style.height = H + "px";
    const g = cv.getContext("2d"); g.setTransform(r, 0, 0, r, 0, 0); g.clearRect(0, 0, W, H);
    pos.textContent = `${contig}:${st.a.toLocaleString()}-${st.b.toLocaleString()} (${fmtBp(st.b - st.a + 1)})`;
    // low-coverage bands
    for (const e of ev) if ((e[0] === "lowcov" || e[0] === "highcov") && e[3] >= st.a && e[2] <= st.b) {
      g.fillStyle = css(e[0] === "lowcov" ? "--low" : "--high"); g.fillRect(X(e[2]), 0, Math.max(2, X(e[3] + 1) - X(e[2])), TRACK.depth); }
    // depth (one value per pixel column: mean of the bins under it)
    g.fillStyle = css("--depth"); g.beginPath(); g.moveTo(0, TRACK.depth);
    for (let x = 0; x < W; x++) {
      const p0 = P(x), p1 = P(x + 1);
      let i0 = Math.floor((p0 - 1) / dp.bin), i1 = Math.max(i0, Math.floor((p1 - 1) / dp.bin));
      i0 = Math.max(0, i0); i1 = Math.min(dp.v.length - 1, i1);
      let s = 0, n = 0; for (let i = i0; i <= i1; i++) { s += dp.v[i]; n++; }
      const v = n ? s / n : 0;
      g.lineTo(x, TRACK.depth - Math.min(v, yMax) / yMax * (TRACK.depth - 6));
    }
    g.lineTo(W, TRACK.depth); g.closePath(); g.fill();
    g.fillStyle = css("--muted"); g.font = "10px system-ui"; g.fillText(`${Math.round(yMax)}×`, 2, 10);
    g.strokeStyle = css("--line"); g.beginPath(); g.moveTo(0, TRACK.depth + .5); g.lineTo(W, TRACK.depth + .5); g.stroke();
    // events lane
    const y0 = TRACK.depth + 3;
    for (const e of ev) {
      if (e[3] < st.a || e[2] > st.b) continue;
      const x0 = X(e[2]), x1 = Math.max(x0 + 2, X(e[3] + 1));
      g.fillStyle = css(EVCOL[e[0]]);
      if (e[0] === "variant") { g.fillRect(x0 - 0.5, y0, 2, 8); }
      else if (e[0] === "breakpoint") { g.globalAlpha = .8; g.fillRect(x0 - 1, 0, 2, TRACK.depth); g.globalAlpha = 1; g.fillRect(x0 - 2, y0, 4, TRACK.ev - 6); }
      else if (e[0] === "lowcov" || e[0] === "highcov") { g.fillRect(x0, y0 + 9, x1 - x0, 4); }
      else { g.fillRect(x0, y0, x1 - x0, TRACK.ev - 6);
        if (e[0] === "cassette" && x1 - x0 > 60) { g.fillStyle = "#fff"; g.font = "10px system-ui"; g.fillText(e[4], x0 + 3, y0 + 9, x1 - x0 - 6); } }
    }
    // genes
    const gy = TRACK.depth + TRACK.ev + 4, span = st.b - st.a;
    if (span <= 300000) {
      g.font = "10px system-ui";
      for (const ge of genes) {
        if (ge[1] < st.a || ge[0] > st.b) continue;
        const x0 = X(ge[0]), x1 = X(ge[1] + 1), w = Math.max(1, x1 - x0);
        g.fillStyle = css("--gene"); g.fillRect(x0, gy + (ge[2] > 0 ? 0 : 7), w, 6);
        if (w > g.measureText(ge[3]).width + 4) { g.fillStyle = css("--ink"); g.fillText(ge[3], x0 + 2, gy + 24); }
      }
    } else { g.fillStyle = css("--muted"); g.font = "10px system-ui"; g.fillText(`${genes.length} genes (zoom below 300 kb to see them)`, 2, gy + 12); }
    // axis
    const ay = H - TRACK.axis + 12, step = niceStep(span / Math.max(2, W / 110));
    g.fillStyle = css("--muted"); g.font = "10px system-ui";
    for (let p = Math.ceil(st.a / step) * step; p <= st.b; p += step) { const x = X(p); g.fillRect(x, ay - 12, 1, 4); g.fillText(fmtBp(p), x + 2, ay); }
    if (sel) { g.fillStyle = "rgba(47,93,138,.18)"; g.fillRect(Math.min(sel[0], sel[1]), 0, Math.abs(sel[1] - sel[0]), H); }
  }
  function niceStep(raw) { const p = Math.pow(10, Math.floor(Math.log10(raw))); for (const m of [1, 2, 5, 10]) if (m * p >= raw) return m * p; return 10 * p; }
  const tip = document.getElementById("tip");
  cv.addEventListener("mousedown", e => { drag = [e.offsetX, e.offsetX]; });
  cv.addEventListener("mousemove", e => {
    if (drag) { drag[1] = e.offsetX; draw(drag); return; }
    const p = Math.round(P(e.offsetX)), i = Math.min(dp.v.length - 1, Math.floor((p - 1) / dp.bin));
    const here = ev.filter(x => x[2] - (st.b - st.a) / W * 3 <= p && x[3] + (st.b - st.a) / W * 3 >= p);
    const gs = genes.filter(x => x[0] <= p && x[1] >= p);
    tip.innerHTML = `<b>${contig}:${p.toLocaleString()}</b><br>depth ≈ ${dp.v[i]}×` +
      gs.map(x => `<br>gene <b>${x[3]}</b> ${x[4] ? "– " + x[4] : ""}`).join("") +
      here.map(x => `<br><span style="color:${css(EVCOL[x[0]])}">■</span> ${x[0]} ${x[2]}-${x[3]} ${x[4]}`).join("");
    tip.style.display = "block"; tip.style.left = Math.min(e.clientX + 14, innerWidth - 370) + "px"; tip.style.top = (e.clientY + 14) + "px";
  });
  cv.addEventListener("mouseleave", () => { tip.style.display = "none"; });
  addEventListener("mouseup", e => {
    if (!drag) return; const [x0, x1] = drag.sort((a, b) => a - b); drag = null;
    if (x1 - x0 > 4) { const a = P(x0), b = P(x1); zoom((a + b) / 2, b - a); } else draw();
  });
  cv.addEventListener("dblclick", () => { st.a = 1; st.b = dp.len; draw(); });
  new ResizeObserver(() => draw()).observe(cv);
  const api = { box, show(a, b) { const w = Math.max(b - a + 1, 500) * 3; zoom((a + b) / 2, w); } };
  if (!cassette) PLOTS[sample + "\u0000" + contig] = api;
  return api;
}

function plots(sample) {
  const dp = D.samples[sample].depth, names = Object.keys(dp);
  if (!names.length) return [$("p", { class: "note" }, "No depth profile.")];
  return names.map(c => plot(sample, c).box);
}

function goTo(sample, contig, a, b) {
  go("sample", sample);
  const p = PLOTS[sample + "\u0000" + contig];
  if (!p) return;
  p.show(a, b); p.box.scrollIntoView({ behavior: "smooth", block: "center" });
}

// ---------------------------------------------------------------- views
function summaryRow(s) {
  const t = D.summary; if (!t) return {};
  const r = t.rows.find(x => x[0] === s); return r ? Object.fromEntries(t.columns.map((c, i) => [c, r[i]])) : {};
}

function tiles(s) {
  const r = summaryRow(s), smp = D.samples[s];
  const big = x => x == null || x === "NA" ? "–" : (num(x) !== null && num(x) >= 1e4 ? Number(num(x)).toLocaleString() : x);
  const cas = r.cassette_insertions && r.cassette_insertions !== "NA" ? (r.cassette_insertions === "none" ? "0" : String(r.cassette_insertions.split(";").length)) : "NA";
  const T = [["Platform", smp.platform], ["Reference", smp.reference], ["Filtered reads", big(r.filt_reads)],
    ["Approx. coverage", r.coverage_approx && r.coverage_approx !== "NA" ? r.coverage_approx + "×" : "–"],
    ["Mapped", r.pct_mapped && r.pct_mapped !== "NA" ? r.pct_mapped + " %" : "–"], ["Mean depth", big(r.mean_depth)],
    ["Min replicon breadth", r.min_replicon_breadth && r.min_replicon_breadth !== "NA" ? r.min_replicon_breadth + " %" : "–"],
    ["Variants", big(r.detected_mutations)], ["Low-coverage regions", big(r.low_cov_regions)],
    ["High-coverage regions", big(r.high_cov_regions)],
    ["Large indels", big(r.large_indels)], ["Breakpoints", big(r.breakpoints)], ["Cassette sites", big(cas)],
    ["Multimer", big(r.multimer)], ["Assembly contigs", big(r.assembly_contigs)]];
  return $("div", { class: "tiles" }, T.filter(([k, v]) => (v !== "NA" && v !== "–") || ["Platform", "Reference"].includes(k))
    .map(([k, v]) => $("div", { class: "tile" }, $("b", {}, v), $("span", {}, k))));
}

function sampleView(s) {
  const smp = D.samples[s], out = [$("h2", {}, s), $("p", { class: "desc" }, `${smp.platform} · reference ${smp.reference} · reads: ${smp.reads.join(", ")}`)];
  if (smp.low_diversity) out.push($("div", { class: "banner" }, "⚠ " + smp.low_diversity));
  out.push(tiles(s));
  for (const sec of D.sections) {
    if (sec.key === "overview") continue;
    if (smp.low_diversity && sec.key !== "qc") continue;
    const blocks = smp.blocks[sec.key] || [];
    out.push($("h3", { id: "s-" + sec.key }, sec.title));
    if (sec.key === "coverage") out.push($("p", { class: "desc" }, SEC.coverage.desc), ...plots(s));
    if (sec.key === "cassettes" && smp.cdepth && Object.keys(smp.cdepth).length) {
      out.push($("h4", {}, "Depth along each cassette"));
      for (const c of Object.keys(smp.cdepth)) out.push(plot(s, c, true).box);
    }
    for (const b of blocks) out.push(...block(b, s));
  }
  return out;
}

function toolView(key) {
  const sec = SEC[key], out = [$("h2", {}, sec.title), $("p", { class: "desc" }, sec.desc)];
  if (key === "overview") {
    const low = D.order.filter(x => D.samples[x].low_diversity);
    if (low.length) out.push($("div", { class: "banner" }, `⚠ Low diversity, not analysed (theoretical depth below min_theoretical_depth): ${low.join(", ")}`));
    if (D.summary) out.push(table(D.summary, { sampleCol: 0, heat: (i, v) => {
      const c = D.summary.columns[i], n = num(v);
      if (c === "status") return v === "ok" ? "ok" : "bad";
      if (c === "min_replicon_breadth" && n !== null) return n < 50 ? "bad" : n < 95 ? "warn" : "ok";
      if (c === "pct_mapped" && n !== null) return n < 50 ? "bad" : n < 90 ? "warn" : null;
      return null; } }));
    out.push($("h3", {}, "Analyses"), $("p", {}, Object.entries(D.toggles).map(([k, v]) => $("span", { class: "chip" }, (v ? "✓ " : "✗ ") + k))));
    out.push($("h3", {}, "Run settings (config)"), table({ columns: ["key", "value"], rows: Object.entries(D.config).map(([k, v]) => [k, String(v)]), total: Object.keys(D.config).length }));
    return out;
  }
  if (key === "coverage" && D.coverage_summary) {
    out.push($("h4", {}, "Every sample × replicon (coverage_summary.tsv)"), table(D.coverage_summary, { sampleCol: 0 }));
  }
  if (key === "plasmids") out.push(...plasmidMatrix());
  // the block titles of this section, in first-seen order
  const titles = [];
  for (const s of D.order) for (const b of D.samples[s].blocks[key] || []) if (!titles.includes(b.title)) titles.push(b.title);
  for (const title of titles) {
    const hasTable = D.order.some(s => (D.samples[s].blocks[key] || []).some(b => b.title === title && b.table));
    out.push($("h3", {}, title));
    if (hasTable) {
      const c = combined(key, title);
      if (c.table) out.push(table(c.table, { sampleCol: 0 }));
      if (c.notes.length) out.push($("p", { class: "note" }, c.notes.join(" · ")));
    } else {
      for (const s of D.order) for (const b of (D.samples[s].blocks[key] || []).filter(b => b.title === title)) {
        out.push($("h4", {}, $("a", { onclick: () => go("sample", s) }, s)));
        const rest = block(b, s); rest.shift(); out.push(...rest);
      }
    }
  }
  if (key === "coverage") for (const s of D.order.filter(x => !D.samples[x].low_diversity)) out.push($("h3", {}, $("a", { onclick: () => go("sample", s) }, s)), ...plots(s));
  return out;
}

function plasmidMatrix() {
  const t = D.coverage_summary; if (!t) return [];
  const iS = t.columns.indexOf("sample"), iSrc = t.columns.indexOf("source"), iR = t.columns.indexOf("replicon"),
        iC = t.columns.indexOf("coverage"), iD = t.columns.indexOf("meandepth");
  const rows = t.rows.filter(r => r[iSrc] === "plasmid_panel"); if (!rows.length) return [];
  const pls = [...new Set(rows.map(r => r[iR]))];
  const cell = {}; for (const r of rows) cell[r[iS] + "\u0000" + r[iR]] = r;
  const tb = $("tbody", {}, D.order.map(s => $("tr", {}, $("td", {}, $("a", { onclick: () => go("sample", s) }, s)),
    pls.map(p => { const r = cell[s + "\u0000" + p]; if (!r) return $("td", {}, "–");
      const b = num(r[iC]); return $("td", { class: "heat " + (b >= 95 ? "ok" : b >= 20 ? "warn" : "bad"), title: `${r[iC]} % covered, ${r[iD]}× mean depth` },
        `${Math.round(b)} % · ${Number(num(r[iD])).toFixed(0)}×`); }))));
  return [$("h3", {}, "Presence matrix (breadth · mean depth)"),
    $("div", { class: "tbl" }, $("table", {}, $("thead", {}, $("tr", {}, $("th", {}, "sample"), pls.map(p => $("th", {}, p)))), tb))];
}

// ---------------------------------------------------------------- navigation
function nav() {
  const n = document.getElementById("nav"); n.replaceChildren();
  if (state.mode === "sample") {
    for (const s of D.order) {
      const r = summaryRow(s);
      const sub = [D.samples[s].platform, r.detected_mutations && r.detected_mutations !== "NA" ? r.detected_mutations + " var" : null,
                   r.low_cov_regions && r.low_cov_regions !== "NA" ? r.low_cov_regions + " lowcov" : null].filter(Boolean).join(" · ");
      const small = D.samples[s].low_diversity ? $("small", { class: "lowdiv" }, "LOW DIVERSITY · skipped") : $("small", {}, sub);
      n.append($("a", { class: s === state.sample ? "on" : "", onclick: () => go("sample", s) }, s, small));
    }
  } else for (const sec of D.sections) n.append($("a", { class: sec.key === state.tool ? "on" : "", onclick: () => go("tool", sec.key) }, sec.title));
}

function go(mode, what) {
  const changed = state.mode !== mode || (mode === "sample" ? state.sample !== what : state.tool !== what);
  state.mode = mode; if (mode === "sample") state.sample = what; else state.tool = what;
  try { localStorage.setItem("seqpipe-report-mode", mode); } catch (e) {}
  history.replaceState(null, "", "#" + mode + "=" + encodeURIComponent(what));
  document.getElementById("m-sample").classList.toggle("on", mode === "sample");
  document.getElementById("m-tool").classList.toggle("on", mode === "tool");
  nav();
  if (!changed && document.getElementById("main").childElementCount) return;
  const m = document.getElementById("main");
  m.replaceChildren(...(mode === "sample" ? sampleView(what) : toolView(what)));
  scrollTo(0, 0);
}

document.getElementById("m-sample").onclick = () => go("sample", state.sample);
document.getElementById("m-tool").onclick = () => go("tool", state.tool);
(function start() {
  const h = decodeURIComponent(location.hash.slice(1)), [m, w] = h.split("=");
  if (m === "sample" && D.samples[w]) return go("sample", w);
  if (m === "tool" && SEC[w]) return go("tool", w);
  let saved = null; try { saved = localStorage.getItem("seqpipe-report-mode"); } catch (e) {}
  go(saved === "tool" ? "tool" : "sample", saved === "tool" ? "overview" : D.order[0]);
})();
</script></body></html>
"""


if "snakemake" in globals():
    _sm = globals()["snakemake"]
    with open(str(_sm.log[0]), "w") as _log:
        sys.stderr = _log
        main(list(_sm.params.samples), str(_sm.params.out_dir), dict(_sm.params.toggles),
             bool(_sm.params.embed_igv), dict(_sm.config), str(_sm.output.html),
             list(_sm.params.get("cassettes", [])))
