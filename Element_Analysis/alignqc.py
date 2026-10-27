#!/usr/bin/env python3
"""alignqc -- paired-end (Element AVITI / Illumina) alignment QC pipeline.

For each sample in a sample sheet, against a GenBank reference (e.g. a plasmid,
phage or bacterial genome), the pipeline:

  index     reference.gb -> FASTA + bwa-mem2 / samtools indexes
  align     bwa-mem2 + samtools sort/index/flagstat          (per sample)
  lowcov    coverage summary + low-coverage regions          (per sample)
  indels    large CIGAR indels + soft-clip breakpoints       (per sample)
  annotate  overlap all result tables with GenBank features  (all samples)

Usage (see README.md for the full tutorial):

  python alignqc.py setup-env              # once, on a login node
  python alignqc.py run                    # submit everything to SLURM
  python alignqc.py run --local            # or run everything on this machine
  python alignqc.py run --from-step lowcov # re-run from a later step
  python alignqc.py align --row 2          # run one step for one sample by hand

Every option can also be set through an environment variable of the same name
in UPPER_CASE (e.g. MIN_MAPQ=30), preserving the old config conventions.

Coordinates in all output tables are 0-based half-open (BED-like).
Python >= 3.8; needs biopython + pysam (installed by `setup-env`) and
bwa-mem2 + samtools on PATH.
"""
import argparse
import os
import shlex
import statistics
import subprocess
import sys
from array import array
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT = Path(__file__).resolve()
STEPS = ["index", "align", "lowcov", "indels", "annotate"]

# =============================================================================
# 1. CONFIGURATION
#    (name, type, default, help). One table drives the CLI, the env-var
#    overrides, and the options forwarded to SLURM jobs.
# =============================================================================
OPTIONS = [
    # ---- inputs / outputs
    ("sample_sheet", str, "samples.tsv", "TSV with header: sample<TAB>R1<TAB>R2"),
    ("reference_gb", str, "reference.gb", "GenBank reference (record IDs = contig names)"),
    ("outdir", str, "results", "output directory"),
    # ---- cluster
    ("slurm_account", str, "def-rodrigu1", "SLURM account"),
    ("modules", str, "StdEnv/2023 bwa-mem2 samtools python", "space-separated `module load` list"),
    ("venv", str, str(Path.home() / "envs/alignqc"), "python virtualenv (made by setup-env)"),
    # ---- alignment
    ("platform", str, "ILLUMINA", "read-group PL tag; use ELEMENT for AVITI data"),
    ("min_mapq", int, 20, "minimum mapping quality used in steps lowcov / indels"),
    # ---- low coverage: a base is 'low' if depth < max(low_depth, low_frac * median)
    ("low_depth", int, 10, "absolute depth floor for 'low'"),
    ("low_frac", float, 0.25, "fraction of median depth below which a base is 'low'"),
    ("low_min_len", int, 20, "minimum low-coverage region length (bp)"),
    ("low_merge_gap", int, 10, "merge low regions closer than this (bp)"),
    # ---- indels / breakpoints
    ("min_indel", int, 10, "minimum CIGAR indel size (bp)"),
    ("min_clip", int, 20, "minimum soft-clip length to call a breakpoint (bp)"),
    ("min_support", int, 3, "minimum reads supporting an event"),
    ("cluster_tol", int, 5, "position tolerance when clustering events (bp)"),
    # ---- annotation
    ("skip_types", str, "source,gene", "GenBank feature types to ignore"),
]

# SLURM resources per step: (cpus, mem, time, is_array)
RESOURCES = {
    "index":    (1, "8G",  "00:30:00", False),
    "align":    (8, "16G", "02:00:00", True),
    "lowcov":   (1, "8G",  "00:30:00", True),
    "indels":   (1, "8G",  "01:00:00", True),
    "annotate": (1, "4G",  "00:30:00", False),
}


def add_config_args(parser):
    """Add every OPTIONS entry as --kebab-case flag; default from env var."""
    for name, typ, default, helptext in OPTIONS:
        env = os.environ.get(name.upper())
        parser.add_argument("--" + name.replace("_", "-"), dest=name, type=typ,
                            default=typ(env) if env is not None else default,
                            help=f"{helptext} [{default}]")


def config_to_argv(cfg):
    """Serialise config back to CLI flags (used to forward to SLURM jobs)."""
    out = []
    for name, *_ in OPTIONS:
        out += ["--" + name.replace("_", "-"), str(getattr(cfg, name))]
    return out


# =============================================================================
# 2. SHARED HELPERS
# =============================================================================
def ref_fasta(cfg):
    return Path(cfg.outdir) / "reference" / "reference.fa"


def read_samples(cfg):
    """Parse the sample sheet -> list of dicts with keys sample, R1, R2.
    Blank lines and lines starting with '#' are ignored."""
    rows, header = [], None
    with open(cfg.sample_sheet) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                missing = {"sample", "R1", "R2"} - set(header)
                if missing:
                    sys.exit(f"Sample sheet needs columns sample,R1,R2 (missing {missing}); "
                             "is it tab-separated?")
                continue
            rows.append(dict(zip(header, fields)))
    if not rows:
        sys.exit(f"No samples found in {cfg.sample_sheet}")
    return rows


def pick_rows(cfg, row):
    """Rows to process: one 1-based `row` (CLI or SLURM array index), else all."""
    samples = read_samples(cfg)
    if row is None:
        row = os.environ.get("SLURM_ARRAY_TASK_ID")
    if row is None:
        return samples
    row = int(row)
    if not 1 <= row <= len(samples):
        sys.exit(f"--row {row} out of range (1..{len(samples)})")
    return [samples[row - 1]]


class Progress:
    """Minimal progress bar on stderr (no extra dependency). Redraws in place
    on a terminal; in SLURM logs (not a tty) prints one line per 10%."""
    def __init__(self, label, total):
        self.label, self.total, self.last = label, max(total, 1), -1
        self.tty = sys.stderr.isatty()

    def set(self, n):
        pct = min(100, 100 * n // self.total)
        shown = pct if self.tty else pct // 10 * 10
        if shown != self.last:
            self.last = shown
            print(f"{self.label} [{'#' * (pct // 4):<25}] {pct:3d}%  ({n:,})",
                  end="\r" if self.tty else "\n", file=sys.stderr, flush=True)

    def close(self, n):
        self.set(self.total if n >= self.total else n)
        if self.tty:
            print(file=sys.stderr)


def sh(cmd, **kw):
    """Run a command (list) and fail loudly."""
    print("+ " + " ".join(shlex.quote(str(c)) for c in cmd), file=sys.stderr)
    return subprocess.run([str(c) for c in cmd], check=True, **kw)


# =============================================================================
# 3. PIPELINE STEPS
# =============================================================================
def step_index(cfg, **_):
    """reference.gb -> reference.fa (+ bwa-mem2 and samtools faidx indexes)."""
    from Bio import SeqIO
    fa = ref_fasta(cfg)
    fa.parent.mkdir(parents=True, exist_ok=True)
    # Contig names = record.id (e.g. NC_000913.3); the annotate step matches on it.
    n = SeqIO.write(SeqIO.parse(cfg.reference_gb, "genbank"), fa, "fasta")
    print(f"Wrote {n} sequence(s) to {fa}")
    sh(["bwa-mem2", "index", fa])
    sh(["samtools", "faidx", fa])


def step_align(cfg, row=None, **_):
    """Align paired-end R1/R2 with bwa-mem2, coordinate-sort, index, flagstat."""
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))
    for s in pick_rows(cfg, row):
        name, r1, r2 = s["sample"], s["R1"], s["R2"]
        if not all(os.path.isfile(f) and os.path.getsize(f) for f in (r1, r2)):
            sys.exit(f"Missing/empty FASTQ: {r1} / {r2}")
        bam = Path(cfg.outdir) / "bam" / f"{name}.sorted.bam"
        stats = Path(cfg.outdir) / "stats"
        bam.parent.mkdir(parents=True, exist_ok=True)
        stats.mkdir(parents=True, exist_ok=True)
        tmp = os.environ.get("SLURM_TMPDIR", "/tmp")
        rg = f"@RG\\tID:{name}\\tSM:{name}\\tPL:{cfg.platform}"
        print(f"[{name}] aligning {r1} + {r2}")
        # bwa-mem2 | samtools sort, streamed (no intermediate SAM on disk)
        bwa = subprocess.Popen(
            ["bwa-mem2", "mem", "-t", str(threads), "-R", rg, str(ref_fasta(cfg)), r1, r2],
            stdout=subprocess.PIPE)
        sh(["samtools", "sort", "-@", "4", "-T", f"{tmp}/{name}", "-o", bam, "-"],
           stdin=bwa.stdout)
        bwa.stdout.close()
        if bwa.wait() != 0:
            sys.exit(f"[{name}] bwa-mem2 failed")
        sh(["samtools", "index", bam])
        with open(stats / f"{name}.flagstat.txt", "w") as fh:
            sh(["samtools", "flagstat", bam], stdout=fh)
        print(f"[{name}] done -> {bam}")


def _median_from_hist(hist, n):
    seen, half = 0, n / 2
    for k in sorted(hist):
        seen += hist[k]
        if seen >= half:
            return k


def step_lowcov(cfg, row=None, **_):
    """Per-contig coverage summary + merged low-coverage regions.

    Depth comes from `samtools depth -a` (every position, including zeros) and
    is held in memory at 4 bytes/base: fine for bacteria/plasmids/phage.
    """
    for s in pick_rows(cfg, row):
        name = s["sample"]
        cov = Path(cfg.outdir) / "coverage"
        cov.mkdir(parents=True, exist_ok=True)
        bam = Path(cfg.outdir) / "bam" / f"{name}.sorted.bam"

        # ---- pass 1: read depths, contig -> array (index = pos-1)
        # -aa (not -a): also report contigs with no reads at all, e.g. a lost plasmid
        depth = {}
        with open(f"{ref_fasta(cfg)}.fai") as fai:
            genome_len = sum(int(l.split("\t")[1]) for l in fai)
        prog = Progress(f"[{name}] reading depth", genome_len)
        proc = subprocess.Popen(["samtools", "depth", "-aa", "-Q", str(cfg.min_mapq), str(bam)],
                                stdout=subprocess.PIPE, text=True)
        n = 0
        for n, line in enumerate(proc.stdout, 1):
            if n % 500_000 == 0:
                prog.set(n)
            contig, _pos, d = line.rstrip("\n").split("\t")
            depth.setdefault(contig, array("I")).append(int(d))
        prog.close(n)
        if proc.wait() != 0 or not depth:
            sys.exit(f"[{name}] samtools depth failed or returned nothing")

        hist = Counter()
        for arr in depth.values():
            hist.update(arr)
        median = _median_from_hist(hist, sum(hist.values()))
        # A base is "low" if depth < max(low_depth, low_frac * median depth)
        thr = max(cfg.low_depth, cfg.low_frac * median)
        print(f"[{name}] median depth={median}  low-coverage threshold=<{thr:.1f}")

        # ---- summary table
        with open(cov / f"{name}.coverage_summary.tsv", "w") as out:
            out.write("#contig\tlength\tmean_depth\tmedian_depth\tpct_zero\tpct_low\n")
            for contig, arr in depth.items():
                n, h = len(arr), Counter(arr)
                low = sum(v for k, v in h.items() if k < thr)
                out.write(f"{contig}\t{n}\t{sum(arr)/n:.1f}\t{_median_from_hist(h, n)}\t"
                          f"{100*h[0]/n:.2f}\t{100*low/n:.2f}\n")

        # ---- pass 2: low regions, merging gaps <= low_merge_gap
        with open(cov / f"{name}.low_coverage.tsv", "w") as out:
            out.write("#contig\tstart\tend\tlength\tmean_depth\tmin_depth\n")
            for contig, arr in depth.items():
                regions = []                    # [start, end, sum, n, min]
                for i, d in enumerate(arr):
                    if d >= thr:
                        continue
                    if regions and i - regions[-1][1] <= cfg.low_merge_gap:
                        r = regions[-1]
                        r[1], r[2], r[3], r[4] = i + 1, r[2] + d, r[3] + 1, min(r[4], d)
                    else:
                        regions.append([i, i + 1, d, 1, d])
                for start, end, tot, n, mn in regions:
                    if end - start >= cfg.low_min_len:
                        out.write(f"{contig}\t{start}\t{end}\t{end-start}\t{tot/n:.1f}\t{mn}\n")
        print(f"[{name}] -> {cov}/{name}.low_coverage.tsv")


def _cluster(events, tol):
    """Group position-sorted tuples (pos, ...) whose positions are within tol."""
    events.sort()
    clusters, cur = [], [events[0]]
    for e in events[1:]:
        if e[0] - cur[-1][0] <= tol:
            cur.append(e)
        else:
            clusters.append(cur)
            cur = [e]
    clusters.append(cur)
    return clusters


def step_indels(cfg, row=None, **_):
    """Large indels and structural breakpoints from a sorted, indexed BAM.

    Two evidence types (PE75 reads cannot span big indels, so both are needed):
      1. CIGAR indels >= min_indel bp, clustered by position -> .indels.tsv
      2. Soft-clip pile-ups (>= min_clip bp) = breakpoints from large
         deletions / insertions / rearrangements; n_split counts reads that
         also carry an SA (supplementary alignment) tag      -> .breakpoints.tsv
    """
    import pysam
    for s in pick_rows(cfg, row):
        name = s["sample"]
        outdir = Path(cfg.outdir) / "variants"
        outdir.mkdir(parents=True, exist_ok=True)
        prefix = outdir / name

        indel_events = defaultdict(list)            # (contig, 'INS'|'DEL') -> [(pos, len)]
        clip_events = defaultdict(lambda: [0, 0])   # (contig, pos, side) -> [n_reads, n_split]

        bam = pysam.AlignmentFile(str(Path(cfg.outdir) / "bam" / f"{name}.sorted.bam"), "rb")
        prog = Progress(f"[{name}] scanning reads", bam.mapped)   # from the .bai, instant
        n = 0
        for n, r in enumerate(bam.fetch(), 1):
            if n % 200_000 == 0:
                prog.set(n)
            # unmapped|secondary|qcfail|duplicate|supplementary in one test
            if r.flag & 0xF04 or r.mapping_quality < cfg.min_mapq:
                continue
            ct = r.cigartuples
            if len(ct) == 1:                        # plain e.g. 75M: no indel, no clip
                continue
            contig, ref = r.reference_name, r.reference_start
            for op, ln in ct:
                if op in (0, 7, 8, 3):              # M, =, X, N consume reference
                    ref += ln
                elif op == 1 and ln >= cfg.min_indel:           # insertion
                    indel_events[(contig, "INS")].append((ref, ln))
                elif op == 2:                                   # deletion
                    if ln >= cfg.min_indel:
                        indel_events[(contig, "DEL")].append((ref, ln))
                    ref += ln
            split = 1 if r.has_tag("SA") else 0
            if ct[0][0] == 4 and ct[0][1] >= cfg.min_clip:      # left soft clip
                e = clip_events[(contig, r.reference_start, "L")]
                e[0] += 1; e[1] += split
            if ct[-1][0] == 4 and ct[-1][1] >= cfg.min_clip:    # right soft clip
                e = clip_events[(contig, r.reference_end - 1, "R")]
                e[0] += 1; e[1] += split
        prog.close(n)
        print(f"[{name}] writing indel / breakpoint tables")

        def depth_at(contig, pos):
            return bam.count(contig, max(pos, 0), max(pos, 0) + 1, read_callback="all")

        # ---- indels
        with open(f"{prefix}.indels.tsv", "w") as out:
            out.write("#contig\tstart\tend\ttype\tn_reads\tmedian_len\tmin_len\tmax_len\t"
                      "depth\tfrac_reads\n")
            for (contig, kind), ev in sorted(indel_events.items()):
                for c in _cluster(ev, cfg.cluster_tol):
                    if len(c) < cfg.min_support:
                        continue
                    start = int(statistics.median(x[0] for x in c))
                    lens = [x[1] for x in c]
                    med = int(statistics.median(lens))
                    end = start + (med if kind == "DEL" else 1)
                    d = depth_at(contig, start)
                    out.write(f"{contig}\t{start}\t{end}\t{kind}\t{len(c)}\t{med}\t"
                              f"{min(lens)}\t{max(lens)}\t{d}\t{len(c)/max(d,1):.2f}\n")

        # ---- breakpoints (merge jittered neighbouring positions on the same side)
        by_side = defaultdict(list)
        for (contig, pos, side), (n, ns) in clip_events.items():
            by_side[(contig, side)].append((pos, n, ns))
        with open(f"{prefix}.breakpoints.tsv", "w") as out:
            out.write("#contig\tstart\tend\tside\tn_reads\tn_split\tdepth\tfrac_reads\n")
            for (contig, side), ev in sorted(by_side.items()):
                for c in _cluster(ev, cfg.cluster_tol):
                    n = sum(x[1] for x in c)
                    if n < cfg.min_support:
                        continue
                    pos = max(c, key=lambda x: x[1])[0]         # strongest position
                    d = depth_at(contig, pos)
                    out.write(f"{contig}\t{pos}\t{pos+1}\t{side}\t{n}\t{sum(x[2] for x in c)}\t"
                              f"{d}\t{n/max(d,1):.2f}\n")
        bam.close()
        print(f"[{name}] done -> {prefix}.{{indels,breakpoints}}.tsv")


def annotate_table(gb, tsv_in, tsv_out, skip_types="source,gene"):
    """Append a 'features' column to any TSV whose first 3 columns are
    contig/start/end (0-based half-open, '#' header) by overlapping with
    GenBank features. Reusable on any table from lowcov / indels."""
    from Bio import SeqIO
    skip = set(skip_types.split(","))
    feats = defaultdict(list)                       # contig -> [(start, end, label)]
    for rec in SeqIO.parse(gb, "genbank"):
        for f in rec.features:
            if f.type in skip:
                continue
            q = f.qualifiers
            nm = (q.get("gene") or q.get("locus_tag") or q.get("label") or ["?"])[0]
            prod = (q.get("product") or q.get("note") or [""])[0]
            label = f"{f.type}:{nm}" + (f"({prod})" if prod else "")
            feats[rec.id].append((int(f.location.start), int(f.location.end), label))
    with open(tsv_in) as fin, open(tsv_out, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if line.startswith("#"):
                fout.write(line + "\tfeatures\n")
                continue
            c = line.split("\t")
            contig, s, e = c[0], int(c[1]), int(c[2])
            hits = [lab for fs, fe, lab in feats.get(contig, []) if fs < e and fe > s]
            fout.write(line + "\t" + (";".join(hits) if hits else ".") + "\n")


def step_annotate(cfg, **_):
    """Annotate every low-coverage / indel / breakpoint table (all samples)."""
    out = Path(cfg.outdir) / "annotated"
    out.mkdir(parents=True, exist_ok=True)
    o = Path(cfg.outdir)
    tables = (sorted(o.glob("coverage/*.low_coverage.tsv"))
              + sorted(o.glob("variants/*.indels.tsv"))
              + sorted(o.glob("variants/*.breakpoints.tsv")))
    for f in tables:
        annotate_table(cfg.reference_gb, f, out / f"{f.stem}.annotated.tsv", cfg.skip_types)
    print(f"Annotated {len(tables)} table(s) in {out}/")


# =============================================================================
# 4. ORCHESTRATION (SLURM or local) AND ENVIRONMENT SETUP
# =============================================================================
def step_setup_env(cfg, **_):
    """One-time: python virtualenv with biopython + pysam (Alliance wheels).
    Run on a login node."""
    cmd = (f"module load {cfg.modules} && mkdir -p {shlex.quote(str(Path(cfg.venv).parent))} "
           f"&& virtualenv --no-download {shlex.quote(cfg.venv)} "
           f"&& source {shlex.quote(cfg.venv)}/bin/activate "
           "&& pip install --no-index --upgrade pip && pip install --no-index biopython pysam")
    sh(["bash", "-c", cmd])
    print(f"Environment ready: {cfg.venv}")


def step_run(cfg, local=False, from_step="index", **_):
    """Run the whole pipeline: submit to SLURM with dependencies, or --local."""
    samples = read_samples(cfg)                      # validates the sheet early
    steps = STEPS[STEPS.index(from_step):]
    print(f"Samples: {len(samples)}   Reference: {cfg.reference_gb}   Output: {cfg.outdir}")

    if local:
        for st in steps:
            print(f"=== {st} ===")
            STEP_FUNCS[st](cfg)                      # row=None -> all samples
        return

    Path("logs").mkdir(exist_ok=True)
    env_setup = (f"module load {cfg.modules}; "
                 f"[ -f {shlex.quote(cfg.venv)}/bin/activate ] && source {shlex.quote(cfg.venv)}/bin/activate; ")
    flags = " ".join(shlex.quote(a) for a in config_to_argv(cfg))
    jobs = {}

    def submit(step, deps):
        cpus, mem, time, is_array = RESOURCES[step]
        cmd = ["sbatch", "--parsable", f"--account={cfg.slurm_account}", f"--job-name={step}",
               f"--cpus-per-task={cpus}", f"--mem={mem}", f"--time={time}",
               f"--output=logs/{step}_%{'A_%a' if is_array else 'j'}.out"]
        if is_array:
            cmd.append(f"--array=1-{len(samples)}")
        if deps:
            cmd.append("--dependency=afterok:" + ":".join(deps))
        cmd.append(f"--wrap={env_setup}python {shlex.quote(str(SCRIPT))} {step} {flags}")
        jid = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip()
        print(f"{step:9s}: job {jid}")
        return jid

    # Dependency graph: index -> align -> {lowcov, indels} -> annotate.
    # The array jobs of one step wait for the whole previous step (afterok).
    if "index" in steps:
        jobs["index"] = submit("index", [])
    if "align" in steps:
        jobs["align"] = submit("align", [jobs["index"]] if "index" in jobs else [])
    prev = [jobs["align"]] if "align" in jobs else []
    for st in ("lowcov", "indels"):
        if st in steps:
            jobs[st] = submit(st, prev)
    if "annotate" in steps:
        jobs["annotate"] = submit("annotate", [jobs[s] for s in ("lowcov", "indels") if s in jobs])
    print("Monitor with: squeue -u $USER")


STEP_FUNCS = {"index": step_index, "align": step_align, "lowcov": step_lowcov,
              "indels": step_indels, "annotate": step_annotate}


def main():
    sys.stdout.reconfigure(line_buffering=True)     # SLURM logs show progress live
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("setup-env", help="one-time virtualenv creation (login node)")
    add_config_args(p)
    p = sub.add_parser("run", help="run all steps (SLURM by default)")
    add_config_args(p)
    p.add_argument("--local", action="store_true", help="run on this machine, no SLURM")
    p.add_argument("--from-step", choices=STEPS, default="index",
                   help="start at this step (earlier outputs must already exist)")
    for st in STEPS:
        p = sub.add_parser(st, help=f"run only the '{st}' step")
        add_config_args(p)
        if RESOURCES[st][3]:
            p.add_argument("--row", type=int,
                           help="1-based sample-sheet row (default: SLURM array index, else all)")

    args = ap.parse_args()
    if args.cmd == "setup-env":
        step_setup_env(args)
    elif args.cmd == "run":
        step_run(args, local=args.local, from_step=args.from_step)
    else:
        STEP_FUNCS[args.cmd](args, row=getattr(args, "row", None))


if __name__ == "__main__":
    main()
