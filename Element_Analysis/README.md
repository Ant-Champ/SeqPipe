# alignqc — paired-end alignment QC for a GenBank reference

A single script, `alignqc.py`, that takes paired-end FASTQs (Element AVITI or Illumina — both are
plain R1/R2 FASTQ) and a GenBank reference, and reports for every sample:

* **Coverage**: per-contig summary and the low-coverage regions
* **Large indels**: from CIGAR strings
* **Breakpoints**: from soft-clip pile-ups (deletions, insertions, rearrangements)
* **Annotation**: every region above overlapped with the GenBank features (CDS, promoters, …)

Intended for small references (plasmids, phage, bacterial genomes).

## Quick start (Alliance / Compute Canada cluster)

Copy-paste on a **login node**, from the pipeline folder. Replace `def-<pi>` with your group's account.

```bash
# --- once: install -----------------------------------------------------------
conda deactivate                          # if your prompt shows (base): conda must be OFF
conda config --set auto_activate_base false   # optional: stop conda starting in every shell
module load StdEnv/2023 python
python alignqc.py setup-env               # builds ~/envs/alignqc with biopython + pysam

# --- each analysis -----------------------------------------------------------
cp samples.example.tsv samples.tsv        # then edit: sample<TAB>R1<TAB>R2, absolute paths
                                          # reference.gb must be in this folder

conda deactivate                          # (each new shell, if conda is still active)
module load StdEnv/2023 python
python alignqc.py run --slurm-account def-<pi>   # submits all steps as SLURM jobs

squeue -u $USER                           # watch the jobs
tail -f logs/indels_*.out                 # progress of a step
ls results/annotated/                     # final tables, once 'annotate' has finished
```

You do **not** activate `~/envs/alignqc` yourself for `run`: each SLURM job loads the modules and
activates it on its own.

**Alternative — everything in one interactive job** (never `--local` on a login node):

```bash
conda deactivate
salloc --account=def-<pi> --time=2:00:00 --mem=16G --cpus-per-task=8
module load StdEnv/2023 python bwa-mem2 samtools
source ~/envs/alignqc/bin/activate        # required with --local
python alignqc.py run --local
```

Long steps (`lowcov`, `indels`) show a progress bar — on screen with `--local`, and one line per 10 %
in `logs/` for SLURM jobs.

The rest of this README explains each part in detail.

**Contents**

1. [Pipeline overview](#1-pipeline-overview)
2. [Installation](#2-installation)
3. [Preparing your data](#3-preparing-your-data)
4. [Running the pipeline](#4-running-the-pipeline)
5. [Parameters](#5-parameters)
6. [Output](#6-output)
7. [Troubleshooting](#7-troubleshooting)
8. [Extending](#8-extending)
9. [Provenance](#9-provenance)

---

## 1. Pipeline overview

```
reference.gb ──► index ──► align ──┬─► lowcov ──┐
                  (once)  (per sample)          ├─► annotate ──► results/annotated/
samples.tsv ───────────────────────┴─► indels ──┘     (all samples)
```

| Step       | What it does                                                            | Runs        |
|------------|-------------------------------------------------------------------------|-------------|
| `index`    | GenBank → FASTA, `bwa-mem2 index`, `samtools faidx`                      | once        |
| `align`    | `bwa-mem2 mem` → `samtools sort` → index + flagstat                      | per sample  |
| `lowcov`   | `samtools depth -a`; flags bases below `max(LOW_DEPTH, LOW_FRAC×median)`; merges nearby ones | per sample |
| `indels`   | CIGAR indels ≥ `MIN_INDEL`, soft-clips ≥ `MIN_CLIP`, clustered by position | per sample |
| `annotate` | adds a `features` column to every table from `lowcov` / `indels`         | once        |

Large indels and breakpoints are found by two kinds of evidence because 75 bp reads cannot span big
indels: reads that *contain* a ≥10 bp indel in their CIGAR, and reads that are *soft-clipped* at the
junction. `n_split` counts clipped reads that also have an `SA` tag (a supplementary alignment on the
other side of the junction), which is the strongest evidence of a true breakpoint.

---

## 2. Installation

The pipeline needs:

| Component | Used by | Where it comes from |
|-----------|---------|---------------------|
| Python ≥ 3.8 | everything | cluster module / system |
| `biopython` (imported as `Bio`) | `index`, `annotate` | pip |
| `pysam` | `indels` | pip |
| `bwa-mem2` | `index`, `align` | cluster module / conda |
| `samtools` | `index`, `align`, `lowcov` | cluster module / conda |
| SLURM | `run` (optional) | the cluster |

> **Windows is not supported natively** — `pysam`, `bwa-mem2` and `samtools` have no Windows builds.
> Use the cluster, or WSL (Ubuntu) on Windows and follow 2B.

### 2A. On a Digital Research Alliance (Compute Canada) cluster — recommended

Do this once, on a **login node** (Narval, Rorqual, Fir, Nibi, Béluga, …). Do **not** use conda on
Alliance clusters; everything comes from modules and the Alliance's prebuilt Python wheels.

**1. Copy the pipeline to the cluster** (from your own computer):

```bash
scp -r alignqc/ <username>@narval.alliancecan.ca:~/projects/def-<pi>/<username>/
```

or `git clone` it there if it lives in a repository.

**2. Check the required modules exist** on that cluster:

```bash
module load StdEnv/2023
module spider bwa-mem2
module spider samtools
```

If `module spider` lists specific versions or says another module must be loaded first, adjust the
`--modules` option accordingly (see step 4 and [Parameters](#5-parameters)).

**3. Create the Python environment.** `setup-env` loads the modules, creates a virtualenv in
`~/envs/alignqc` and installs `biopython` + `pysam` from the Alliance wheels:

```bash
cd ~/projects/def-<pi>/<username>/alignqc
module load StdEnv/2023 python        # so that `python` exists to run setup-env
python alignqc.py setup-env
```

<details>
<summary>Manual equivalent (if <code>setup-env</code> fails)</summary>

```bash
module load StdEnv/2023 python bwa-mem2 samtools
virtualenv --no-download ~/envs/alignqc
source ~/envs/alignqc/bin/activate
pip install --no-index --upgrade pip
pip install --no-index biopython pysam
```
</details>

**4. Set your SLURM account.** The default is `def-rodrigu1`. If yours differs, either pass
`--slurm-account def-<pi>` on every `run`, or set it once in your `~/.bashrc`:

```bash
echo 'export SLURM_ACCOUNT=def-<pi>' >> ~/.bashrc
```

(Any option can be given as an upper-case environment variable — see [Parameters](#5-parameters).)
The same applies to `MODULES` if you had to change the module list in step 2.

**5. Check the installation:**

```bash
module load StdEnv/2023 python bwa-mem2 samtools
source ~/envs/alignqc/bin/activate
python -c "import Bio, pysam; print('biopython', Bio.__version__, '| pysam', pysam.__version__)"
bwa-mem2 version
samtools --version | head -1
```

All four lines must print a version. You do **not** need to keep the environment activated to submit
jobs: `run` loads the modules and activates the virtualenv inside each SLURM job.

### 2B. On your own Linux / macOS / WSL machine

Use conda (via [Miniforge](https://github.com/conda-forge/miniforge)):

```bash
conda create -n alignqc -c conda-forge -c bioconda python biopython pysam bwa-mem2 samtools
conda activate alignqc
```

Then always run with `--local` (no SLURM; see [Running](#4-running-the-pipeline)). Skip `setup-env`.

---

## 3. Preparing your data

You need two input files in the pipeline folder (or anywhere, passing their paths as options):

### The reference: `reference.gb`

A GenBank file (one or more records). Contig names in all outputs are the GenBank `record.id`
(the accession.version, e.g. `NC_000913.3`) — this is also how features are matched during
annotation. Use another path with `--reference-gb`.

### The sample sheet: `samples.tsv`

Start from the example:

```bash
cp samples.example.tsv samples.tsv
```

and edit it: **tab-separated**, header required, one row per sample (blank lines and `#` lines are
ignored):

```
sample	R1	R2
strainA	/project/def-pi/me/fastq/A_R1.fastq.gz	/project/def-pi/me/fastq/A_R2.fastq.gz
strainB	/project/def-pi/me/fastq/B_R1.fastq.gz	/project/def-pi/me/fastq/B_R2.fastq.gz
```

* Use **absolute paths** to the FASTQs (gzipped or not).
* Sample names become file names — no spaces or slashes.
* Many editors turn tabs into spaces; if you get `Sample sheet needs columns sample,R1,R2`,
  check with `cat -A samples.tsv` (tabs show as `^I`).

On a cluster, keep raw FASTQs in `/project` and consider writing results to `/scratch`
(`--outdir /scratch/<username>/alignqc_results`), then copy the `annotated/` tables back.

---

## 4. Running the pipeline

Always run from the pipeline folder (SLURM logs are written to `./logs/`).

### Whole pipeline

```bash
python alignqc.py run              # cluster: submits all steps to SLURM with dependencies
python alignqc.py run --local      # one machine: runs all steps sequentially
```

On the cluster `run` only submits jobs and returns immediately; it can be run from the login node
(load `python` first: `module load StdEnv/2023 python`). It prints one job ID per step. Then:

```bash
squeue -u $USER                    # running / pending jobs
ls logs/                           # logs/<step>_<jobid>.out (array jobs: <step>_<jobid>_<row>.out)
sacct -j <jobid>                   # state and resource usage of a finished job
```

When the `annotate` job has finished, results are in `results/annotated/` (see [Output](#6-output)).

For **Element AVITI** data add `--platform ELEMENT` (it only sets the read-group `PL` tag).

### Running parts of the pipeline

```bash
python alignqc.py run --from-step lowcov   # BAMs already exist; redo lowcov, indels, annotate
python alignqc.py run --from-step annotate # only re-annotate

python alignqc.py align --row 2            # one step, one sample (row 2 of the sheet), in this process
python alignqc.py lowcov                   # one step, all samples, in this process
```

`--from-step` submits a fresh dependency chain starting there; it does **not** check that the earlier
outputs exist. Single-step commands run directly in your shell, so on a cluster use them only inside
an interactive job (`salloc --account=def-<pi> --time=1:00:00 --mem=8G --cpus-per-task=4`) with the
modules loaded and the virtualenv activated — not on a login node.

---

## 5. Parameters

Every parameter is a command-line flag *and* an upper-case environment variable:

```bash
python alignqc.py run --min-mapq 30 --low-depth 20 --outdir run2
MIN_SUPPORT=5 SAMPLE_SHEET=other.tsv python alignqc.py run --local
python alignqc.py run --help             # full list with defaults
```

| Flag | Default | Meaning |
|------|---------|---------|
| `--sample-sheet` | `samples.tsv` | sample sheet |
| `--reference-gb` | `reference.gb` | GenBank reference |
| `--outdir` | `results` | output directory |
| `--slurm-account` | `def-rodrigu1` | SLURM account — **change to yours** |
| `--modules` | `StdEnv/2023 bwa-mem2 samtools python` | `module load` list (cluster-specific) |
| `--venv` | `~/envs/alignqc` | virtualenv created by `setup-env` |
| `--platform` | `ILLUMINA` | read-group platform; `ELEMENT` for AVITI |
| `--min-mapq` | 20 | min mapping quality for `lowcov` and `indels` |
| `--low-depth` | 10 | absolute depth floor for "low" |
| `--low-frac` | 0.25 | fraction of median depth for "low" |
| `--low-min-len` | 20 | min low-coverage region length (bp) |
| `--low-merge-gap` | 10 | merge low regions closer than this (bp) |
| `--min-indel` | 10 | min CIGAR indel size (bp) |
| `--min-clip` | 20 | min soft-clip length for a breakpoint (bp) |
| `--min-support` | 3 | min reads supporting an event |
| `--cluster-tol` | 5 | position tolerance when clustering (bp) |
| `--skip-types` | `source,gene` | GenBank feature types ignored when annotating |

SLURM resources per step (CPUs, memory, time) are in the `RESOURCES` table near the top of
`alignqc.py`:

| Step | CPUs | Memory | Time | Array (one task per sample) |
|------|------|--------|------|-----------------------------|
| `index` | 1 | 8G | 0:30 | no |
| `align` | 8 | 16G | 2:00 | yes |
| `lowcov` | 1 | 8G | 0:30 | yes |
| `indels` | 1 | 8G | 1:00 | yes |
| `annotate` | 1 | 4G | 0:30 | no |

Raise them there if jobs hit `TIMEOUT` or `OUT_OF_MEMORY` (check with `sacct -j <jobid>`).

---

## 6. Output

```
results/
  reference/reference.fa(.fai, bwa-mem2 index files)
  bam/<sample>.sorted.bam(.bai)
  stats/<sample>.flagstat.txt
  coverage/<sample>.coverage_summary.tsv
  coverage/<sample>.low_coverage.tsv
  variants/<sample>.indels.tsv
  variants/<sample>.breakpoints.tsv
  annotated/<sample>.{low_coverage,indels,breakpoints}.annotated.tsv   <- start here
logs/<step>_<jobid>.out                                                <- SLURM logs
```

All coordinates are **0-based, half-open** (BED-like). All tables start with a `#` header line.

| Table | Columns |
|-------|---------|
| `coverage_summary` | contig, length, mean_depth, median_depth, pct_zero, pct_low |
| `low_coverage` | contig, start, end, length, mean_depth, min_depth |
| `indels` | contig, start, end, type (`INS`/`DEL`), n_reads, median_len, min_len, max_len, depth, frac_reads |
| `breakpoints` | contig, start, end, side (`L`/`R` soft clip), n_reads, n_split, depth, frac_reads |
| `annotated/*` | the above + `features` (`type:name(product)` joined by `;`, or `.` if none) |

**Reading the results**

* Check `stats/<sample>.flagstat.txt` first: a low "properly paired" or mapped percentage means
  the wrong reference or contaminated reads.
* `pct_low` in the coverage summary is a quick per-sample quality check.
* `frac_reads` ≈ 1 → the whole population carries the event (e.g. a clonal deletion); ≈ 0.1 → a minor subpopulation.
* A real deletion typically shows an `L` and an `R` breakpoint on either side of a `low_coverage` region and/or a `DEL` entry.
* Low-coverage with `min_depth` 0 over a feature suggests that feature is deleted; GC-extreme regions show shallow but non-zero dips.

---

## 7. Troubleshooting

| Symptom | Cause / fix |
|---------|-------------|
| `python: command not found` (cluster) | `module load StdEnv/2023 python` first |
| `ModuleNotFoundError: No module named 'Bio'` / `'pysam'` | Environment not created or not active. `setup-env` builds it but does **not** activate it in your shell: with `--local` run `conda deactivate` then `source ~/envs/alignqc/bin/activate` first (or use `run` without `--local`, which activates it in each job); on your own machine `conda activate alignqc` |
| Prompt shows `(base)` on the cluster | Conda is active and can shadow the module Python: `conda deactivate` (and `conda config --set auto_activate_base false` to make it permanent) |
| A step seems stuck | `lowcov` / `indels` print a progress bar (`tail -f logs/<step>_*.out` for SLURM jobs); with no progress for a long time, check the job with `sacct -j <jobid>` |
| `pip install` tries to download / fails on the cluster | Use `pip install --no-index …` (Alliance wheels), never plain `pip install` |
| `module: command not found` | Not on a cluster: use `--local` with the tools already on `PATH` (2B) |
| `module load` error for `bwa-mem2` / `samtools` | Name/version differs on this cluster: check `module spider <name>` and pass `--modules "…"` |
| `sbatch: error: Invalid account` | Set `--slurm-account` / `SLURM_ACCOUNT` to your group's account (`def-<pi>`) |
| `bwa-mem2: command not found` / `samtools: …` | Tools not loaded: `module load bwa-mem2 samtools` (cluster) or install them (2B) |
| `Sample sheet needs columns sample,R1,R2` | File isn't tab-separated or the header is missing |
| `Missing/empty FASTQ` | Wrong path in the sheet; use absolute paths |
| Features all `.` in annotated tables | FASTA contig name ≠ GenBank record id; re-run `index` after changing the reference |
| `--from-step` job fails immediately | The earlier step's output doesn't exist yet |
| Later jobs stay pending with `DependencyNeverSatisfied` | An earlier step failed: read its log in `logs/`, `scancel` the rest, fix, re-run with `--from-step` |
| Job state `TIMEOUT` / `OUT_OF_MEMORY` | Raise that step's entry in `RESOURCES` |

---

## 8. Extending

* **Annotate another table**: `annotate_table(gb, tsv_in, tsv_out)` in `alignqc.py` works on any TSV whose
  first three columns are contig/start/end.
* **Add a step**: write a `step_xxx(cfg, row=None, **_)` function, register it in `STEPS`,
  `STEP_FUNCS`, `RESOURCES` and wire its dependencies in `step_run`.
* **Large genomes**: `lowcov` keeps the depth of every base in memory (4 bytes/base) — fine up to a few Mb,
  switch to streaming for anything much larger.

## 9. Provenance

This replaces the earlier multi-file prototype (`00_config.sh`, `01_…` – `05_…`, `run_pipeline.sh`, and
the helper `.py` files), kept unchanged in `legacy/`. Logic and defaults are the same; the shell
scripts, config file and helper scripts are merged into `alignqc.py`.
