# Sequencing analysis pipeline: ONT + Element / Illumina

One Snakemake pipeline for **Oxford Nanopore (ONT) long reads** and **paired-end short reads
(Element AVITI or Illumina)**, mapped against **GenBank references** (plasmids, phages,
bacterial genomes). ONT and short-read samples can be mixed in one run. Each analysis can be
switched on or off in `config.yaml`.

It merges the former `ONT Analysis` (NanoForge) pipeline and the `Element_Analysis` (alignqc)
pipeline. It is meant to run **interactively inside a `salloc` allocation**, not with `sbatch`.

---

## Quick reminder (when you come back after a few months)

Run these on the cluster. Replace `def-rodrigu1` if your account differs.

```bash
# 1. login node: go to the pipeline and activate the driver environment
cd ~/Pipeline                       # wherever you copied this folder
conda activate seqpipe

# 2. login node: prepare a run folder
mkdir -p ~/runs/MY_RUN && cd ~/runs/MY_RUN
cp ~/Pipeline/template/config.yaml ~/Pipeline/template/sample_sheet.tsv .
nano sample_sheet.tsv               # one row per sample (see "The sample sheet")
nano config.yaml                    # paths + which analyses to run

# 3. login node: check (dry run: lists the jobs, runs nothing)
~/Pipeline/run.sh config.yaml -n

# 4. get a compute allocation inside tmux, then run
tmux new -s pipeline
salloc --account=def-rodrigu1 --time=6:00:00 --cpus-per-task=16 --mem=64G
~/Pipeline/run.sh ~/runs/MY_RUN/config.yaml

# 5. results: start with <out_dir>/summary.tsv
```

If the connection drops, log back in to **the same login node** and run `tmux attach -t pipeline`.
If the allocation ends before the run finishes, start a new `salloc` and run the same command
again. Finished steps are not redone.

---

## Contents

1. [What the pipeline does](#1-what-the-pipeline-does)
2. [Installation (once)](#2-installation-once)
3. [Preparing a run](#3-preparing-a-run)
4. [Running](#4-running)
5. [Output](#5-output)
6. [Configuration reference](#6-configuration-reference)
7. [Troubleshooting](#7-troubleshooting)
8. [Changing the pipeline](#8-changing-the-pipeline)
9. [Provenance](#9-provenance)

---

## 1. What the pipeline does

```
               ┌─ ONT: filtlong ─┐   ┌─ ONT: minimap2 ─┐
FASTQ ─ stage ─┤                 ├───┤                 ├─ sorted.bam ─┬─ variants (medaka | bcftools) ─ GATC + gene annotation
               └─ SR : fastp ────┘   └─ SR : bwa-mem2 ─┘              ├─ coverage per replicon + low-coverage regions ─ annotation
GenBank ─ FASTA + indexes ────────────────────────────────────────────┤─ large indels + breakpoints ─ annotation
                                                                      ├─ multimers (ONT only)
filtered reads ─┬─ de novo assembly (Flye | Unicycler)                └──────────────┐
                └─ plasmid panel screen ─ cassette hotspots ─ IGV snapshots          ├─ summary.tsv
```

| Step | config toggle | ONT | Element / Illumina (SR) |
|---|---|---|---|
| Read filtering | `run_read_filtering` | filtlong (length/quality) | fastp (adapters, 3' quality) |
| Alignment to the sample's reference | always | minimap2 `map-ont` | bwa-mem2 |
| Read statistics, flagstat, coverage per replicon | always | ✓ | ✓ |
| SNPs / small indels → `simple.tsv` | `run_variant_calling` | medaka | bcftools (haploid) |
| GATC site + gene annotation of those variants | `run_gatc_annotation` | ✓ | ✓ |
| Low-coverage regions (deleted genes), annotated | `run_coverage_analysis` | ✓ | ✓ |
| Large indels + breakpoints, annotated | `run_sv_scan` | ✓ | ✓ |
| Multimer / concatemer detection | `run_multimer_analysis` | ✓ | skipped (reads too short) |
| De novo assembly | `run_denovo_assembly` | Flye | Unicycler |
| Coverage of a shared plasmid panel | `run_plasmid_screen` | ✓ | ✓ |
| Where cassette-carrying reads land on the reference | `run_cassette_hotspots` | ✓ | ✓ (both mates kept) |
| IGV pictures of each hotspot | `run_igv_snapshots` | ✓ | ✓ |
| `summary.tsv`, one row per sample | always | ✓ | ✓ |

**How a deletion is found.** It shows up in three tables:
- in `low_coverage`, as a region with `min_depth` 0 over the deleted gene;
- in `breakpoints`, as an `R` breakpoint at its left edge and an `L` breakpoint at its right edge (soft-clipped reads, with `n_split` reads that align on both sides);
- in `indels`, as a `DEL`, when reads span it (long reads).

For a gene replaced by a resistance cassette (Keio-style), the plasmid panel screen also finds
the cassette, and `cassette_hotspots` places it at the gene.

---

## 2. Installation (once)

The only requirement is **conda** working on the cluster. Snakemake builds every tool
(minimap2, medaka, bwa-mem2, …) in its own conda environment the first time.

### 2.1 Copy the pipeline to the cluster

From your own computer, in the folder that contains `Pipeline/` (PowerShell, Git Bash or a Linux
terminal):

```bash
scp -r Pipeline <username>@<cluster-login-address>:~/
```

The pipeline is then in `~/Pipeline` on the cluster. All the commands below assume that location.
If you put it somewhere else, adapt the paths.

Then, on the cluster:

```bash
chmod +x ~/Pipeline/run.sh
```

### 2.2 Check conda and set the channels

On a **login node**:

```bash
conda --version                       # must print a version. If not, load/initialise conda first
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### 2.3 Create the driver environment (snakemake)

```bash
conda env create -f ~/Pipeline/envs/snakemake.yaml       # creates the environment "seqpipe"
conda activate seqpipe
snakemake --version                                       # must print 8.x
```

### 2.4 Build the tool environments (login node, needs internet)

Compute nodes often have **no internet**, so build all environments once, from the login node:

```bash
conda activate seqpipe
~/Pipeline/run.sh setup
```

This builds every per-step environment into `~/Pipeline/.conda_envs/`, and all your runs share
them. It takes 15–40 min, and medaka is the largest. When it finishes it prints
`All conda environments are built in …`.

> Re-run `~/Pipeline/run.sh setup` whenever a file in `envs/` changes. A run that needs a
> missing environment would otherwise try to download it from a compute node.

### 2.5 Test the installation on the bundled data (≈ 5 min)

```bash
tmux new -s pipeline
salloc --account=def-rodrigu1 --time=1:00:00 --cpus-per-task=4 --mem=16G
conda activate seqpipe
~/Pipeline/run.sh test
column -t -s $'\t' ~/Pipeline/test/results/summary.tsv
exit                                     # leave the allocation
```

The test data is synthetic: a 20 kb chromosome plus a 4 kb plasmid. The strain has **2 SNPs in
geneA**, and **geneB (bp 9001–10200) is replaced by a 1 kb cassette**. For both test samples
(`test_ont`, `test_element`) the summary should show:
- `low_cov_regions` = 1;
- `breakpoints` = 2;
- `cassette_hotspots` = 2;
- `detected_mutations` = 2.

`test/results/` can be deleted afterwards.

---

## 3. Preparing a run

Use one folder per run. Keep the raw FASTQs where they are; the pipeline only links to them.

```bash
mkdir -p ~/runs/WGS_KEIO_2026_11 && cd ~/runs/WGS_KEIO_2026_11
cp ~/Pipeline/template/config.yaml ~/Pipeline/template/sample_sheet.tsv .
```

### 3.1 The sample sheet

The sheet is **tab-separated**, with a header line and one row per sample:

| column | required | meaning |
|---|---|---|
| `sample_name` | yes | becomes file and folder names. No spaces or slashes |
| `platform` | yes | `ont`, `element` or `illumina` |
| `reads1` | yes | **ONT**: a barcode (→ `{fastq_dir}{fastq_base_name}{barcode}.fastq.gz`) or a FASTQ path. **SR**: the R1 FASTQ |
| `reads2` | SR only | the R2 FASTQ (leave empty for ONT) |
| `reference` | see below | GenBank file name in `ref_dir`, with or without `.gb` (or an absolute path). Empty → `default_reference` from the config |

FASTQ paths can be absolute or relative to `fastq_dir`.

Example (`template/sample_sheet.tsv`):

```
sample_name	platform	reads1	reads2	reference
KEIO_ylcJ_ont	ont	barcode05		ecoli_BW25113
KEIO_ylcI_ont	ont	/full/path/to/barcode06.fastq.gz		ecoli_BW25113
ylcJ	element	/full/path/to/ylcJ_R1.fastq.gz	/full/path/to/ylcJ_R2.fastq.gz	ecoli_BW25113
acnB_acnA	element	/full/path/to/acnB_acnA_R1.fastq.gz	/full/path/to/acnB_acnA_R2.fastq.gz	
```

The last row has no reference, so it uses `default_reference`.

- Many editors replace tabs with spaces. Check with `cat -A sample_sheet.tsv`: tabs show as `^I`.
- Lines starting with `#` are ignored.
- **Old ONT sheets still work.** A sheet with the columns `sample_name barcode ref_fasta` is
  read as ONT samples.

### 3.2 The references

Put the GenBank files (`.gb`) in `ref_dir`. They are the reference for alignment **and** the
source of the gene annotation, so they must contain the features (CDS, genes, …).

Replicon names in every output table are the GenBank `VERSION` (e.g. `NC_000913.3`). References
with several records (chromosome + plasmids) are fine.

### 3.3 The config

Edit `config.yaml`. The first part is all you normally touch:

```yaml
out_dir      : "/home/<you>/runs/WGS_KEIO_2026_11/results/"
sample_sheet : "/home/<you>/runs/WGS_KEIO_2026_11/sample_sheet.tsv"
ref_dir      : "/home/<you>/runs/WGS_KEIO_2026_11/genbank_files/"
fastq_dir    : "/net/.../fastq_pass/"       # for ONT barcodes / relative FASTQ paths
fastq_base_name   : ""                       # e.g. "FAX12345_pass_" if files are FAX12345_pass_barcode05.fastq.gz
default_reference : "ecoli_BW25113"          # used where the sheet's reference column is empty

run_read_filtering    : true
run_variant_calling   : true
run_gatc_annotation   : true
run_coverage_analysis : true
run_sv_scan           : true
run_multimer_analysis : true
run_denovo_assembly   : false
run_plasmid_screen    : false
run_cassette_hotspots : true     # only used if run_plasmid_screen is true
run_igv_snapshots     : false
```

Every other parameter is explained in the file itself and in
[section 6](#6-configuration-reference).

---

## 4. Running

### 4.1 Dry run first (login node is fine)

```bash
conda activate seqpipe
~/Pipeline/run.sh ~/runs/WGS_KEIO_2026_11/config.yaml -n
```

A dry run lists the jobs it would run and catches errors in the sheet or config: missing
columns, a FASTQ or reference not found, an unknown platform, and so on. It runs nothing.

### 4.2 The real run, inside salloc

Always start `tmux` **before** `salloc`, so the run survives a dropped connection:

```bash
tmux new -s pipeline
salloc --account=def-rodrigu1 --time=6:00:00 --cpus-per-task=16 --mem=64G
conda activate seqpipe
~/Pipeline/run.sh ~/runs/WGS_KEIO_2026_11/config.yaml
```

- `run.sh` uses **all CPUs and memory of the allocation** and runs as many jobs at once as fit.
  Each step's own threads/memory are in the config (`align_threads`, `medaka_mem_mb`, …).
- **How big an allocation?**
  - 16 CPUs / 64 G is comfortable for about 10–20 bacterial samples.
  - Plasmid-only runs need much less (`--cpus-per-task=8 --mem=16G`).
  - Allow time: a few hours for many ONT samples with medaka and assembly.
- `run.sh` refuses to start a real run outside an allocation. Never run analyses on a login node.
- tmux reminders:
  - `Ctrl-b d` detaches (the run continues);
  - `tmux attach -t pipeline` comes back;
  - `tmux ls` lists sessions.
  - tmux sessions live on one login node, so note which one (`hostname`).

When it finishes it prints `Done. Start with: <out_dir>/summary.tsv`.

### 4.3 Re-running, partial runs

Anything after the config file is passed to snakemake:

```bash
RUN=~/runs/WGS_KEIO_2026_11/config.yaml
~/Pipeline/run.sh $RUN                         # continue / complete a run: only missing outputs are made
~/Pipeline/run.sh $RUN --forcerun lowcov       # redo one step (and everything after it) for all samples
~/Pipeline/run.sh $RUN --config lowcov_min_depth=5    # override one config value for this run
~/Pipeline/run.sh $RUN --until align_sr        # stop after a step
~/Pipeline/run.sh $RUN --unlock                # after a crash: "Directory cannot be locked"
~/Pipeline/run.sh $RUN -n                      # dry run: lists each job and the reason it would run
```

- **To add samples**, add rows to the sheet and run again. Only the new samples are processed,
  plus the run-level tables.
- **To switch an analysis on later**, set its toggle to `true` and run again.
- **After changing a parameter**, Snakemake re-runs the affected steps automatically.

Rule (step) names for `--forcerun` / `--until`:
- **Staging and filtering:** `stage_reads_ont`, `stage_reads_sr`, `filter_reads_ont`, `filter_reads_sr`.
- **Alignment:** `align_ont`, `align_sr`.
- **Variants:** `variant_calling_ont`, `variant_calling_sr`, `filter_variants_ont`, `filter_variants_sr`, `id_gatc_variants`, `annotate_gatc_variants`.
- **Coverage and SV scan:** `lowcov`, `sv_scan`, `annotate_low_coverage`, `annotate_sv`.
- **Other analyses:** `detect_multimers`, `assemble_ont`, `assemble_sr`, `plasmid_screen_align`, `cassette_hotspots`.
- **Summary:** `run_summary`.

---

## 5. Output

```
<out_dir>/
  summary.tsv                         <- start here: one row per sample
  coverage_summary.tsv                <- every sample x replicon (reference + plasmid panel): breadth, depth
  <sample>/
    <sample>.sorted.bam(.bai)         alignment to the sample's reference (open in IGV with <sample>.ref.fasta)
    <sample>.flagstat.txt             % mapped (wrong reference or contamination if low)
    <sample>.ref.gb / .ref.fasta      the reference used
    qc/<sample>.read_stats.tsv        reads / bases / N50, raw and filtered
    qc/<sample>.fastp.html            short reads: trimming report
    <sample>.simple.tsv               SNPs / small indels
    <sample>.annotated.gatc.simple.tsv  ... + NearGATC + gene annotation (+ multimer column, ONT)
    coverage/<sample>.coverage_summary.tsv   per replicon: mean/median depth, % zero, % low
    annotated/<sample>.low_coverage.annotated.tsv
    annotated/<sample>.indels.annotated.tsv
    annotated/<sample>.breakpoints.annotated.tsv
    <sample>.multimer.summary.tsv     ONT: multimer verdict
    <sample>.assembly.fasta / .assembly_info.txt / .assembly_graph.gfa
    <sample>.plasmid_panel.sorted.bam, <sample>.plasmid_coverage.tsv
    plasmidic_reads/<sample>.cassette_hotspots.tsv / .bed, igv_cassette_hotspots/*.png
    logs/<step>.log                   one log per step: look here when something fails
```

### summary.tsv columns

`NA` means the analysis is off or does not apply to the platform. `ERR` means the file was
there but could not be read, and the reason is in `logs/run_summary.log`.

| column | meaning |
|---|---|
| `raw_reads`, `raw_bases`, `filt_reads`, `filt_bases` | before / after filtering (R1+R2 summed for short reads) |
| `filt_N50` | ONT read N50 after filtering |
| `ref_length`, `ref_records` | total reference length (all records), number of records |
| `coverage_approx` | filtered bases / reference length |
| `pct_mapped` | from flagstat |
| `mean_depth` | mean aligned depth over the whole reference |
| `min_replicon_breadth` | lowest % of a replicon covered. **A lost plasmid shows here as ~0** |
| `detected_mutations` | rows in `simple.tsv` |
| `low_cov_regions`, `large_indels`, `breakpoints` | rows in the annotated tables |
| `multimer`, `multimer_copies`, `multimer_pct` | ONT multimer call |
| `assembly_contigs`, `assembly_circular` | number of contigs / circular contigs (`skipped` if too little data) |
| `cassette_hotspots` | number of places where the cassette reads land |

### The region tables (low coverage, indels, breakpoints)

All coordinates are **0-based, half-open** (BED-like). Each table starts with a `#` header line.
The last column, `features`, lists the overlapping GenBank features as `type:name(product)`,
or `.` if there are none.

| table | columns |
|---|---|
| `low_coverage` | contig, start, end, length, mean_depth, min_depth |
| `indels` | contig, start, end, type (`INS`/`DEL`), n_reads, median_len, min_len, max_len, depth, frac_reads |
| `breakpoints` | contig, start, end, side (`L`/`R` soft clip), n_reads, n_split, depth, frac_reads |

How to read them:
- **`frac_reads`** ≈ 1 means the whole population carries the event (clonal), and ≈ 0.1 means a
  minor subpopulation.
- **A real deletion** shows as an `R` breakpoint and an `L` breakpoint on either side of a
  `low_coverage` region with `min_depth` 0.
- **GC-extreme regions** give shallow but non-zero dips.
- **Low-coverage threshold.** A base is "low" when its depth is below
  `max(lowcov_min_depth, lowcov_frac × median depth)`.

`simple.tsv` columns are `CHROM POS REF ALT %Reads MutationType` (1-based `POS`). `%Reads` is the
percentage of reads carrying the ALT allele.

---

## 6. Configuration reference

Everything is in `config.yaml`, which has a comment on each line. The main parameters:

| key | default | meaning |
|---|---|---|
| `stage_mode` | `symlink` | `copy` puts the raw FASTQs inside `out_dir` |
| `min_mapq` | 20 | min mapping quality for low coverage, SV scan, bcftools |
| `filter_mode` | `best95per` | ONT: `best95per`, `best95per_max100x` (also caps at 100× the reference length) or `none` |
| `filtlong_min_length`, `filtlong_keep_percent` | 500, 95 | ONT filtering |
| `fastp_min_length`, `fastp_min_qual` | 30, 20 | short-read trimming |
| `simplify_min_qual` | 50 | ONT variants: keep those in more than this % of reads (despite the name, not a QUAL) |
| `bcftools_min_qual` | 30 | SR variants: minimum QUAL |
| `lowcov_min_depth`, `lowcov_frac` | 10, 0.25 | "low" = depth < max(min_depth, frac × median) |
| `lowcov_min_len`, `lowcov_merge_gap` | 20, 10 | minimum region length; merge regions closer than this |
| `sv_scan.short.*`, `sv_scan.ont.*` | see file | min indel / soft-clip size, min supporting reads, clustering tolerance per read type |
| `sv_ignore_contig_ends` | 10 | ignore breakpoints at the very ends of a replicon (circular origin) |
| `multimer_*` | see file | multimer calling thresholds (ONT) |
| `flye_read_mode` | `--nano-raw` | `--nano-hq` for recent (R10, SUP) basecalls |
| `flye_min_coverage` | 10 | skip assembly below this × reference length (clearly labelled, does not fail) |
| `unicycler_mode` | `normal` | `conservative` / `normal` / `bold` |
| `plasmid_panel_dir` | — | folder with the panel `.fasta` files |
| `cassette_reference` | `Chlor_cassette` | record name of the cassette inside the panel |
| `cassette_hotspot_*` | 1000, 500, 2 | hotspot clustering gap, flank, min reads |
| `igv_exe`, `igv_heap_mb` | — | IGV launcher. Snapshots also need `xvfb-run` on the node |
| `*_threads`, `*_mem_mb` | see file | per-job resources inside the allocation |

---

## 7. Troubleshooting

| symptom | cause / fix |
|---|---|
| `conda: command not found` | load / initialise conda on that cluster (as when you installed) |
| `conda environment 'seqpipe' not found` | do [2.3](#23-create-the-driver-environment-snakemake) |
| run stuck or failing on "Creating conda environment" on a compute node | the environments were not built: run `~/Pipeline/run.sh setup` on a login node ([2.4](#24-build-the-tool-environments-login-node-needs-internet)) |
| `ERROR: not inside a SLURM allocation` | start `salloc` first (inside tmux), or add `-n` for a dry run |
| `sample_sheet is missing required column(s)` | the sheet is not tab-separated, or the header is misspelled |
| `MissingInputException … fastq.gz` | wrong FASTQ path / barcode / `fastq_base_name`. Check one with `ls` |
| `no reference for [...]` | fill the `reference` column or set `default_reference` |
| `Directory cannot be locked` | a previous run crashed: `run.sh <config> --unlock`, then run again |
| a step failed | read `<out_dir>/<sample>/logs/<step>.log`. Fix, then run the same command again |
| allocation time ran out | new `salloc`, same command. Unfinished jobs are redone (`--rerun-incomplete` is on) |
| `features` all `.` in annotated tables | the GenBank file has no features, or they are only of a skipped type (`annotate_skip_types`) |
| very low `pct_mapped` | wrong reference for that sample, or contamination |
| assembly says `# SKIPPED` | not enough data (< `flye_min_coverage` × reference length). Lower it to force a try |
| IGV snapshots fail | check `igv_exe` and that `xvfb-run` exists on the node (`which xvfb-run`) |
| `OUT_OF_MEMORY` / killed jobs | ask more `--mem` in `salloc`, or lower `*_threads` |

---

## 8. Changing the pipeline

```
Pipeline/
  Snakefile          entry point: config, sample sheet, toggles, rule all
  rules/*.smk        one file per analysis (common, qc, align, variants, coverage, sv,
                     multimer, assembly, plasmid_screen, summary)
  envs/*.yaml        conda environment per step (+ snakemake.yaml = driver env)
  scripts/*.py       analysis scripts
  template/          config.yaml + sample_sheet.tsv to copy for a new run
  test/              synthetic test data (make_test_data.py), test config + sheet
  run.sh             launcher
```

- **Read type.** ONT and short reads are told apart with rule-level
  `wildcard_constraints: sample = ONT_RE` or `SR_RE`. Steps that differ have an `_ont` and an
  `_sr` rule that write the same output file.
- **Try a new version of a script.** Copy it as e.g. `scripts/simplify.5.py` and put
  `script: {simplify: "simplify.5.py"}` in a run's config. The keys are the first argument of
  `script_path(...)` in the rules.
- **Add a step.**
  1. Write the rule in a `rules/*.smk` file, with a `conda:` env.
  2. Add its toggle in the `Snakefile` (`RUN_…`).
  3. Add its final output to `_rule_all_inputs()`.
  4. Optionally add a column in `scripts/summary.py`.
- **After any change**, run `./run.sh test` (inside salloc). If an `envs/` file changed, run
  `./run.sh setup` first.

---

## 9. Provenance

Merged on 2026-10-27 from:
- `ONT Analysis/` (NanoForge, top-level version, plus the Flye minimum-coverage check from `V2/`);
- `Element_Analysis/alignqc.py`, whose low-coverage, indel/breakpoint and annotation steps are
  now `scripts/lowcov.py`, `sv_scan.py` and `annotate_regions.py`.

Changes from those pipelines:
- **Read types:** short reads are supported throughout. A new bcftools variant path writes the
  same `simple.tsv` as medaka.
- **Alignment:** ONT alignment is now its own minimap2 step that keeps supplementary
  alignments. medaka's own BAM stays in `<sample>/medaka/`.
- **Coverage:** `coverage.4.py` was dropped and replaced by the low-coverage analysis.
- **Annotation fix:** `id_GATC.2.py` and `annotate.2.py` now use the variant's contig, which
  matters for multi-record references, and fix an off-by-one in feature overlap.
- **Cassette naming:** "Chlor" naming is generalised to `cassette_*`. The old `chlor_*` config
  keys are still read.
- **Assembly bug fixed:** a hard-coded 4.6 Mb minimum that skipped every plasmid assembly.
- **Execution:** no SLURM job submission any more (`salloc` + `run.sh`). Environments are
  shared by all runs.
