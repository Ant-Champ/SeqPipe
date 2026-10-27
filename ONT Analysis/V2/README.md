# NanoForge

ONT analysis for plasmid / small-genome samples: reference-based variant
calling and/or de novo assembly.

## Layout
```
NanoForgeFolder/
├── Snakefile            # entry point (rule all + integration rules)
├── rules/               # step specific rules (QC, alignment, coverage,
│                         #   multimer, assembly, summary)
├── envs/                # per-rule pinned conda environments
└── scripts/             # .py scripts
CurrentRunFolder/
├── config.yaml
├── sample_sheet.tsv
├── fastq/
└── genbank_files/
```

## Installation

The only thing you install globally is Snakemake itself (plus conda).
Every bioinformatics tool is installed automatically per-rule from `envs/*.yaml`
and stored in a **shared environment prefix** so that it is built only once and
then reused by every analysis (see [Environment prefix](#environment-prefix)).

### 1. Install conda
If you don't already have it, install Miniconda:
```bash
# Linux x86_64
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### 2. Create a Snakemake environment
```bash
conda create -n snakemake -c conda-forge -c bioconda \
    python=3.12 'snakemake>=7' biopython
conda activate snakemake
```

### 3. Verify
```bash
snakemake --version          # should print 7.x or newer
conda --version
```

## Environment prefix

By default Snakemake creates the per-rule tool environments under
`.snakemake/conda/` **inside whatever folder you launch the run from**. That
means every new analysis re-solves and re-downloads the entire tool stack and
keeps its own copy

Instead, you can build the environments **once** into a shared location and point every
run at it with `--conda-prefix`. Snakemake identifies each environment by a hash
of its `envs/*.yaml` file, so as long as the pipeline's env files are unchanged,
the existing environments are reused instead of rebuilt.

Pick a stable, absolute path for the shared prefix. Storing it next to the
pipeline install keeps the environments together with the `envs/*.yaml` files
they are built from:
```bash
snake_dir="/home/def-rodrigu1/antoine9/ont/NanoForgeFolder"
conda_prefix="${snake_dir}conda_envs"
```

### One-time setup: build all environments

Run this once after installation (or whenever an `envs/*.yaml` file changes).
It solves and downloads every tool environment without running the workflow:
```bash

cd /net/nfs-bio/home/def-rodrigu1/.../your_run_directory # must use the /net/nfs-bio/home/ path, not just /home/...

snakemake --snakefile ${snake_dir}/Snakefile \
          --configfile config.yaml \
          --use-conda \
          --conda-prefix ${conda_prefix} \
          --conda-create-envs-only \
          --cores 1
```
This step is slow — it solves and downloads everything — but you only pay for it
once. Every subsequent run (below) reuses the environments in `${conda_prefix}`.


## Run (install envs again)

```bash
cd /net/nfs-bio/home/def-rodrigu1/antoine9/ont/Runs/WGS_KEIO_CHLOR  
newgrp def-rodrigu1 # needed for the module/conda to work properly with salloc
salloc --time=2:00:00 -p bio --cpus-per-task=1 --mem=2gb

snake_dir="/home/def-rodrigu1/antoine9/ont/NanoForgeFolder"
conda_prefix="${snake_dir}conda_envs"

conda activate snakemake

snakemake --snakefile ${snake_dir}Snakefile \
          --configfile config.yaml \
          --use-conda \
          --conda-prefix ${conda_prefix} \
          --cores 8
```

## Toggling analyses

Three analysis branches run off the same QC-filtered reads and are
switched on/off for the whole run (not per sample) in `config.yaml`:

* `run_alignment_analysis` -- medaka variant calling against each sample's
  reference. Three sub-analyses inside this branch have their own flags,
  each only meaningful when `run_alignment_analysis` is true:
  * `run_coverage_analysis` -- samtools depth + low-coverage detection,
    plus a `samtools coverage` breadth/depth summary for every replicon
    in the sample's own reference
  * `run_multimer_analysis` -- multimer / concatemer detection
  * `run_gatc_annotation`   -- GATC-site ID + gene annotation (produces
    the final `<sample>.annotated.gatc.simple.tsv` table)
* `run_denovo_assembly` -- Flye de novo assembly, no reference needed.
  Samples containing more than one replicon (e.g. a chromosome + a
  plasmid) are expected to assemble into several contigs; check
  `<sample>.assembly_info.txt` for the per-contig length/coverage/
  circularity breakdown. Samples with too little filtered-read depth for
  their reference size are skipped cleanly (empty fasta, a `# SKIPPED`
  line in `assembly_info.txt`) instead of crashing the run -- tune the
  bar with `flye_min_coverage`.
* `run_plasmid_screen` -- maps every sample against a shared panel of
  candidate plasmids (`plasmid_panel_dir`, one or more fasta files) to
  check which ones are present, independent of each sample's own
  reference.

Enable any combination -- at least one of `run_alignment_analysis` /
`run_denovo_assembly` / `run_plasmid_screen` must stay `true`.

### Coverage summary

Whenever `run_coverage_analysis` and/or `run_plasmid_screen` is on, every
sample's reference replicon(s) and/or screened plasmids land in one
run-level table, **`coverage_summary.tsv`** -- long format, one row per
sample x replicon:

```
sample  source          replicon             ...  coverage  meandepth  ...
S1      reference       BW25113_chromosome   ...  99.96     42.7
S1      plasmid_panel   pCP20                ...  99.44     15.2
S1      plasmid_panel   pKD46                ...  0.63      0.01
```

`source` is `reference` (from the sample's own `ref_fasta`) or
`plasmid_panel` (from the shared panel), so you can tell at a glance
which mapping produced a row -- high `coverage`/`meandepth` means
present, near-zero means absent. The two contributing tables
(`reference_coverage_summary.tsv`, `plasmid_coverage_summary.tsv`) are
still written as intermediates if you want to inspect one source alone.

Note: even with `run_alignment_analysis: false`, every sample still needs
a working reference (`ref_fasta` + matching GenBank file), because the
shared QC/filtering step stages and reads it regardless of the toggles.