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

# 5. results: copy <out_dir>/report.html to your PC and open it in a browser
#    (or read <out_dir>/summary.tsv on the cluster)
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
                ├─ plasmid panel screen                                              ├─ summary.tsv
                └─ genome + all cassettes ─ junctions ─ insertion sites ─ IGV ───────┴─ report.html
```

| Step | config toggle | ONT | Element / Illumina (SR) |
|---|---|---|---|
| Read filtering | `run_read_filtering` | filtlong (length/quality) | fastp (adapters, 3' quality) |
| Alignment to the sample's reference | always | minimap2 `map-ont` | bwa-mem2 |
| Read statistics, flagstat, coverage per replicon | always | ✓ | ✓ |
| SNPs / small indels → `simple.tsv` | `run_variant_calling` | medaka | bcftools (haploid) |
| GATC site + gene annotation of those variants | `run_gatc_annotation` | ✓ | ✓ |
| Depth check: below `min_theoretical_depth` (5×) the sample is "low diversity" and skipped | always (0 = off) | ✓ | ✓ |
| Low-coverage regions (deleted genes), annotated | `run_coverage_analysis` | ✓ | ✓ (repeats masked: `run_mappability_filter`) |
| Large indels + breakpoints, annotated | `run_sv_scan` | ✓ | ✓ |
| Multimer / concatemer detection | `run_multimer_analysis` | ✓ | skipped (reads too short) |
| De novo assembly | `run_denovo_assembly` | Flye | Unicycler |
| Coverage of a shared plasmid panel | `run_plasmid_screen` | ✓ | ✓ |
| Where the cassettes of `cassette_dir` are inserted | `run_cassette_search` | ✓ | ✓ (mates resolve identity) |
| IGV pictures of each insertion site | `run_igv_snapshots` | ✓ | ✓ |
| `summary.tsv`, one row per sample | always | ✓ | ✓ |
| `report.html`: every result, by sample or by analysis, with coverage plots and the IGV pictures | `run_html_report` (default on) | ✓ | ✓ |

**How a deletion is found.** It shows up in three tables:
- in `low_coverage`, as a region with `min_depth` 0 over the deleted gene;
- in `breakpoints`, as an `R` breakpoint at its left edge and an `L` breakpoint at its right edge (soft-clipped reads, with `n_split` reads that align on both sides);
- in `indels`, as a `DEL`, when reads span it (long reads).

**How a cassette insertion is found** (`run_cassette_search`, see [3.4](#34-cassettes)). The
reads are mapped once against the sample's genome plus **all** cassettes. Reads that are
partly genome and partly cassette mark each junction to the base. A left and a right junction
together make one insertion site. For a Keio-style replacement, the bases between the two
junctions are the replaced gene.

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

The test data is synthetic: a 20 kb chromosome plus a 4 kb plasmid. The strain has:
- **2 SNPs in geneA**;
- **geneB (bp 9001–10200) replaced by `TEST_Kan_Cassette`**;
- **geneC (bp 15001–16800) replaced by `TEST_Chlor_Cassette`**, in reverse orientation.

The two test cassettes share identical 150 bp ends, like real Keio cassettes. For both test
samples (`test_ont`, `test_element`) the summary should show:
- `detected_mutations` = 2;
- `low_cov_regions` = 2;
- `breakpoints` = 4;
- `cassette_insertions` = `TEST_Kan_Cassette@testchr.1:9002-10199[geneB];TEST_Chlor_Cassette@testchr.1:15000-16800[geneC]`.

The Kan site is 1–2 bp inside the true ends (9000–10200) because the random test sequence
happens to share a couple of bases with the cassette at those junctions (microhomology). An
aligner cannot place a junction more precisely than that, and real data behaves the same way.

`test/results/` can be deleted afterwards.

---

## 3. Preparing a run

Use one folder per run. Keep the raw FASTQs where they are; the pipeline only links to them.

```bash
mkdir -p ~/runs/WGS_KEIO_2026_11 && cd ~/runs/WGS_KEIO_2026_11
cp ~/Pipeline/template/config.yaml ~/Pipeline/template/sample_sheet.tsv .
```

### 3.1 The sample sheet

The sheet is **tab-separated**, with a header line and one row per sample. All five columns
must be present:

| column | meaning |
|---|---|
| `sample_name` | becomes file and folder names. No spaces or slashes |
| `platform` | `ont`, `element` or `illumina` |
| `reads1` | **full path** of the FASTQ: the ONT reads, or R1 for short reads (`.fastq`/`.fq`, optionally `.gz`) |
| `reads2` | **full path** of R2 for short reads. **Leave empty for ONT** |
| `reference` | GenBank file name in `ref_dir`, with or without `.gb`, or a full path to a `.gb` file |

Example (`template/sample_sheet.tsv`):

```
sample_name	platform	reads1	reads2	reference
KEIO_ylcJ_ont	ont	/full/path/to/fastq_pass/barcode05.fastq.gz		ecoli_BW25113
KEIO_ylcI_ont	ont	/full/path/to/fastq_pass/barcode06.fastq.gz		ecoli_BW25113
ylcJ	element	/full/path/to/ylcJ_R1.fastq.gz	/full/path/to/ylcJ_R2.fastq.gz	ecoli_BW25113
acnB_acnA	element	/full/path/to/acnB_acnA_R1.fastq.gz	/full/path/to/acnB_acnA_R2.fastq.gz	ecoli_BW25113
```

**Nothing has a default.** Before anything runs, even a dry run (`-n`), the pipeline checks
every row and **stops with the full list of problems**. It stops on:
- an empty or unknown `platform`;
- an empty `reads1`, or an empty `reads2` for a short-read sample;
- a path that is not a FASTQ file, or a FASTQ file that doesn't exist;
- a `reads2` value on an ONT sample;
- an empty `reference`, or a reference file not found;
- a duplicated `sample_name`.

Use full paths (`/home/...`, `/net/...`). A relative path is read from the run folder (the
folder that contains `config.yaml`).

To fill the ONT column quickly, list the files with full paths:
```bash
ls -1 /full/path/to/fastq_pass/*.fastq.gz
```

- Many editors replace tabs with spaces. Check with `cat -A sample_sheet.tsv`: tabs show as `^I`.
- Lines starting with `#` are ignored.
- Old ONT sheets with `barcode` / `ref_fasta` columns are no longer accepted. The pipeline stops
  and says so.

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

run_read_filtering    : true
run_variant_calling   : true
run_gatc_annotation   : true
run_coverage_analysis : true
run_sv_scan           : true
run_multimer_analysis : true
run_denovo_assembly   : false
run_plasmid_screen    : false
run_cassette_search   : false    # true: locate every cassette of cassette_dir, in every sample
run_igv_snapshots     : false
cassette_dir          : "/home/antoine9/Pipeline/references/cassettes/"
```

Every other parameter is explained in the file itself and in
[section 6](#6-configuration-reference).

### 3.4 Cassettes

To find where antibiotic cassettes are inserted:

1. Put **one file per cassette** in a folder: GenBank (`.gb`, `.gbk`) or FASTA (`.fasta`, `.fa`).
   The file name, without its extension, becomes the cassette name in every table:
   ```
   /home/antoine9/Pipeline/references/cassettes/
       KEIO_Kan_Cassette.gb
       KEIO_Chlor_sfGFP_Cassette.gb
   ```
2. In `config.yaml`, set `run_cassette_search : true` and point `cassette_dir` at that folder.

When the search is on, **every sample is searched for every cassette in the folder**. To add a
cassette, drop its file in the folder. To leave one out, move its file elsewhere.

**When cassettes share sequence.** Keio cassettes have identical FRT / priming-site ends, and
those ends are exactly where the junction reads fall. The pipeline handles this as follows:
- **One combined search.** All cassettes are searched together, so a read that fits two
  cassettes equally well gets mapping quality 0 instead of being given to one of them at random.
- **The position is always reported.** It is reliable even when the read only covers the shared
  end.
- **The cassette's name needs unique evidence.** It is only given from evidence with mapping
  quality ≥ `cassette_min_mapq`: the read itself (ONT reads usually reach the unique middle of the
  cassette), or for short reads its mate.
- **Otherwise the site is labelled honestly:** `ambiguous(A|B)` when nothing unique is seen, or
  `mixed(A:n,B:m)` when both cassettes have unique support at the same site. `mixed` can mean a
  mixed population, or two cassettes sitting next to each other.
- **Read counts per cassette** (`cassette_counts.tsv`) count uniquely mapped reads only, so a
  shared end cannot make an absent cassette look present.
  - These counts are primary alignments only. With ONT, most reads that cross a cassette
    have their primary alignment on the genome, so the counts stay small (a few reads) even
    for a cassette that is clearly there.
  - For ONT, judge presence from the insertion sites, not from these counts.

The search does not depend on the plasmid panel screen (`run_plasmid_screen`). The two can be
switched on independently.

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
- **Coverage and SV scan:** `lowcov` (low and high coverage), `sv_scan`, `annotate_coverage_regions`, `annotate_sv`.
- **Other analyses:** `detect_multimers`, `assemble_ont`, `assemble_sr`, `plasmid_screen_align`, `cassette_align_ont`, `cassette_align_sr`, `find_cassette_insertions`.
- **Summary and report:** `run_summary`, `depth_profile`, `html_report`.

---

## 5. Output

```
<out_dir>/
  report.html                         <- start here: everything in one page (see 5.1)
  summary.tsv                         <- one row per sample
  coverage_summary.tsv                <- every sample x replicon (reference + plasmid panel): breadth, depth
  depth_check.tsv                     <- theoretical depth of every sample, ok / low_diversity
  mappability/<ref>.pe75_ins300_q20.unmappable.bed   short reads: repeats ignored by the low-coverage scan
  <sample>/
    <sample>.sorted.bam(.bai)         alignment to the sample's reference (open in IGV with <sample>.ref.fasta)
    <sample>.flagstat.txt             % mapped (wrong reference or contamination if low)
    <sample>.ref.gb / .ref.fasta      the reference used
    qc/<sample>.read_stats.tsv        reads / bases / N50, raw and filtered
    qc/<sample>.fastp.html            short reads: trimming report
    qc/<sample>.depth_check.tsv       theoretical depth (a low_diversity sample stops here)
    <sample>.simple.tsv               SNPs / small indels
    <sample>.annotated.gatc.simple.tsv  ... + NearGATC + gene annotation (+ multimer column, ONT)
    coverage/<sample>.coverage_summary.tsv   per replicon: mean/median depth, % zero, % low, % masked
    annotated/<sample>.low_coverage.annotated.tsv    deletions (with the genes hit)
    annotated/<sample>.high_coverage.annotated.tsv   duplications / amplifications (fold vs replicon median)
    annotated/<sample>.indels.annotated.tsv
    annotated/<sample>.breakpoints.annotated.tsv
    <sample>.multimer.summary.tsv     ONT: multimer verdict
    <sample>.assembly.fasta / .assembly_info.txt / .assembly_graph.gfa
    <sample>.plasmid_panel.sorted.bam, <sample>.plasmid_coverage.tsv
    annotated/<sample>.cassette_insertions.annotated.tsv   insertion sites (+ replaced genes)
    annotated/<sample>.cassette_junctions.annotated.tsv    every genome/cassette junction
    cassettes/<sample>.cassette_counts.tsv                 uniquely mapped reads per cassette
    cassettes/<sample>.genome_cassettes.sorted.bam         reads vs genome + cassettes (IGV), igv_snapshots/*.png
    coverage/<sample>.depth_bins.tsv  binned depth along each replicon (for the report plots)
    logs/<step>.log                   one log per step: look here when something fails
```

### 5.1 The HTML report

`report.html` is one file: copy it to your PC and open it in any browser. It needs no internet
connection. To copy it, use WinSCP, or from WSL:

```bash
scp narval:~/runs/MY_RUN/results/report.html /mnt/c/Users/Ant/Downloads/
```

The two buttons at the top switch between two views:
- **By sample:** pick a sample on the left. Its page shows key numbers, then every analysis for
  that sample: read QC, coverage plots, variants, low coverage, indels and breakpoints,
  cassettes with their IGV pictures, plasmid panel, multimers and assembly.
- **By analysis:** pick an analysis on the left. Its page shows that analysis for every sample,
  in one table with a `sample` column. *Overview* holds `summary.tsv` and the run settings.
  *Plasmid panel* shows a sample × plasmid presence matrix.

The coverage plots show depth along each replicon, with coloured marks:

| Colour | Meaning |
|---|---|
| red band | low-coverage region |
| teal band | high-coverage region (duplication) |
| green | cassette insertion |
| orange line | breakpoint |
| purple | large indel |
| black tick | variant |
| grey | genes (shown below 300 kb of view) |

To use the plots:
- Drag across a plot to zoom; double-click it to reset.
- Hover over a plot to see the depth, gene and event at that position.
- **Click a row of any table that has a position** (variant, low-coverage region, breakpoint,
  cassette site) to jump to that place on the sample's plot.

Tables can be sorted (click a column header) and filtered (type in the box). A table keeps its
first 5000 rows; the TSV on disk is complete.

The IGV pictures appear when `run_igv_snapshots` is on. They are embedded in the file;
set `report_embed_igv: false` to link to the PNG files instead and keep `report.html` small. The
links to the fastp reports and the assemblies are relative, so they only work if `report.html`
stays in `<out_dir>`, or if you copy the whole `<out_dir>`.

### summary.tsv columns

`NA` means the analysis is off or does not apply to the platform. `ERR` means the file was
there but could not be read, and the reason is in `logs/run_summary.log`.

| column | meaning |
|---|---|
| `status` | `ok`, or `low_diversity (1.97x < 5x)`: too few bases for `min_theoretical_depth`, the sample was not aligned or analysed (other columns `NA`) |
| `raw_reads`, `raw_bases`, `filt_reads`, `filt_bases` | before / after filtering (R1+R2 summed for short reads) |
| `filt_N50` | ONT read N50 after filtering |
| `ref_length`, `ref_records` | total reference length (all records), number of records |
| `coverage_approx` | filtered bases / reference length |
| `pct_mapped` | from flagstat |
| `mean_depth` | mean aligned depth over the whole reference |
| `min_replicon_breadth` | lowest % of a replicon covered. **A lost plasmid shows here as ~0** |
| `detected_mutations` | rows in `simple.tsv` |
| `low_cov_regions`, `high_cov_regions`, `large_indels`, `breakpoints` | rows in the annotated tables |
| `multimer`, `multimer_copies`, `multimer_pct` | ONT multimer call |
| `assembly_contigs`, `assembly_circular` | number of contigs / circular contigs (`skipped` if too little data) |
| `cassette_insertions` | `cassette@contig:start-end[replaced genes]` per site, `;`-separated, `none` if no site. One-sided sites read `cassette@contig:pos(left_only)` |

### The region tables (low coverage, indels, breakpoints)

All coordinates are **0-based, half-open** (BED-like). Each table starts with a `#` header line.
The last column, `features`, lists the overlapping GenBank features as `type:name(product)`,
or `.` if there are none.

| table | columns |
|---|---|
| `low_coverage` | contig, start, end, length, mean_depth, min_depth, masked_bp, rel_depth (mean / median) |
| `high_coverage` | contig, start, end, length, mean_depth, replicon_median, fold |
| `indels` | contig, start, end, type (`INS`/`DEL`), n_reads, median_len, min_len, max_len, depth, frac_reads |
| `breakpoints` | contig, start, end, side (`L`/`R` soft clip), n_reads, n_split, depth, frac_reads |

How to read them:
- **`frac_reads`** ≈ 1 means the whole population carries the event (clonal), and ≈ 0.1 means a
  minor subpopulation.
- **A real deletion** shows as an `R` breakpoint and an `L` breakpoint on either side of a
  `low_coverage` region with `min_depth` 0.
- **GC-extreme regions** give shallow but non-zero dips.
- **Low-coverage threshold.** It is relative, with a fixed floor: a base is "low" when its
  depth, averaged over `lowcov_smooth` bp, is below `max(lowcov_min_depth, lowcov_frac × median depth)`,
  the median being taken over the whole reference. With the defaults: median 80× → low below 20×;
  median 30× → low below 10× (the floor). The threshold is capped at `lowcov_max_thr_frac × median`
  (0.5): at a median of 8×, a 5× floor would otherwise call ordinary depth dips "low" all over
  the genome. Low stretches closer than `lowcov_merge_gap` are merged; a region is kept only if it
  is at least `lowcov_min_len` bp **and** its mean depth is at most `lowcov_region_max_frac × median`
  (0.25): a real deletion is ~0× all along, a chain of dips is not. The log
  (`<sample>/logs/lowcov.log`) prints the threshold used and how many candidates were dropped.
- **Mappability mask (short reads).** Reads from repeats (rRNA operons, IS elements, tRNA
  clusters) get MAPQ 0 and look like a deletion. With `run_mappability_filter` (default on), perfect
  PE reads (`mappability_read_len`, `mappability_insert`) are simulated from each reference,
  mapped back with bwa-mem2, and the bases they cannot cover at MAPQ ≥ `min_mapq` are ignored by
  the low-coverage scan (and left out of the median). The mask is
  `<out_dir>/mappability/<reference>.pe75_ins300_q20.unmappable.bed` (once per reference; set
  `mappability_dir` to a shared folder to reuse it between runs). `masked_bp` in the table says how
  much of a region is masked; `pct_masked` is in `coverage/<sample>.coverage_summary.tsv`.
- **A duplication** shows as a `high_coverage` region with `fold` ≈ 2 (≈ 3 for two extra
  copies), often with breakpoints at its ends. The threshold is relative to the median of the
  replicon the region is on, so a multicopy plasmid is not flagged as a whole (its copy number
  is in `coverage_summary.tsv`). Short reads give noisier depth than ONT; raise
  `highcov_min_len` if small GC-rich regions come up.

`simple.tsv` columns are `CHROM POS REF ALT %Reads MutationType` (1-based `POS`). `%Reads` is the
percentage of reads carrying the ALT allele.

### The cassette tables

Coordinates are on the genome, 0-based like the other region tables. Each table ends with the
`features` column.

| table | columns |
|---|---|
| `cassette_insertions` | contig, start, end, length, status, cassette, orientation, n_reads_left, n_reads_right, left_junction, right_junction |
| `cassette_junctions` | contig, start, end, genome_side, cassette, cassette_end, orientation, n_reads, n_unique, n_ambiguous, support, junction_id |
| `cassette_counts` | cassette, length, reads_unique, reads_ambiguous, breadth_unique_pct, mean_depth_unique |

The columns, in plain words:
- **`status`**:
  - `paired` means a junction was found on both sides. `start`–`end` is then the replaced
    stretch (`length` bp), and `features` lists the replaced genes.
  - `left_only` / `right_only` means only one junction was found, for example at low depth or
    when one side lies in a repeat.
- **`cassette`**: the cassette's name, `ambiguous(A|B)` (only reads in the shared ends), or
  `mixed(A:n,B:m)`.
- **`orientation`**: `+` means the cassette reads in the same direction as the genome.
- **`support`**: uniquely assigned reads per cassette at that junction.
- **`genome_side`**: `left` means the genome is on the left of the junction and the cassette
  continues to the right.

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
| `min_theoretical_depth` | 5 | filtered bases / reference length below this → `low_diversity`: sample skipped after QC, red banner in the report. 0 = off |
| `lowcov_min_depth`, `lowcov_frac` | 10, 0.25 | "low" = smoothed depth < max(min_depth, frac × median) |
| `lowcov_max_thr_frac` | 0.5 | ...but the threshold never exceeds this × the median |
| `lowcov_smooth` | 50 | bp over which depth is averaged before calling "low" (1 = per base) |
| `lowcov_min_len`, `lowcov_merge_gap` | 20, 10 | minimum region length (raise it, e.g. 500, to keep only large deletions); merge regions closer than this |
| `lowcov_region_max_frac` | 0.25 | drop a region whose mean depth is above this × the median |
| `run_mappability_filter` | true | short reads: ignore unmappable (repeat) bases in the low-coverage scan |
| `mappability_read_len`, `mappability_insert` | 75, 300 | read length and fragment size of the simulated pairs |
| `mappability_step`, `mappability_min_frac` | 2, 0.5 | one pair every N bp; unmappable below this fraction of the expected depth |
| `mappability_dir` | `<out_dir>/mappability/` | where the masks go (a shared folder reuses them between runs) |
| `highcov_frac`, `highcov_window` | 1.7, 100 | "high" = a window whose mean depth > frac × the median of its own replicon (duplications; fold ≈ 2 per extra copy) |
| `highcov_min_len`, `highcov_merge_gap` | 1000, 300 | minimum high-coverage region length; merge regions closer than this |
| `sv_scan.short.*`, `sv_scan.ont.*` | see file | min indel / soft-clip size, min supporting reads, clustering tolerance per read type |
| `sv_ignore_contig_ends` | 10 | ignore breakpoints at the very ends of a replicon (circular origin) |
| `multimer_*` | see file | multimer calling thresholds (ONT) |
| `flye_read_mode` | `--nano-raw` | `--nano-hq` for recent (R10, SUP) basecalls |
| `flye_min_coverage` | 10 | skip assembly below this × reference length (clearly labelled, does not fail) |
| `unicycler_mode` | `normal` | `conservative` / `normal` / `bold` |
| `plasmid_panel_dir` | — | folder with the panel plasmids: `.gb/.gbk/.genbank/.fasta/.fa/.fna`, one plasmid per file (named after the file) or several records per file (`file_recordid`). Use an absolute path: a relative one is read from the config's folder. Checked at start-up |
| `cassette_dir` | — | folder with one `.gb`/`.fasta` per cassette. All of them are searched |
| `cassette_min_mapq` | 20 | evidence below this mapping quality cannot name the cassette |
| `cassette_min_reads` | 2 | minimum reads for a junction |
| `cassette_max_span` | 20000 | max distance between the two junctions of one site |
| `cassette_min_segment_*`, `cassette_cluster_tol_*` | 100/30, 20/10 | min aligned bp on each side of a junction; position tolerance (ONT / short) |
| `igv_exe`, `igv_heap_mb`, `igv_flank` | —, 4000, 1000 | IGV launcher, Java heap, bp shown around a site. Snapshots also need `xvfb-run` on the node |
| `run_html_report` | true | write `<out_dir>/report.html` (see 5.1) |
| `report_max_bins` | 20000 | depth points per replicon in the report plots (a 4.6 Mb genome → 230 bp per point) |
| `report_embed_igv` | true | put the IGV pictures inside `report.html` (false: link to the PNG files) |
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
| `N problem(s) in sample_sheet …` | the list below it says which sample and column. Fix them all, then re-run |
| `config key(s) fastq_dir, … are no longer used` | an old config: delete those lines and put full FASTQ paths and the reference in the sheet |
| `old ONT sample sheet (barcode / ref_fasta columns)` | rewrite the sheet with the five columns ([3.1](#31-the-sample-sheet)) |
| `Directory cannot be locked` | a previous run crashed: `run.sh <config> --unlock`, then run again |
| a step failed | read `<out_dir>/<sample>/logs/<step>.log`. Fix, then run the same command again |
| allocation time ran out | new `salloc`, same command. Unfinished jobs are redone (`--rerun-incomplete` is on) |
| `features` all `.` in annotated tables | the GenBank file has no features, or they are only of a skipped type (`annotate_skip_types`) |
| very low `pct_mapped` | wrong reference for that sample, or contamination |
| assembly says `# SKIPPED` | not enough data (< `flye_min_coverage` × reference length). Lower it to force a try |
| IGV snapshots fail | check `igv_exe` and that `xvfb-run` exists on the node (`which xvfb-run`) |
| no IGV pictures, no error | `run_igv_snapshots` is `false` (the report says "IGV snapshots are off"). Set it to `true` (and `run_cassette_search: true`) |
| a sample has a red "Low diversity" banner | too few bases for `min_theoretical_depth`: re-sequence, or lower `min_theoretical_depth` to force the analysis |
| low-coverage regions all over a low-depth sample | check `<sample>/logs/lowcov.log` (threshold, median). Lower `lowcov_min_depth`, raise `lowcov_min_len` (e.g. 500), or lower `lowcov_region_max_frac` |
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
- **Cassette search (2026-10-27):** the old method (pull out the reads that hit a cassette in the
  plasmid panel, re-map them, cluster them into "hotspots") is replaced by `run_cassette_search`.
  The new search uses junction reads against genome + all cassettes, and handles homology
  between cassettes. The keys `run_cassette_hotspots`, `cassette_reference`,
  `cassette_hotspot_*` and `chlor_*` are no longer used.
- **Assembly bug fixed:** a hard-coded 4.6 Mb minimum that skipped every plasmid assembly.
- **Execution:** no SLURM job submission any more (`salloc` + `run.sh`). Environments are
  shared by all runs.
