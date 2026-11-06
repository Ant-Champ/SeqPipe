# SeqPipe

A Snakemake pipeline for whole-genome sequencing of bacteria, plasmids and phages, from
**Oxford Nanopore** long reads and **paired-end short reads** (Element AVITI, Illumina), mapped
against **GenBank references**. ONT and short-read samples can be mixed in one run.

For every sample it reports:

- read QC and filtering (filtlong / fastp);
- SNPs and small indels (medaka / bcftools), annotated with genes and nearby GATC sites;
- low- and high-coverage regions (deletions, duplications), large indels and breakpoints;
- multimers / concatemers (ONT) and de novo assemblies (Flye / Unicycler);
- coverage of a panel of plasmids;
- where antibiotic cassettes (e.g. Keio) are inserted, with a picture of each insertion site.

Everything ends up in one `summary.tsv` and a self-contained `report.html` that opens offline
in any browser, with zoomable coverage plots and the aligned reads.

## Quick start

Requires Linux and [conda](https://conda-forge.org/download/). Every tool is installed by the
pipeline in its own conda environment.

```bash
git clone https://github.com/Ant-Champ/SeqPipe.git
cd SeqPipe/Pipeline
conda env create -f envs/snakemake.yaml     # the driver environment, "seqpipe"
conda activate seqpipe
./run.sh setup                              # build every tool environment (once, needs internet)
SEQPIPE_NO_SLURM=1 ./run.sh test            # run the bundled synthetic test data (a few minutes)
```

On a SLURM cluster, run analyses inside a `salloc` allocation instead of setting
`SEQPIPE_NO_SLURM`.

**Full documentation:** [`Pipeline/README.md`](Pipeline/README.md): installation, the sample
sheet, every option, the output files and troubleshooting.

## Repository layout

```
Pipeline/
  README.md     full documentation
  run.sh        launcher
  Snakefile     workflow entry point
  rules/        one .smk file per analysis
  scripts/      analysis and report scripts
  envs/         conda environment of each step
  template/     config.yaml and sample_sheet.tsv to copy for a new run
  test/         synthetic test data and its config
```

## License

SeqPipe is free software, released under the [GNU Affero General Public License v3.0](LICENSE).
