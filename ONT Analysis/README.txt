NanoForge Chlor_cassette -> reference -> hotspot -> IGV snapshot additions
=========================================================================

Files
-----
Snakefile
    Modified version of the user's supplied main Snakefile. It adds the
    Chlor-cassette targets to rule all and requires the reference alignment
    branch for IGV snapshots.

rules/plasmid_screen.smk
    Complete replacement for the plasmid-screen rule module, preserving the
    original panel/coverage rules and adding:
      - unique Chlor_cassette FASTQ extraction
      - reference alignment
      - FASTA indexing
      - hotspot detection
      - IGV snapshot generation

scripts/extract_chlor_reads.py
    Collects unique read names from Chlor_cassette alignments and extracts
    those reads once from the original filtered FASTQ.

scripts/find_chlor_hotspots.py
    Clusters Chlor-derived read alignments and expands each cluster by 500 bp
    by default.

scripts/run_igv_snapshots.py
    Converts the sample GenBank to GFF3, creates an IGV genome JSON and batch
    file, then launches IGV to produce one PNG per hotspot.

scripts/genbank_to_gff.py
    Standalone GenBank -> GFF3 converter. The IGV script contains the same
    conversion internally, so this file is optional utility/reference code.

config_additions.yaml
    Entries to merge into the existing config.yaml.

plasmid_screen_env_additions.yaml
    Python packages that need to be present in the existing plasmid_screen
    environment: biopython and pysam.

Important
---------
1. The IGV rule assumes the all-read reference BAM is:
       OUT_DIR/{sample}/{sample}.sorted.bam
   which matches the user's variant_calling rule.

2. It assumes the sample reference GenBank is:
       OUT_DIR/{sample}/{sample}.ref.gb
   and FASTA is:
       OUT_DIR/{sample}/{sample}.ref.fasta

3. It assumes the Chlor cassette FASTA record is named:
       Chlor_cassette
   Change chlor_cassette_reference in config.yaml if `samtools idxstats`
   shows a different contig name.

4. The panel alignment intentionally does NOT use --secondary=no. This is
   important when a read has an alignment to Chlor_cassette as well as another
   panel sequence.

5. IGV must be installed separately. Set igv_exe in config.yaml to the path
   of igv.sh if it is not on PATH.

6. On a headless Linux/SLURM node, the script automatically uses xvfb-run
   if that command is available.

Expected output
---------------
OUT_DIR/
  sample/
    sample.sorted.bam
    sample.sorted.bam.bai
    sample.ref.fasta
    sample.ref.gb
    plasmidic_reads/
      Chlor_cassette.fastq.gz
      sample.Chlor_cassette.sorted.bam
      sample.Chlor_cassette.sorted.bam.bai
      sample.chlor_hotspots.tsv
      sample.chlor_hotspots.bed
      igv_chlor_hotspots/
        igv_genome.json
        igv_batch.txt
        <one PNG per hotspot>

Run
---
snakemake --snakefile Snakefile \
          --configfile config.yaml \
          --use-conda \
          --cores N

If you only want to force the new IGV branch after previous outputs exist:
snakemake --snakefile Snakefile \
          --configfile config.yaml \
          --use-conda \
          --cores N \
          --rerun-incomplete
