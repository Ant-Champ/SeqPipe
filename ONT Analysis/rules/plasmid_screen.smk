# ======================================================================
# rules/plasmid_screen.smk :: plasmid-panel presence/coverage screening
# + Chlor_cassette extraction, reference alignment, hotspot detection,
#   and automated IGV snapshots
# ======================================================================

rule build_plasmid_panel:
    """Concatenate every fasta in plasmid_panel_dir into one combined reference."""
    output:
        panel = OUT_DIR + "plasmid_panel/combined_panel.fasta",
    log:
        OUT_DIR + "logs/build_plasmid_panel.log",
    params:
        panel_dir = PLASMID_PANEL_DIR,
    run:
        import glob
        import os

        files = sorted(
            glob.glob(os.path.join(params.panel_dir, "*.fasta"))
            + glob.glob(os.path.join(params.panel_dir, "*.fa"))
        )

        if not files:
            raise WorkflowError(
                f"no .fasta/.fa files found in plasmid_panel_dir ({params.panel_dir})"
            )

        os.makedirs(os.path.dirname(output.panel), exist_ok=True)

        with open(output.panel, "w") as out_fh, open(log[0], "w") as log_fh:
            log_fh.write("plasmid panel built from:\n")
            for f in files:
                log_fh.write(f"  {f}\n")
                with open(f) as in_fh:
                    out_fh.write(in_fh.read())


rule plasmid_screen_align:
    input:
        reads = OUT_DIR + "{sample}/{sample}.filt.fastq.gz",
        panel = OUT_DIR + "plasmid_panel/combined_panel.fasta",
    output:
        bam = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam",
        bai = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam.bai",
    log:
        OUT_DIR + "{sample}/logs/plasmid_screen_align.log",
    threads: config.get("plasmid_screen_threads", 4)
    resources:
        mem_mb = config.get("plasmid_screen_mem_mb", 4000),
        runtime = config.get("plasmid_screen_runtime", 30),
    conda:
        "../envs/plasmid_screen.yaml"
    shell:
        r"""
        set -euo pipefail

        minimap2 -t {threads} -ax map-ont \
            "{input.panel}" "{input.reads}" 2> "{log}" \
        | samtools sort -@ {threads} -o "{output.bam}" - 2>> "{log}"

        samtools index "{output.bam}" 2>> "{log}"
        """


rule extract_chlor_cassette_reads:
    """
    Extract unique original FASTQ reads having an alignment to Chlor_cassette.
    Read names are collected from the panel BAM, then the original FASTQ is
    filtered so supplementary/secondary alignments do not duplicate a read.
    """
    input:
        bam = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam",
        reads = OUT_DIR + "{sample}/{sample}.filt.fastq.gz",
    output:
        reads = OUT_DIR + "plasmidic_reads/{sample}/Chlor_cassette.fastq.gz",
    log:
        OUT_DIR + "{sample}/logs/extract_chlor_cassette_reads.log",
    params:
        exe = script("extract_chlor_reads", "extract_chlor_reads.py"),
        chlor_ref = config.get("chlor_cassette_reference", "Chlor_cassette"),
    conda:
        "../envs/plasmid_screen.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.reads}")"

        # First isolate alignments to the requested reference sequence and
        # let the helper script collect unique query names.
        samtools view -b "{input.bam}" "{params.chlor_ref}" \
            2> "{log}" > "{output.reads}.tmp.bam"

        python "{params.exe}" \
            --bam "{output.reads}.tmp.bam" \
            --reads "{input.reads}" \
            --out "{output.reads}" \
            >> "{log}" 2>&1

        rm -f "{output.reads}.tmp.bam"
        """


rule index_reference_for_chlor_igv:
    """Create the FASTA index required by minimap2/IGV for the sample reference."""
    input:
        ref = OUT_DIR + "{sample}/{sample}.ref.fasta",
    output:
        fai = OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
    log:
        OUT_DIR + "{sample}/logs/index_reference_for_chlor_igv.log",
    conda:
        "../envs/plasmid_screen.yaml"
    shell:
        r"""
        set -euo pipefail
        samtools faidx "{input.ref}" > "{log}" 2>&1
        """


rule align_chlor_cassette_reads_to_reference:
    """
    Align Chlor_cassette-derived reads against the sample reference genome.
    The resulting BAM is intended for IGV and hotspot detection.
    """
    input:
        reads = OUT_DIR + "plasmidic_reads/{sample}/Chlor_cassette.fastq.gz",
        ref = OUT_DIR + "{sample}/{sample}.ref.fasta",
        ref_fai = OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
    output:
        bam = OUT_DIR + "{sample}/plasmidic_reads/{sample}.Chlor_cassette.sorted.bam",
        bai = OUT_DIR + "{sample}/plasmidic_reads/{sample}.Chlor_cassette.sorted.bam.bai",
    log:
        OUT_DIR + "{sample}/logs/align_chlor_cassette_to_reference.log",
    threads: config.get("plasmid_screen_threads", 4)
    resources:
        mem_mb = config.get("plasmid_screen_mem_mb", 4000),
        runtime = config.get("plasmid_screen_runtime", 30),
    conda:
        "../envs/plasmid_screen.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.bam}")"

        minimap2 -t {threads} -ax map-ont \
            "{input.ref}" "{input.reads}" 2> "{log}" \
        | samtools sort -@ {threads} -o "{output.bam}" - 2>> "{log}"

        samtools index "{output.bam}" 2>> "{log}"
        """


rule chlor_hotspots:
    """
    Find clusters of Chlor_cassette-read alignments on the reference genome.
    Each cluster is expanded by chlor_hotspot_flank bp for IGV visualization.
    """
    input:
        bam = OUT_DIR + "{sample}/plasmidic_reads/{sample}.Chlor_cassette.sorted.bam",
    output:
        tsv = OUT_DIR + "{sample}/plasmidic_reads/{sample}.chlor_hotspots.tsv",
        bed = OUT_DIR + "{sample}/plasmidic_reads/{sample}.chlor_hotspots.bed",
    log:
        OUT_DIR + "{sample}/logs/chlor_hotspots.log",
    params:
        exe = script("find_chlor_hotspots", "find_chlor_hotspots.py"),
        max_gap = config.get("chlor_hotspot_max_gap", 1000),
        flank = config.get("chlor_hotspot_flank", 500),
        min_reads = config.get("chlor_hotspot_min_reads", 2),
    conda:
        "../envs/plasmid_screen.yaml"
    shell:
        r"""
        set -euo pipefail

        python "{params.exe}" \
            --bam "{input.bam}" \
            --out-tsv "{output.tsv}" \
            --out-bed "{output.bed}" \
            --max-gap {params.max_gap} \
            --flank {params.flank} \
            --min-reads {params.min_reads} \
            > "{log}" 2>&1
        """


rule chlor_igv_snapshots:
    """
    Generate one IGV PNG per Chlor_cassette hotspot.

    The snapshot contains:
      - reference sequence
      - GenBank-derived annotations
      - all-read BAM
      - Chlor_cassette-read BAM
      - hotspot +/- chlor_hotspot_flank bp
    """
    input:
        hotspots = OUT_DIR + "{sample}/plasmidic_reads/{sample}.chlor_hotspots.tsv",
        all_bam = OUT_DIR + "{sample}/{sample}.sorted.bam",
        all_bai = OUT_DIR + "{sample}/{sample}.sorted.bam.bai",
        chlor_bam = OUT_DIR + "{sample}/plasmidic_reads/{sample}.Chlor_cassette.sorted.bam",
        chlor_bai = OUT_DIR + "{sample}/plasmidic_reads/{sample}.Chlor_cassette.sorted.bam.bai",
        ref = OUT_DIR + "{sample}/{sample}.ref.fasta",
        ref_fai = OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
        gb = OUT_DIR + "{sample}/{sample}.ref.gb",
    output:
        snapshots = directory(
            OUT_DIR + "{sample}/plasmidic_reads/igv_chlor_hotspots"
        ),
    log:
        OUT_DIR + "{sample}/logs/chlor_igv_snapshots.log",
    params:
        exe = script("run_igv_snapshots", "run_igv_snapshots.py"),
        igv = config.get("igv_exe", "igv.sh"),
    resources:
        mem_mb = config.get("igv_mem_mb", 6000),
        runtime = config.get("igv_runtime", 30),
    conda:
        "../envs/plasmid_screen.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{output.snapshots}"

        python "{params.exe}" \
            --hotspots "{input.hotspots}" \
            --all-bam "{input.all_bam}" \
            --chlor-bam "{input.chlor_bam}" \
            --reference "{input.ref}" \
            --genbank "{input.gb}" \
            --outdir "{output.snapshots}" \
            --igv "{params.igv}" \
            --memory 4000 \
            > "{log}" 2>&1
        """


rule plasmid_screen_coverage:
    """Per-plasmid breadth/depth for one sample (one row per panel plasmid)."""
    input:
        bam = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam",
    output:
        tsv = OUT_DIR + "{sample}/{sample}.plasmid_coverage.tsv",
    log:
        OUT_DIR + "{sample}/logs/plasmid_screen_coverage.log",
    conda:
        "../envs/plasmid_screen.yaml"
    shell:
        r"""samtools coverage "{input.bam}" > "{output.tsv}" 2> "{log}" """


rule aggregate_plasmid_coverage:
    """Run-level long-format table: one row per sample x plasmid."""
    input:
        tsvs = expand(
            OUT_DIR + "{sample}/{sample}.plasmid_coverage.tsv",
            sample=SAMPLE_NAMES
        ),
    output:
        summary = OUT_DIR + "plasmid_coverage_summary.tsv",
    log:
        OUT_DIR + "logs/aggregate_plasmid_coverage.log",
    params:
        exe = script("coverage_summary", "aggregate_coverage.1.py"),
        samples = SAMPLE_NAMES,
    shell:
        "python {params.exe} -in {input.tsvs} -samples {params.samples} "
        "-name_col replicon -source plasmid_panel -out {output.summary} > {log} 2>&1"
