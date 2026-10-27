# ======================================================================
# rules/plasmid_screen.smk :: plasmid-panel presence/coverage screening
#
# Maps every sample's filtered reads against a shared panel of candidate
# plasmid sequences (config: plasmid_panel_dir -- one or more .fasta/.fa
# files; each fasta record's header becomes that plasmid's name in the
# output). Unlike the single per-sample reference used for variant
# calling, this panel is the same for every sample: it answers "which of
# these plasmids shows up in this sample's reads", not "how does this
# sample compare to its one expected reference".
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
            glob.glob(os.path.join(params.panel_dir, "*.fasta")) +
            glob.glob(os.path.join(params.panel_dir, "*.fa"))
        )
        if not files:
            raise WorkflowError(
                f"no .fasta/.fa files found in plasmid_panel_dir ({params.panel_dir})"
            )
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
        mem_mb  = config.get("plasmid_screen_mem_mb", 4000),
        runtime = config.get("plasmid_screen_runtime", 30),   # minutes (SLURM)
    conda:
        "../envs/plasmid_screen.yaml"
    shell:
        r"""
        set -euo pipefail
        minimap2 -t {threads} -ax map-ont --secondary=yes \
                 "{input.panel}" "{input.reads}" 2> "{log}" \
            | samtools sort -@ {threads} -o "{output.bam}" - 2>> "{log}"
        samtools index "{output.bam}" 2>> "{log}"
        """

rule extract_chlor_cassette_reads:
    """
    Extract reads that align to Chlor_cassette from the plasmid panel BAM.
    """
    input:
        bam = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam"

    output:
        reads = OUT_DIR + "plasmidic_reads/{sample}/Chlor_cassette.fastq.gz"

    log:
        OUT_DIR + "{sample}/logs/extract_chlor_cassette_reads.log"

    conda:
        "../envs/plasmid_screen.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.reads}")"
        samtools view \
            -b \
            "{input.bam}" \
            "Chlor_cassette" \
            2> "{log}" \
        | samtools fastq \
            - \
            2>> "{log}" \
        | gzip > "{output.reads}"
        """

rule align_chlor_cassette_reads_to_reference:
    """
    Align Chlor_cassette-derived reads against the sample reference genome.
    Produces a sorted/indexed BAM suitable for IGV visualization.
    """
    input:
        reads = OUT_DIR + "plasmidic_reads/{sample}/Chlor_cassette.fastq.gz",
        ref   = OUT_DIR + "{sample}/{sample}.ref.fasta"

    output:
        bam = OUT_DIR + "{sample}/plasmidic_reads/{sample}.Chlor_cassette.sorted.bam",
        bai = OUT_DIR + "{sample}/plasmidic_reads/{sample}.Chlor_cassette.sorted.bam.bai"

    log:
        OUT_DIR + "{sample}/logs/align_chlor_cassette_to_reference.log"

    threads: config.get("plasmid_screen_threads", 4)

    resources:
        mem_mb = config.get("plasmid_screen_mem_mb", 4000),
        runtime = config.get("plasmid_screen_runtime", 30)

    conda:
        "../envs/plasmid_screen.yaml"

    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.bam}")"

        minimap2 \
            -t {threads} \
            -ax map-ont \
            "{input.ref}" \
            "{input.reads}" \
            2> "{log}" \
        | samtools sort \
            -@ {threads} \
            -o "{output.bam}" \
            - 2>> "{log}"

        samtools index \
            "{output.bam}" \
            2>> "{log}"
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
        tsvs = expand(OUT_DIR + "{sample}/{sample}.plasmid_coverage.tsv", sample=SAMPLE_NAMES),
    output:
        summary = OUT_DIR + "plasmid_coverage_summary.tsv",
    log:
        OUT_DIR + "logs/aggregate_plasmid_coverage.log",
    params:
        exe     = script("coverage_summary", "aggregate_coverage.1.py"),
        samples = SAMPLE_NAMES,
    shell:
        "python {params.exe} -in {input.tsvs} -samples {params.samples} "
        "-name_col replicon -source plasmid_panel -out {output.summary} > {log} 2>&1"
