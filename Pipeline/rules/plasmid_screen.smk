# ======================================================================
# rules/plasmid_screen.smk :: plasmid-panel presence/coverage screening
# (ONT and short reads)
#
# Every filtered read set is aligned to the combined plasmid panel
# (minimap2 map-ont for ONT, sr for paired short reads); samtools coverage
# gives breadth/depth per panel plasmid -> plasmid_coverage_summary.tsv.
# (Cassette locations in the genome: see rules/cassettes.smk.)
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
    """Align the filtered reads (ONT, or R1+R2) to the combined plasmid panel."""
    input:
        reads = filt_reads,
        panel = OUT_DIR + "plasmid_panel/combined_panel.fasta",
    output:
        bam = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam",
        bai = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam.bai",
    log:
        OUT_DIR + "{sample}/logs/plasmid_screen_align.log",
    threads: config.get("plasmid_screen_threads", 4)
    resources:
        mem_mb = config.get("plasmid_screen_mem_mb", 4000),
    params:
        preset = minimap2_preset,
    conda:
        "../envs/plasmid_screen.yaml"
    shell:
        r"""
        set -euo pipefail

        minimap2 -t {threads} -ax {params.preset} \
            "{input.panel}" {input.reads} 2> "{log}" \
        | samtools sort -@ {threads} -o "{output.bam}" - 2>> "{log}"

        samtools index "{output.bam}" 2>> "{log}"
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
        exe = script_path("coverage_summary", "aggregate_coverage.1.py"),
        samples = SAMPLE_NAMES,
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.exe} -in {input.tsvs} -samples {params.samples} "
        "-name_col replicon -source plasmid_panel -out {output.summary} > {log} 2>&1"
