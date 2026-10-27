# ======================================================================
# rules/coverage.smk :: per-base depth + low-coverage region detection
# ======================================================================

rule compute_coverage:
    input:
        bam = OUT_DIR + "{sample}/{sample}.sorted.bam",
    output:
        cov = OUT_DIR + "{sample}/{sample}.cov.txt",
    log:
        OUT_DIR + "{sample}/logs/compute_coverage.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""samtools depth -a "{input.bam}" > "{output.cov}" 2> "{log}" """


rule analyze_coverage:
    input:
        cov       = OUT_DIR + "{sample}/{sample}.cov.txt",
        ref_fasta = OUT_DIR + "{sample}/{sample}.ref.fasta",
    output:
        log_file = OUT_DIR + "{sample}/{sample}.cov.log",
    params:
        exe        = script("coverage", "coverage.4.py"),
        low_cov_thr= config.get("coverage_low_threshold", 5),
    shell:
        "python {params.exe} -cov {input.cov} -ref {input.ref_fasta} "
        "-l {params.low_cov_thr} > {output.log_file}"


# ----------------------------------------------------------------------
# samtools-coverage style breadth/depth summary for the sample's own
# reference (same metrics as rules/plasmid_screen.smk's panel screen,
# just against ref.fasta instead of the shared plasmid panel). Reuses the
# BAM already produced by variant_calling -- no extra mapping needed.
# ----------------------------------------------------------------------
rule reference_coverage:
    input:
        bam = OUT_DIR + "{sample}/{sample}.sorted.bam",
    output:
        tsv = OUT_DIR + "{sample}/{sample}.ref_coverage.tsv",
    log:
        OUT_DIR + "{sample}/logs/reference_coverage.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""samtools coverage "{input.bam}" > "{output.tsv}" 2> "{log}" """


rule aggregate_reference_coverage:
    """Run-level long-format table: one row per sample x reference replicon."""
    input:
        tsvs = expand(OUT_DIR + "{sample}/{sample}.ref_coverage.tsv", sample=SAMPLE_NAMES),
    output:
        summary = OUT_DIR + "reference_coverage_summary.tsv",
    log:
        OUT_DIR + "logs/aggregate_reference_coverage.log",
    params:
        exe     = script("coverage_summary", "aggregate_coverage.1.py"),
        samples = SAMPLE_NAMES,
    shell:
        "python {params.exe} -in {input.tsvs} -samples {params.samples} "
        "-name_col replicon -source reference -out {output.summary} > {log} 2>&1"
