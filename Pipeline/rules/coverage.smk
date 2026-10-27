# ======================================================================
# rules/coverage.smk :: coverage per replicon + low-coverage regions
#   reference_coverage : samtools coverage (breadth/depth per replicon),
#                        aggregated into reference_coverage_summary.tsv
#   lowcov             : low-coverage regions (from the Element alignqc
#                        pipeline): a base is "low" if its depth is below
#                        max(lowcov_min_depth, lowcov_frac x median depth)
# Same for ONT and short reads (all from {sample}.sorted.bam).
# ======================================================================

# (rule reference_coverage lives in rules/align.smk: the summary always needs it)

rule aggregate_reference_coverage:
    """Run-level long-format table: one row per sample x reference replicon."""
    input:
        tsvs = expand(OUT_DIR + "{sample}/{sample}.ref_coverage.tsv", sample=SAMPLE_NAMES),
    output:
        summary = OUT_DIR + "reference_coverage_summary.tsv",
    log:
        OUT_DIR + "logs/aggregate_reference_coverage.log",
    params:
        exe     = script_path("coverage_summary", "aggregate_coverage.1.py"),
        samples = SAMPLE_NAMES,
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.exe} -in {input.tsvs} -samples {params.samples} "
        "-name_col replicon -source reference -out {output.summary} > {log} 2>&1"


rule lowcov:
    input:
        bam = OUT_DIR + "{sample}/{sample}.sorted.bam",
        fai = OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
    output:
        summary = OUT_DIR + "{sample}/coverage/{sample}.coverage_summary.tsv",
        low     = OUT_DIR + "{sample}/coverage/{sample}.low_coverage.tsv",
    log:
        OUT_DIR + "{sample}/logs/lowcov.log",
    params:
        exe       = script_path("lowcov", "lowcov.py"),
        min_mapq  = config.get("min_mapq", 20),
        low_depth = config.get("lowcov_min_depth", 10),
        low_frac  = config.get("lowcov_frac", 0.25),
        min_len   = config.get("lowcov_min_len", 20),
        merge_gap = config.get("lowcov_merge_gap", 10),
    conda:
        "../envs/regions.yaml"
    shell:
        "python {params.exe} --bam {input.bam} --fai {input.fai} "
        "--out-summary {output.summary} --out-low {output.low} "
        "--min-mapq {params.min_mapq} --low-depth {params.low_depth} "
        "--low-frac {params.low_frac} --low-min-len {params.min_len} "
        "--low-merge-gap {params.merge_gap} > {log} 2>&1"


rule annotate_low_coverage:
    input:
        tsv = OUT_DIR + "{sample}/coverage/{sample}.low_coverage.tsv",
        gb  = OUT_DIR + "{sample}/{sample}.ref.gb",
    output:
        OUT_DIR + "{sample}/annotated/{sample}.low_coverage.annotated.tsv",
    log:
        OUT_DIR + "{sample}/logs/annotate_low_coverage.log",
    params:
        exe        = script_path("annotate_regions", "annotate_regions.py"),
        skip_types = config.get("annotate_skip_types", "source,gene"),
    conda:
        "../envs/regions.yaml"
    shell:
        "python {params.exe} --genbank {input.gb} --in {input.tsv} --out {output} "
        "--skip-types {params.skip_types} > {log} 2>&1"
