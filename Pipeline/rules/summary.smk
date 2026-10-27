# ======================================================================
# rules/summary.smk :: run-level summary table (one row per sample)
#
# Aggregation rule. Inputs are only the files that are enabled by the RUN_*
# toggles; scripts/summary.py (Snakemake `script:`) reads the small result
# files, never the FASTQs (read_stats holds the seqkit numbers).
# ======================================================================

def _run_summary_inputs():
    d = {
        "read_stats": expand(OUT_DIR + "{sample}/qc/{sample}.read_stats.tsv", sample=SAMPLE_NAMES),
        "flagstats":  expand(OUT_DIR + "{sample}/{sample}.flagstat.txt", sample=SAMPLE_NAMES),
        "ref_cov":    expand(OUT_DIR + "{sample}/{sample}.ref_coverage.tsv", sample=SAMPLE_NAMES),
        "ref_fastas": expand(OUT_DIR + "{sample}/{sample}.ref.fasta", sample=SAMPLE_NAMES),
    }
    if RUN_VARIANTS:
        d["align_snps"] = expand(OUT_DIR + "{sample}/{sample}.simple.tsv", sample=SAMPLE_NAMES)
    if RUN_GATC:
        d["annotated_gatc_align_snps"] = expand(OUT_DIR + "{sample}/{sample}.annotated.gatc.simple.tsv", sample=SAMPLE_NAMES)
    if RUN_COVERAGE:
        d["low_coverage"] = expand(OUT_DIR + "{sample}/annotated/{sample}.low_coverage.annotated.tsv", sample=SAMPLE_NAMES)
    if RUN_SV:
        d["sv"] = expand(OUT_DIR + "{sample}/annotated/{sample}.{table}.annotated.tsv",
                         sample=SAMPLE_NAMES, table=["indels", "breakpoints"])
    if RUN_MULTIMER and ONT_SAMPLES:
        d["multimer"] = expand(OUT_DIR + "{sample}/{sample}.multimer.summary.tsv", sample=ONT_SAMPLES)
    if RUN_DENOVO:
        d["assembly_infos"] = expand(OUT_DIR + "{sample}/{sample}.assembly_info.txt", sample=SAMPLE_NAMES)
    if RUN_HOTSPOTS:
        d["cassette_hotspots"] = expand(
            OUT_DIR + "{sample}/plasmidic_reads/{sample}.cassette_hotspots.tsv", sample=SAMPLE_NAMES
        )
    return d


rule run_summary:
    input:
        **_run_summary_inputs()
    output:
        summary = OUT_DIR + "summary.tsv",
    log:
        OUT_DIR + "logs/run_summary.log",
    params:
        samples = samples[["sample_name", "platform", "reads1", "reads2", "reference"]].to_dict("records"),
        out_dir = OUT_DIR,
        toggles = {
            "variants":  bool(RUN_VARIANTS),
            "gatc":      bool(RUN_GATC),
            "coverage":  bool(RUN_COVERAGE),
            "sv":        bool(RUN_SV),
            "multimer":  bool(RUN_MULTIMER and ONT_SAMPLES),
            "denovo":    bool(RUN_DENOVO),
            "hotspots":  bool(RUN_HOTSPOTS),
        },
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summary.py"
