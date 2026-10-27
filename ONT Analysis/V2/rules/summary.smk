# ======================================================================
# rules/summary.smk :: run-level summary table (one row per sample)
#
# This is an aggregation rule. Only included in rule all's targets when
# RUN_ALIGNMENT is on (see Snakefile) -- everything it reports comes from
# that branch, so its inputs below only ask for what's actually enabled.
# ======================================================================

def _run_summary_inputs():
    d = {
        "sheet":       config["sample_sheet"],
        "ref_fastas":  expand(OUT_DIR + "{sample}/{sample}.ref.fasta",  sample=SAMPLE_NAMES),
        "simple_tsvs": expand(OUT_DIR + "{sample}/{sample}.simple.tsv", sample=SAMPLE_NAMES),
    }
    if RUN_COVERAGE:
        d["cov_txt"]  = expand(OUT_DIR + "{sample}/{sample}.cov.txt", sample=SAMPLE_NAMES)
        d["cov_logs"] = expand(OUT_DIR + "{sample}/{sample}.cov.log", sample=SAMPLE_NAMES)
    if RUN_MULTIMER:
        d["multimer"] = expand(OUT_DIR + "{sample}/{sample}.multimer.summary.tsv", sample=SAMPLE_NAMES)
    if RUN_GATC:
        # not read by the script, but forces the summary to be the last step
        d["annotated"] = expand(OUT_DIR + "{sample}/{sample}.annotated.gatc.simple.tsv", sample=SAMPLE_NAMES)
    return d


rule run_summary:
    input:
        **_run_summary_inputs()
    output:
        summary = OUT_DIR + "summary.tsv",
    log:
        OUT_DIR + "logs/run_summary.log",
    params:
        exe        = script("summary", "summary.2.py"),
        out_dir    = OUT_DIR,
        fastq_dir  = FASTQ_DIR,
        fastq_base = FASTQ_BASE_NAME,
    conda:
        "../envs/summary.yaml"
    shell:
        r"""
        python {params.exe} \
            -out_dir "{params.out_dir}" \
            -sample_sheet "{input.sheet}" \
            -fastq_dir "{params.fastq_dir}" \
            -fastq_base_name "{params.fastq_base}" \
            -out "{output.summary}" > "{log}" 2>&1
        """
