# ======================================================================
# rules/multimer.smk :: head-to-tail multimer (concatemer) detection
# ONT samples only: it needs reads longer than the monomer.
#
# A circular monomer of length L cannot produce a read longer than L, so
# reads clearly longer than L are either true multimers or library
# chimeras.  We re-map the filtered reads to the monomeric reference with
# minimap2 (PAF, split alignments kept) and test, read by read, whether
# the alignment blocks wrap the circle seamlessly and in one orientation.
# ======================================================================

rule multimer_align:
    input:
        reads     = OUT_DIR + "{sample}/{sample}.filt.fastq.gz",
        ref_fasta = OUT_DIR + "{sample}/{sample}.ref.fasta",
    output:
        paf = OUT_DIR + "{sample}/{sample}.multimer.paf",
    wildcard_constraints:
        sample = ONT_RE,
    log:
        OUT_DIR + "{sample}/logs/multimer_align.log",
    threads: config.get("multimer_threads", 4)
    resources:
        mem_mb = config.get("multimer_mem_mb", 4000),
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        set -euo pipefail
        minimap2 -t {threads} -x map-ont -c --secondary=no \
                 "{input.ref_fasta}" "{input.reads}" \
                 > "{output.paf}" 2> "{log}"
        """


rule detect_multimers:
    input:
        paf       = OUT_DIR + "{sample}/{sample}.multimer.paf",
        ref_fasta = OUT_DIR + "{sample}/{sample}.ref.fasta",
    output:
        tsv      = OUT_DIR + "{sample}/{sample}.multimer.tsv",
        summary  = OUT_DIR + "{sample}/{sample}.multimer.summary.tsv",
        log_file = OUT_DIR + "{sample}/{sample}.multimer.log",
    wildcard_constraints:
        sample = ONT_RE,
    params:
        exe       = script_path("multimer", "multimer.2.py"),
        min_frac  = config.get("multimer_min_fraction", 0.02),
        min_reads = config.get("multimer_min_reads", 5),
        tol       = config.get("multimer_junction_tol", 0.02),
        min_qcov  = config.get("multimer_min_query_cov", 0.8),
        cut       = config.get("multimer_length_cut", 1.5),
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.exe} -paf {input.paf} -ref {input.ref_fasta} "
        "-out {output.tsv} -summary {output.summary} -min_frac {params.min_frac} "
        "-min_reads {params.min_reads} -tol {params.tol} "
        "-min_qcov {params.min_qcov} -cut {params.cut} > {output.log_file}"
