# ======================================================================
# rules/variants.smk :: SNP / small-indel calling -> {sample}.simple.tsv
#   ONT          : medaka_variant + simplify.4.py
#   short reads  : bcftools mpileup | call (haploid) + simplify_bcftools.py
# Both tables have the same columns (CHROM POS REF ALT %Reads MutationType),
# so the GATC / gene annotation chain below is shared.
# ======================================================================

rule variant_calling_ont:
    input:
        reads     = OUT_DIR + "{sample}/{sample}.filt.fastq.gz",
        ref_fasta = OUT_DIR + "{sample}/{sample}.ref.fasta",
    output:
        medaka_dir = directory(OUT_DIR + "{sample}/medaka"),
        vcf        = OUT_DIR + "{sample}/{sample}.vcf",
    wildcard_constraints:
        sample = ONT_RE,
    log:
        OUT_DIR + "{sample}/logs/variant_calling.log",
    threads: config.get("medaka_threads", 4)
    resources:
        mem_mb = config.get("medaka_mem_mb", 8000),
    conda:
        "../envs/medaka.yaml"
    shell:
        r"""
        set -euo pipefail
        medaka_variant -i "{input.reads}" -r "{input.ref_fasta}" \
                       -o "{output.medaka_dir}" -t {threads} -s > "{log}" 2>&1
        cp "{output.medaka_dir}/medaka.annotated.vcf" "{output.vcf}"
        """


rule filter_variants_ont:
    input:
        vcf = OUT_DIR + "{sample}/{sample}.vcf",
    output:
        tsv = OUT_DIR + "{sample}/{sample}.simple.tsv",
    wildcard_constraints:
        sample = ONT_RE,
    log:
        OUT_DIR + "{sample}/logs/filter_variants.log",
    params:
        exe      = script_path("simplify", "simplify.4.py"),
        min_qual = config.get("simplify_min_qual", 50),
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.exe} -in {input.vcf} -out {output.tsv} "
        "-min {params.min_qual} > {log} 2>&1"


rule variant_calling_sr:
    input:
        bam       = OUT_DIR + "{sample}/{sample}.sorted.bam",
        ref_fasta = OUT_DIR + "{sample}/{sample}.ref.fasta",
        fai       = OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
    output:
        vcf = OUT_DIR + "{sample}/{sample}.vcf",
    wildcard_constraints:
        sample = SR_RE,
    log:
        OUT_DIR + "{sample}/logs/variant_calling.log",
    params:
        min_mapq  = config.get("min_mapq", 20),
        min_baseq = config.get("bcftools_min_baseq", 20),
        max_depth = config.get("bcftools_max_depth", 10000),
    threads: 2
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        set -euo pipefail
        bcftools mpileup -f "{input.ref_fasta}" -q {params.min_mapq} -Q {params.min_baseq} \
                 -d {params.max_depth} -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP \
                 -Ou "{input.bam}" 2> "{log}" \
          | bcftools call -m -v --ploidy 1 -Ov -o "{output.vcf}" 2>> "{log}"
        """


rule filter_variants_sr:
    input:
        vcf = OUT_DIR + "{sample}/{sample}.vcf",
    output:
        tsv = OUT_DIR + "{sample}/{sample}.simple.tsv",
    wildcard_constraints:
        sample = SR_RE,
    log:
        OUT_DIR + "{sample}/logs/filter_variants.log",
    params:
        exe      = script_path("simplify_bcftools", "simplify_bcftools.py"),
        min_qual = config.get("bcftools_min_qual", 30),
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.exe} -in {input.vcf} -out {output.tsv} "
        "-min {params.min_qual} > {log} 2>&1"


# ----------------------------------------------------------------------
# GATC-site + gene context, then the multimer verdict (ONT only)
# ----------------------------------------------------------------------
if RUN_GATC:

    rule id_gatc_variants:
        input:
            ref_fasta  = OUT_DIR + "{sample}/{sample}.ref.fasta",
            align_snps = OUT_DIR + "{sample}/{sample}.simple.tsv",
        output:
            gatc_align_snps = OUT_DIR + "{sample}/{sample}.gatc.simple.tsv",
        log:
            OUT_DIR + "{sample}/logs/id_gatc_variants.log",
        params:
            exe = script_path("id_gatc", "id_GATC.2.py"),
        conda:
            "../envs/python.yaml"
        shell:
            "python {params.exe} "
            "-ref_fasta {input.ref_fasta} "
            "-variants {input.align_snps} "
            "-out {output.gatc_align_snps} > {log} 2>&1"


    # Intermediate table; the multimer column is added by add_multimer_flag_*,
    # which produces the final filename.
    rule annotate_gatc_variants:
        input:
            ref             = OUT_DIR + "{sample}/{sample}.ref.gb",
            gatc_align_snps = OUT_DIR + "{sample}/{sample}.gatc.simple.tsv",
        output:
            annotated = temp(OUT_DIR + "{sample}/{sample}.annotated.gatc.nomulti.tsv"),
        log:
            OUT_DIR + "{sample}/logs/annotate_gatc_variants.log",
        params:
            exe = script_path("annotate", "annotate.2.py"),
        conda:
            "../envs/python.yaml"
        shell:
            "python {params.exe} "
            "-t {input.gatc_align_snps} "
            "-a {input.ref} "
            "-o {output.annotated} > {log} 2>&1"


    # Samples that get the multimer column: ONT samples when the multimer
    # analysis is on. Everyone else gets the annotated table unchanged.
    MULTIMER_SAMPLES = ONT_SAMPLES if RUN_MULTIMER else []

    rule add_multimer_flag:
        """Append the sample-level multimer verdict (yes/no) to the final table."""
        input:
            annotated = OUT_DIR + "{sample}/{sample}.annotated.gatc.nomulti.tsv",
            summary   = OUT_DIR + "{sample}/{sample}.multimer.summary.tsv",
        output:
            annotated_gatc_align_snps = OUT_DIR + "{sample}/{sample}.annotated.gatc.simple.tsv",
        wildcard_constraints:
            sample = _alternation(MULTIMER_SAMPLES),
        log:
            OUT_DIR + "{sample}/logs/add_multimer_flag.log",
        params:
            exe   = script_path("multimer_flag", "add_multimer_column.1.py"),
            col   = config.get("multimer_column_name", "multimer"),
            extra = "-extra" if config.get("multimer_extra_columns", False) else "",
        conda:
            "../envs/python.yaml"
        shell:
            "python {params.exe} "
            "-in {input.annotated} "
            "-summary {input.summary} "
            "-col {params.col} {params.extra} "
            "-out {output.annotated_gatc_align_snps} > {log} 2>&1"


    rule skip_multimer_flag:
        """No multimer analysis for this sample: pass the annotated table through unchanged."""
        input:
            annotated = OUT_DIR + "{sample}/{sample}.annotated.gatc.nomulti.tsv",
        output:
            annotated_gatc_align_snps = OUT_DIR + "{sample}/{sample}.annotated.gatc.simple.tsv",
        wildcard_constraints:
            sample = _alternation([s for s in SAMPLE_NAMES if s not in MULTIMER_SAMPLES]),
        shell:
            r"""cp -- "{input.annotated}" "{output.annotated_gatc_align_snps}" """
