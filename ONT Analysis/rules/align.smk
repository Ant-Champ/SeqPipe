# ======================================================================
# rules/align.smk :: reference alignment + variant calling
# ======================================================================

rule variant_calling:
    input:
        reads     = OUT_DIR + "{sample}/{sample}.filt.fastq.gz",
        ref_fasta = OUT_DIR + "{sample}/{sample}.ref.fasta",
    output:
        medaka_dir = directory(OUT_DIR + "{sample}/medaka"),
        sorted_bam = OUT_DIR + "{sample}/{sample}.sorted.bam",
        bam_index  = OUT_DIR + "{sample}/{sample}.sorted.bam.bai",
        vcf        = OUT_DIR + "{sample}/{sample}.vcf",
    log:
        OUT_DIR + "{sample}/logs/variant_calling.log",
    threads: config.get("medaka_threads", 4)
    resources:
        mem_mb  = config.get("medaka_mem_mb", 8000),
        runtime = config.get("medaka_runtime", 120),   # minutes (SLURM)
    conda:
        "../envs/medaka.yaml"
    shell:
        r"""
        set -euo pipefail
        medaka_variant -i "{input.reads}" -r "{input.ref_fasta}" \
                       -o "{output.medaka_dir}" -t {threads} -s > "{log}" 2>&1
        mv "{output.medaka_dir}/calls_to_ref.bam"       "{output.sorted_bam}"
        mv "{output.medaka_dir}/calls_to_ref.bam.bai"   "{output.bam_index}"
        mv "{output.medaka_dir}/medaka.annotated.vcf"   "{output.vcf}"
        """


rule filter_variants:
    input:
        vcf = OUT_DIR + "{sample}/{sample}.vcf",
    output:
        tsv = OUT_DIR + "{sample}/{sample}.simple.tsv",
    log:
        OUT_DIR + "{sample}/logs/filter_variants.log",
    params:
        exe     = script("simplify", "simplify.4.py"),
        min_qual= config.get("simplify_min_qual", 50),
    shell:
        "python {params.exe} -in {input.vcf} -out {output.tsv} "
        "-min {params.min_qual} > {log} 2>&1"
