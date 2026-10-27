# ======================================================================
# rules/align.smk :: alignment to the sample's own reference
#   ONT          : minimap2 -x map-ont
#   short reads  : bwa-mem2 mem (paired)
# Both keep supplementary alignments (SA tags) for the breakpoint scan,
# and write {sample}.sorted.bam(.bai) + {sample}.flagstat.txt.
# ======================================================================

rule align_ont:
    input:
        reads     = OUT_DIR + "{sample}/{sample}.filt.fastq.gz",
        ref_fasta = OUT_DIR + "{sample}/{sample}.ref.fasta",
    output:
        bam = OUT_DIR + "{sample}/{sample}.sorted.bam",
        bai = OUT_DIR + "{sample}/{sample}.sorted.bam.bai",
    wildcard_constraints:
        sample = ONT_RE,
    log:
        OUT_DIR + "{sample}/logs/align.log",
    params:
        rg = lambda wc: f"@RG\\tID:{wc.sample}\\tSM:{wc.sample}\\tPL:ONT",
    threads: config.get("align_threads", 8)
    resources:
        mem_mb = config.get("align_mem_mb", 16000),
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        set -euo pipefail
        tmp="${{TMPDIR:-/tmp}}"
        minimap2 -ax map-ont -t {threads} -R '{params.rg}' "{input.ref_fasta}" "{input.reads}" 2> "{log}" \
          | samtools sort -@ {threads} -T "$tmp/{wildcards.sample}.align" -o "{output.bam}" - 2>> "{log}"
        samtools index "{output.bam}"
        """


rule bwa_index:
    input:
        OUT_DIR + "{sample}/{sample}.ref.fasta",
    output:
        multiext(OUT_DIR + "{sample}/{sample}.ref.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    wildcard_constraints:
        sample = SR_RE,
    log:
        OUT_DIR + "{sample}/logs/bwa_index.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""bwa-mem2 index "{input}" > "{log}" 2>&1"""


rule align_sr:
    input:
        r1        = OUT_DIR + "{sample}/{sample}.filt_R1.fastq.gz",
        r2        = OUT_DIR + "{sample}/{sample}.filt_R2.fastq.gz",
        ref_fasta = OUT_DIR + "{sample}/{sample}.ref.fasta",
        index     = multiext(OUT_DIR + "{sample}/{sample}.ref.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        bam = OUT_DIR + "{sample}/{sample}.sorted.bam",
        bai = OUT_DIR + "{sample}/{sample}.sorted.bam.bai",
    wildcard_constraints:
        sample = SR_RE,
    log:
        OUT_DIR + "{sample}/logs/align.log",
    params:
        # read-group platform tag: ELEMENT or ILLUMINA
        rg = lambda wc: f"@RG\\tID:{wc.sample}\\tSM:{wc.sample}\\tPL:{samples.at[wc.sample, 'platform'].upper()}",
    threads: config.get("align_threads", 8)
    resources:
        mem_mb = config.get("align_mem_mb", 16000),
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        set -euo pipefail
        tmp="${{TMPDIR:-/tmp}}"
        bwa-mem2 mem -t {threads} -R '{params.rg}' "{input.ref_fasta}" "{input.r1}" "{input.r2}" 2> "{log}" \
          | samtools sort -@ {threads} -T "$tmp/{wildcards.sample}.align" -o "{output.bam}" - 2>> "{log}"
        samtools index "{output.bam}"
        """


rule flagstat:
    input:
        OUT_DIR + "{sample}/{sample}.sorted.bam",
    output:
        OUT_DIR + "{sample}/{sample}.flagstat.txt",
    conda:
        "../envs/align.yaml"
    shell:
        r"""samtools flagstat "{input}" > "{output}" """


rule reference_coverage:
    """samtools coverage: breadth / mean depth per reference replicon."""
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
