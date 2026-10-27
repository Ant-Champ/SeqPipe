# ======================================================================
# rules/qc.smk :: read length/quality filtering
# ======================================================================

def target_bases(wildcards):
    """100x of the reference length, for depth-capped filtering."""
    return ref_length(wildcards.sample) * 100


rule initial_trimming:
    input:
        fastq     = OUT_DIR + "{sample}/{sample}.fastq.gz",
        ref_fasta = OUT_DIR + "{sample}/{sample}.ref.fasta",
    output:
        OUT_DIR + "{sample}/{sample}.filt.fastq.gz",
    log:
        OUT_DIR + "{sample}/logs/initial_trimming.log",
    params:
        mode         = FILTER_MODE,
        min_length   = config.get("filtlong_min_length", 500),
        keep_percent = config.get("filtlong_keep_percent", 95),
        target_bases = target_bases,
    threads: 1
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        set -euo pipefail
        case "{params.mode}" in
          best95per)
            filtlong --min_length {params.min_length} \
                     --keep_percent {params.keep_percent} \
                     "{input.fastq}" 2> "{log}" | gzip > "{output}"
            ;;
          best95per_max100x)
            filtlong --min_length {params.min_length} \
                     --target_bases {params.target_bases} \
                     --keep_percent {params.keep_percent} \
                     "{input.fastq}" 2> "{log}" | gzip > "{output}"
            ;;
          none)
            cp -- "{input.fastq}" "{output}"
            ;;
          *)
            echo "ERROR: unknown filter_mode '{params.mode}' \
(expected: best95per | best95per_max100x | none)" >&2
            exit 1
            ;;
        esac
        """
