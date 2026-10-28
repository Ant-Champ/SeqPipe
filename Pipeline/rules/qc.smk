# ======================================================================
# rules/qc.smk :: read filtering / trimming + read statistics
#   ONT          : filtlong (length / quality), as in the ONT pipeline
#   short reads  : fastp (adapter + quality trimming)
# run_read_filtering: false passes the reads through unchanged.
# depth_check: samples with too few bases for min_theoretical_depth are
# 'low_diversity' and skipped; depth_checks (checkpoint, one per run) gathers
# them so ok_samples() can decide what to run next.
# ======================================================================

FILTER_MODE = config.get("filter_mode", "best95per") if RUN_FILTERING else "none"
MIN_THEORETICAL_DEPTH = config.get("min_theoretical_depth", 5) or 0


def target_bases(wildcards):
    """100x of the reference length, for depth-capped filtering."""
    return ref_length(wildcards.sample) * 100


rule filter_reads_ont:
    input:
        fastq     = OUT_DIR + "{sample}/{sample}.fastq.gz",
        ref_fasta = OUT_DIR + "{sample}/{sample}.ref.fasta",
    output:
        OUT_DIR + "{sample}/{sample}.filt.fastq.gz",
    wildcard_constraints:
        sample = ONT_RE,
    log:
        OUT_DIR + "{sample}/logs/filter_reads.log",
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
            cp -L -- "{input.fastq}" "{output}"
            ;;
          *)
            echo "ERROR: unknown filter_mode '{params.mode}' \
(expected: best95per | best95per_max100x | none)" >&2
            exit 1
            ;;
        esac
        """


rule filter_reads_sr:
    input:
        r1 = OUT_DIR + "{sample}/{sample}_R1.fastq.gz",
        r2 = OUT_DIR + "{sample}/{sample}_R2.fastq.gz",
    output:
        r1   = OUT_DIR + "{sample}/{sample}.filt_R1.fastq.gz",
        r2   = OUT_DIR + "{sample}/{sample}.filt_R2.fastq.gz",
        json = OUT_DIR + "{sample}/qc/{sample}.fastp.json",
        html = OUT_DIR + "{sample}/qc/{sample}.fastp.html",
    wildcard_constraints:
        sample = SR_RE,
    log:
        OUT_DIR + "{sample}/logs/filter_reads.log",
    params:
        run        = "yes" if RUN_FILTERING else "no",
        min_length = config.get("fastp_min_length", 30),
        min_qual   = config.get("fastp_min_qual", 20),
    threads: config.get("fastp_threads", 4)
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        set -euo pipefail
        if [ "{params.run}" = "yes" ]; then
            fastp -i "{input.r1}" -I "{input.r2}" -o "{output.r1}" -O "{output.r2}" \
                  --detect_adapter_for_pe --cut_right --cut_right_mean_quality {params.min_qual} \
                  --length_required {params.min_length} --thread {threads} \
                  --json "{output.json}" --html "{output.html}" > "{log}" 2>&1
        else
            cp -L -- "{input.r1}" "{output.r1}"
            cp -L -- "{input.r2}" "{output.r2}"
            echo '{{"skipped": "run_read_filtering is false"}}' > "{output.json}"
            echo 'read filtering skipped' > "{output.html}"
        fi
        """


rule read_stats:
    """Read count / bases / N50 of raw and filtered reads (seqkit stats) -> summary."""
    input:
        raw  = lambda wc: staged_reads(wc.sample),
        filt = filt_reads,
    output:
        OUT_DIR + "{sample}/qc/{sample}.read_stats.tsv",
    log:
        OUT_DIR + "{sample}/logs/read_stats.log",
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        r"""seqkit stats -a -T -j {threads} {input.raw} {input.filt} > "{output}" 2> "{log}" """


rule depth_check:
    """Theoretical depth = filtered bases / reference length. Below
    min_theoretical_depth the sample is 'low_diversity': it is not aligned or
    analysed further (red banner in the report). Off when min_theoretical_depth is 0."""
    input:
        stats = OUT_DIR + "{sample}/qc/{sample}.read_stats.tsv",
        fai   = OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
    output:
        OUT_DIR + "{sample}/qc/{sample}.depth_check.tsv",
    params:
        filt      = lambda wc: [os.path.basename(p) for p in filt_reads(wc)],
        min_depth = MIN_THEORETICAL_DEPTH,
    run:
        import csv
        with open(input.stats, newline="") as fh:
            rows = {os.path.basename(r["file"]): r for r in csv.DictReader(fh, delimiter="\t")}
        bases = sum(int(rows[f]["sum_len"].replace(",", "")) for f in params.filt)
        with open(input.fai) as fh:
            ref_len = sum(int(l.split("\t")[1]) for l in fh if l.strip())
        depth = bases / ref_len if ref_len else 0.0
        status = "ok" if depth >= params.min_depth else "low_diversity"
        with open(output[0], "w") as out:
            out.write("sample\tfilt_bases\tref_length\ttheoretical_depth\tmin_theoretical_depth\tstatus\n")
            out.write(f"{wildcards.sample}\t{bases}\t{ref_len}\t{depth:.2f}\t{params.min_depth}\t{status}\n")


checkpoint depth_checks:
    """All depth_check tables in one: the DAG is rebuilt from it (ok_samples())."""
    input:
        expand(OUT_DIR + "{sample}/qc/{sample}.depth_check.tsv", sample=SAMPLE_NAMES),
    output:
        OUT_DIR + "depth_check.tsv",
    run:
        with open(output[0], "w") as out:
            for i, path in enumerate(input):
                with open(path) as fh:
                    lines = fh.read().splitlines()
                out.write("\n".join(lines if i == 0 else lines[1:]) + "\n")
