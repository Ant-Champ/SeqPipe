# ======================================================================
# rules/sv.smk :: large indels + breakpoints (from the Element alignqc
# pipeline), for deleted / inserted / rearranged regions.
#   indels      : CIGAR insertions/deletions >= min_indel bp
#   breakpoints : soft-clip pile-ups >= min_clip bp; n_split counts reads
#                 that also have a supplementary alignment (SA tag)
# Thresholds are per read type (config sv_scan.ont / sv_scan.short):
# ONT reads have long noisy ends, so they need larger clips.
# ======================================================================

SV_DEFAULTS = {
    "ont":   {"min_indel": 30, "min_clip": 200, "min_support": 3, "cluster_tol": 20},
    "short": {"min_indel": 10, "min_clip": 20,  "min_support": 3, "cluster_tol": 5},
}


def sv_param(key):
    def _get(wildcards):
        kind = "ont" if is_ont(wildcards.sample) else "short"
        return (config.get("sv_scan", {}) or {}).get(kind, {}).get(key, SV_DEFAULTS[kind][key])
    return _get


rule sv_scan:
    input:
        bam = OUT_DIR + "{sample}/{sample}.sorted.bam",
        bai = OUT_DIR + "{sample}/{sample}.sorted.bam.bai",
    output:
        indels      = OUT_DIR + "{sample}/variants/{sample}.indels.tsv",
        breakpoints = OUT_DIR + "{sample}/variants/{sample}.breakpoints.tsv",
    log:
        OUT_DIR + "{sample}/logs/sv_scan.log",
    params:
        exe         = script_path("sv_scan", "sv_scan.py"),
        min_mapq    = config.get("min_mapq", 20),
        min_indel   = sv_param("min_indel"),
        min_clip    = sv_param("min_clip"),
        min_support = sv_param("min_support"),
        cluster_tol = sv_param("cluster_tol"),
        edge        = config.get("sv_ignore_contig_ends", 10),
    conda:
        "../envs/regions.yaml"
    shell:
        "python {params.exe} --bam {input.bam} "
        "--out-indels {output.indels} --out-breakpoints {output.breakpoints} "
        "--min-mapq {params.min_mapq} --min-indel {params.min_indel} "
        "--min-clip {params.min_clip} --min-support {params.min_support} "
        "--cluster-tol {params.cluster_tol} --edge {params.edge} > {log} 2>&1"


rule annotate_sv:
    input:
        tsv = OUT_DIR + "{sample}/variants/{sample}.{table}.tsv",
        gb  = OUT_DIR + "{sample}/{sample}.ref.gb",
    output:
        OUT_DIR + "{sample}/annotated/{sample}.{table}.annotated.tsv",
    wildcard_constraints:
        table = "indels|breakpoints",
    log:
        OUT_DIR + "{sample}/logs/annotate_{table}.log",
    params:
        exe        = script_path("annotate_regions", "annotate_regions.py"),
        skip_types = config.get("annotate_skip_types", "source,gene"),
    conda:
        "../envs/regions.yaml"
    shell:
        "python {params.exe} --genbank {input.gb} --in {input.tsv} --out {output} "
        "--skip-types {params.skip_types} > {log} 2>&1"
