# ======================================================================
# rules/coverage.smk :: coverage per replicon + low- and high-coverage regions
#   reference_coverage : samtools coverage (breadth/depth per replicon),
#                        aggregated into reference_coverage_summary.tsv
#   lowcov             : low-coverage regions (from the Element alignqc
#                        pipeline): a base is "low" if its depth is below
#                        max(lowcov_min_depth, lowcov_frac x median depth);
#                        and high-coverage regions (duplications): windows
#                        above highcov_frac x the median of their replicon
#   mappability        : (short reads, run_mappability_filter) regions of
#                        each reference where PE reads of mappability_read_len
#                        cannot map uniquely; ignored by the low-coverage scan
# Same for ONT and short reads (all from {sample}.sorted.bam).
# ======================================================================

RUN_MAPPABILITY = config.get("run_mappability_filter", True)
MAPPABILITY_DIR = os.path.join(config.get("mappability_dir") or OUT_DIR + "mappability", "")
MAP_READ_LEN    = config.get("mappability_read_len", 75)
MAP_INSERT      = config.get("mappability_insert", 300)
MAP_MIN_MAPQ    = config.get("min_mapq", 20)

# one mask per reference GenBank file, named after the file (+ the read settings,
# so changing them makes a new mask; the folder can be shared between runs)
REF_KEYS = {s: re.sub(r"[^\w.-]", "_", os.path.splitext(os.path.basename(ref_gb_path(s)))[0])
            for s in SR_SAMPLES}
REF_KEY_PATH = {}
for _s, _k in REF_KEYS.items():
    if REF_KEY_PATH.setdefault(_k, ref_gb_path(_s)) != ref_gb_path(_s):
        raise WorkflowError(f"two different reference files are both named {_k}: rename one")


def mappability_bed(sample):
    """Unmappable-regions BED of a sample's reference, or None (ONT / filter off)."""
    if not RUN_MAPPABILITY or is_ont(sample):
        return None
    return MAPPABILITY_DIR + f"{REF_KEYS[sample]}.pe{MAP_READ_LEN}_ins{MAP_INSERT}_q{MAP_MIN_MAPQ}.unmappable.bed"


if RUN_MAPPABILITY and SR_SAMPLES:

    rule mappability:
        """Simulated perfect PE reads from the reference, mapped back with bwa-mem2:
        bases they cannot cover uniquely (MAPQ < min_mapq) = unmappable."""
        input:
            gb = lambda wc: REF_KEY_PATH[wc.ref],
        output:
            bed = MAPPABILITY_DIR + "{ref}.pe{read_len}_ins{ins}_q{mapq}.unmappable.bed",
        wildcard_constraints:
            ref = _alternation(sorted(REF_KEY_PATH)),
            read_len = r"\d+", ins = r"\d+", mapq = r"\d+",
        log:
            MAPPABILITY_DIR + "logs/{ref}.pe{read_len}_ins{ins}_q{mapq}.log",
        params:
            exe      = script_path("mappability", "mappability.py"),
            step     = config.get("mappability_step", 2),
            min_frac = config.get("mappability_min_frac", 0.5),
        threads: config.get("align_threads", 8)
        resources:
            mem_mb = config.get("align_mem_mb", 16000),
        conda:
            "../envs/mappability.yaml"
        shell:
            "python {params.exe} --genbank {input.gb} --out-bed {output.bed} "
            "--read-len {wildcards.read_len} --insert {wildcards.ins} --min-mapq {wildcards.mapq} "
            "--step {params.step} --min-frac {params.min_frac} --threads {threads} > {log} 2>&1"

# (rule reference_coverage lives in rules/align.smk: the summary always needs it)

rule aggregate_reference_coverage:
    """Run-level long-format table: one row per sample x reference replicon."""
    input:
        tsvs = lambda wc: expand(OUT_DIR + "{sample}/{sample}.ref_coverage.tsv", sample=ok_samples()),
    output:
        summary = OUT_DIR + "reference_coverage_summary.tsv",
    log:
        OUT_DIR + "logs/aggregate_reference_coverage.log",
    params:
        exe     = script_path("coverage_summary", "aggregate_coverage.1.py"),
        samples = lambda wc: ok_samples(),
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.exe} -in {input.tsvs} -samples {params.samples} "
        "-name_col replicon -source reference -out {output.summary} > {log} 2>&1"


rule lowcov:
    input:
        bam = OUT_DIR + "{sample}/{sample}.sorted.bam",
        fai = OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
        mask = lambda wc: mappability_bed(wc.sample) or [],
    output:
        summary = OUT_DIR + "{sample}/coverage/{sample}.coverage_summary.tsv",
        low     = OUT_DIR + "{sample}/coverage/{sample}.low_coverage.tsv",
        high    = OUT_DIR + "{sample}/coverage/{sample}.high_coverage.tsv",
    log:
        OUT_DIR + "{sample}/logs/lowcov.log",
    params:
        exe       = script_path("lowcov", "lowcov.py"),
        min_mapq  = config.get("min_mapq", 20),
        low_depth = config.get("lowcov_min_depth", 10),
        low_frac  = config.get("lowcov_frac", 0.25),
        min_len   = config.get("lowcov_min_len", 20),
        merge_gap = config.get("lowcov_merge_gap", 10),
        smooth    = config.get("lowcov_smooth", 50),
        max_thr   = config.get("lowcov_max_thr_frac", 0.5),
        region_max = config.get("lowcov_region_max_frac", 0.25),
        mask      = lambda wc: f"--mask {mappability_bed(wc.sample)}" if mappability_bed(wc.sample) else "",
        high_frac      = config.get("highcov_frac", 1.7),
        high_window    = config.get("highcov_window", 100),
        high_min_len   = config.get("highcov_min_len", 1000),
        high_merge_gap = config.get("highcov_merge_gap", 300),
    conda:
        "../envs/regions.yaml"
    shell:
        "python {params.exe} --bam {input.bam} --fai {input.fai} "
        "--out-summary {output.summary} --out-low {output.low} --out-high {output.high} "
        "--min-mapq {params.min_mapq} --low-depth {params.low_depth} "
        "--low-frac {params.low_frac} --low-min-len {params.min_len} "
        "--low-merge-gap {params.merge_gap} --low-smooth {params.smooth} "
        "--low-max-thr-frac {params.max_thr} --low-region-max-frac {params.region_max} {params.mask} "
        "--high-frac {params.high_frac} "
        "--high-window {params.high_window} --high-min-len {params.high_min_len} "
        "--high-merge-gap {params.high_merge_gap} > {log} 2>&1"


rule annotate_coverage_regions:
    """Genes hit by the low- / high-coverage regions."""
    input:
        tsv = OUT_DIR + "{sample}/coverage/{sample}.{level}_coverage.tsv",
        gb  = OUT_DIR + "{sample}/{sample}.ref.gb",
    output:
        OUT_DIR + "{sample}/annotated/{sample}.{level}_coverage.annotated.tsv",
    wildcard_constraints:
        level = "low|high",
    log:
        OUT_DIR + "{sample}/logs/annotate_{level}_coverage.log",
    params:
        exe        = script_path("annotate_regions", "annotate_regions.py"),
        skip_types = config.get("annotate_skip_types", "source,gene"),
    conda:
        "../envs/regions.yaml"
    shell:
        "python {params.exe} --genbank {input.gb} --in {input.tsv} --out {output} "
        "--skip-types {params.skip_types} > {log} 2>&1"
