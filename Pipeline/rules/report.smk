# ======================================================================
# rules/report.smk :: one self-contained HTML report for the run
# (run_html_report; default on)
#
# depth_profile bins the per-base depth of each sample (samtools depth -aa)
# into at most report_max_bins bins per replicon, for the interactive
# coverage plots. html_report (scripts/report.py) reads the same small
# result files as the summary, plus the depth bins, the GenBank genes and
# the hotspot pictures (run_hotspot_pictures), and writes <out_dir>/report.html:
# one file, no internet needed to open it, two views (by sample / by
# analysis). Links to bigger files (fastp reports, assemblies) are relative,
# so keep report.html inside <out_dir>.
# ======================================================================

rule depth_profile:
    """Mean depth per bin along every replicon (contigs with no reads included)."""
    input:
        bam = OUT_DIR + "{sample}/{sample}.sorted.bam",
        bai = OUT_DIR + "{sample}/{sample}.sorted.bam.bai",
        fai = OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
    output:
        OUT_DIR + "{sample}/coverage/{sample}.depth_bins.tsv",
    log:
        OUT_DIR + "{sample}/logs/depth_profile.log",
    params:
        max_bins = config.get("report_max_bins", 20000),
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        set -euo pipefail
        samtools depth -aa "{input.bam}" 2> "{log}" \
        | awk -v maxb={params.max_bins} 'BEGIN {{ OFS = "\t" }}
            NR == FNR {{ order[++nc] = $1; len[$1] = $2
                         bs[$1] = int(($2 + maxb - 1) / maxb); if (bs[$1] < 1) bs[$1] = 1; next }}
            {{ k = $1 SUBSEP int(($2 - 1) / bs[$1]); s[k] += $3; n[k]++ }}
            END {{ print "#contig", "start", "end", "mean_depth"
                   for (i = 1; i <= nc; i++) {{
                       c = order[i]; nb = int((len[c] - 1) / bs[c]) + 1
                       for (b = 0; b < nb; b++) {{
                           k = c SUBSEP b; e = (b + 1) * bs[c]; if (e > len[c]) e = len[c]
                           printf "%s\t%d\t%d\t%.1f\n", c, b * bs[c] + 1, e, (n[k] ? s[k] / n[k] : 0)
                       }}
                   }} }}' "{input.fai}" - > "{output}"
        """


RUN_REPORT_READS = config.get("report_reads", True)


def _read_event_tables(wildcards):
    """Tables whose positions get reads in the report (larger references)."""
    s, ann = wildcards.sample, OUT_DIR + "{0}/annotated/{0}".format(wildcards.sample)
    t = []
    if RUN_VARIANTS:
        t.append(OUT_DIR + f"{s}/{s}.simple.tsv")
    if RUN_SV:
        t += [ann + ".indels.annotated.tsv", ann + ".breakpoints.annotated.tsv"]
    if RUN_COVERAGE:
        t += [ann + ".low_coverage.annotated.tsv", ann + ".high_coverage.annotated.tsv"]
    if RUN_CASSETTES:
        t += [ann + ".cassette_hotspots.annotated.tsv", ann + ".cassette_insertions.annotated.tsv"]
    return t


rule read_snippets:
    """A browser-style sample of the aligned reads for the report: the whole
    reference when small (plasmids...), else windows around every event."""
    input:
        bam    = OUT_DIR + "{sample}/{sample}.sorted.bam",
        bai    = OUT_DIR + "{sample}/{sample}.sorted.bam.bai",
        fasta  = OUT_DIR + "{sample}/{sample}.ref.fasta",
        fai    = OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
        events = _read_event_tables,
    output:
        OUT_DIR + "{sample}/coverage/{sample}.reads.json",
    log:
        OUT_DIR + "{sample}/logs/read_snippets.log",
    params:
        exe       = script_path("read_snippets", "read_snippets.py"),
        whole_max = config.get("report_reads_whole_ref_max", 200000),
        flank     = config.get("report_reads_flank", 500),
        max_rows  = config.get("report_reads_max_rows", 40),
        min_mapq  = config.get("min_mapq", 20),
    conda:
        "../envs/regions.yaml"
    shell:
        "python {params.exe} --bam {input.bam} --fasta {input.fasta} --events {input.events} "
        "--out {output} --whole-ref-max {params.whole_max} --flank {params.flank} "
        "--max-rows {params.max_rows} --min-mapq {params.min_mapq} > {log} 2>&1"


def _report_inputs(wildcards=None):
    d = dict(_run_summary_inputs())
    OK = ok_samples()
    d["summary"] = OUT_DIR + "summary.tsv"
    d["depth"] = expand(OUT_DIR + "{sample}/coverage/{sample}.depth_bins.tsv", sample=OK)
    if RUN_REPORT_READS:
        d["reads"] = expand(OUT_DIR + "{sample}/coverage/{sample}.reads.json", sample=OK)
    d["ref_gb"] = expand(OUT_DIR + "{sample}/{sample}.ref.gb", sample=SAMPLE_NAMES)
    if RUN_COVERAGE or RUN_PLASMID_SCREEN:
        d["coverage_summary"] = OUT_DIR + "coverage_summary.tsv"
    if RUN_PLASMID_SCREEN:
        d["plasmid_cov"] = expand(OUT_DIR + "{sample}/{sample}.plasmid_coverage.tsv", sample=OK)
    if RUN_CASSETTES:
        d["cassette_more"] = expand(OUT_DIR + "{sample}/annotated/{sample}.cassette_junctions.annotated.tsv",
                                    sample=OK) \
                           + expand(OUT_DIR + "{sample}/annotated/{sample}.cassette_hotspots.annotated.tsv",
                                    sample=OK) \
                           + expand(OUT_DIR + "{sample}/cassettes/{sample}.cassette_counts.tsv", sample=OK) \
                           + expand(OUT_DIR + "{sample}/cassettes/{sample}.cassette_depth.tsv", sample=OK)
    if RUN_PICTURES:
        d["pictures"] = expand(OUT_DIR + "{sample}/cassettes/pictures", sample=OK)
    return d


rule html_report:
    input:
        unpack(_report_inputs)
    output:
        html = OUT_DIR + "report.html",
    log:
        OUT_DIR + "logs/html_report.log",
    params:
        samples = samples[["sample_name", "platform", "reads1", "reads2", "reference"]].to_dict("records"),
        out_dir = OUT_DIR,
        toggles = {
            "filtering":     bool(RUN_FILTERING),
            "variants":      bool(RUN_VARIANTS),
            "gatc":          bool(RUN_GATC),
            "coverage":      bool(RUN_COVERAGE),
            "sv":            bool(RUN_SV),
            "multimer":      bool(RUN_MULTIMER and ONT_SAMPLES),
            "denovo":        bool(RUN_DENOVO),
            "plasmid_panel": bool(RUN_PLASMID_SCREEN),
            "cassettes":     bool(RUN_CASSETTES),
            "pictures":      bool(RUN_PICTURES),
        },
        embed_pictures = config.get("report_embed_pictures", True),
        cassettes = CASSETTES if RUN_CASSETTES else [],      # their features under the cassette plots
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/report.py"
