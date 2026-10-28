# ======================================================================
# rules/cassettes.smk :: where are the cassettes inserted in the genome?
# (run_cassette_search; ONT and short reads; every cassette in cassette_dir)
#
# Each sample's reads are mapped once against its genome + ALL cassettes
# (one combined reference), so that a read matching the region two
# cassettes share (Keio FRT / priming-site ends) gets MAPQ 0 instead of
# being given to one cassette at random. Reads that are part genome / part
# cassette give each junction to the base; a left + right junction pair is
# one insertion site (for a Keio replacement, the span is the deleted gene).
# Short reads: soft-clipped genome reads whose clip continues into a cassette
# are junction evidence too (most junction reads have no split alignment).
# The cassette IDENTITY of a site is only called from uniquely mapped
# evidence (the read itself, or its mate); otherwise it is reported as
# ambiguous(A|B) -- the position is still valid.
# ======================================================================

CASSETTE_SEARCH_PARAMS = {
    "ont":   {"min_segment": 100, "cluster_tol": 20},
    "short": {"min_segment": 30,  "cluster_tol": 10},
}


def _cassette_param(key):
    def _get(wildcards):
        kind = "ont" if is_ont(wildcards.sample) else "short"
        return config.get(f"cassette_{key}_{kind}", CASSETTE_SEARCH_PARAMS[kind][key])
    return _get


rule build_cassettes:
    """All cassettes of cassette_dir -> one FASTA (record names = cassette names)."""
    input:
        [c["path"] for c in CASSETTES],
    output:
        fasta = OUT_DIR + "cassettes/cassettes.fasta",
    log:
        OUT_DIR + "logs/build_cassettes.log",
    run:
        from Bio import SeqIO
        with open(output.fasta, "w") as out_fh, open(log[0], "w") as log_fh:
            for c in CASSETTES:
                rec = list(SeqIO.parse(c["path"], c["format"]))[c["index"]]
                out_fh.write(f">{c['name']}\n{str(rec.seq).upper()}\n")
                log_fh.write(f"{c['name']}\t{len(rec.seq)} bp\t{c['path']}\n")


rule cassette_reference:
    """Sample genome + all cassettes, as one reference."""
    input:
        genome    = OUT_DIR + "{sample}/{sample}.ref.fasta",
        cassettes = OUT_DIR + "cassettes/cassettes.fasta",
    output:
        OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.fasta",
    run:
        from Bio import SeqIO
        genome_ids = {r.id for r in SeqIO.parse(input.genome, "fasta")}
        clash = genome_ids & {c["name"] for c in CASSETTES}
        if clash:
            raise WorkflowError(f"cassette name(s) {sorted(clash)} also used by a genome record: rename the cassette file(s)")
        with open(output[0], "w") as out_fh:
            for path in (input.genome, input.cassettes):
                with open(path) as in_fh:
                    out_fh.write(in_fh.read())


rule cassette_align_ont:
    input:
        reads = OUT_DIR + "{sample}/{sample}.filt.fastq.gz",
        ref   = OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.fasta",
    output:
        bam = OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.sorted.bam",
        bai = OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.sorted.bam.bai",
    wildcard_constraints:
        sample = ONT_RE,
    log:
        OUT_DIR + "{sample}/logs/cassette_align.log",
    threads: config.get("align_threads", 8)
    resources:
        mem_mb = config.get("align_mem_mb", 16000),
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        set -euo pipefail
        tmp="${{TMPDIR:-/tmp}}"
        minimap2 -ax map-ont -t {threads} "{input.ref}" "{input.reads}" 2> "{log}" \
          | samtools sort -@ {threads} -T "$tmp/{wildcards.sample}.cassette" -o "{output.bam}" - 2>> "{log}"
        samtools index "{output.bam}"
        """


rule cassette_bwa_index:
    input:
        OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.fasta",
    output:
        multiext(OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    wildcard_constraints:
        sample = SR_RE,
    log:
        OUT_DIR + "{sample}/logs/cassette_bwa_index.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""bwa-mem2 index "{input}" > "{log}" 2>&1"""


rule cassette_align_sr:
    input:
        r1    = OUT_DIR + "{sample}/{sample}.filt_R1.fastq.gz",
        r2    = OUT_DIR + "{sample}/{sample}.filt_R2.fastq.gz",
        ref   = OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.fasta",
        index = multiext(OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        bam = OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.sorted.bam",
        bai = OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.sorted.bam.bai",
    wildcard_constraints:
        sample = SR_RE,
    log:
        OUT_DIR + "{sample}/logs/cassette_align.log",
    threads: config.get("align_threads", 8)
    resources:
        mem_mb = config.get("align_mem_mb", 16000),
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        set -euo pipefail
        tmp="${{TMPDIR:-/tmp}}"
        bwa-mem2 mem -t {threads} "{input.ref}" "{input.r1}" "{input.r2}" 2> "{log}" \
          | samtools sort -@ {threads} -T "$tmp/{wildcards.sample}.cassette" -o "{output.bam}" - 2>> "{log}"
        samtools index "{output.bam}"
        """


rule find_cassette_insertions:
    input:
        bam = OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.sorted.bam",
        bai = OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.sorted.bam.bai",
        ref = OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.fasta",
    output:
        junctions  = OUT_DIR + "{sample}/cassettes/{sample}.cassette_junctions.tsv",
        insertions = OUT_DIR + "{sample}/cassettes/{sample}.cassette_insertions.tsv",
        counts     = OUT_DIR + "{sample}/cassettes/{sample}.cassette_counts.tsv",
        depth      = OUT_DIR + "{sample}/cassettes/{sample}.cassette_depth.tsv",
    log:
        OUT_DIR + "{sample}/logs/find_cassette_insertions.log",
    params:
        exe         = script_path("find_cassette_insertions", "find_cassette_insertions.py"),
        cassettes   = " ".join(c["name"] for c in CASSETTES),
        platform    = lambda wc: "ont" if is_ont(wc.sample) else "short",
        min_mapq    = config.get("cassette_min_mapq", 20),
        min_reads   = config.get("cassette_min_reads", 2),
        max_span    = config.get("cassette_max_span", 20000),
        min_segment = _cassette_param("min_segment"),
        cluster_tol = _cassette_param("cluster_tol"),
        min_clip    = config.get("cassette_min_clip", 12),
        mate_window = config.get("cassette_mate_window", 1000),
    conda:
        "../envs/regions.yaml"
    shell:
        "python {params.exe} --bam {input.bam} --cassettes {params.cassettes} "
        "--platform {params.platform} --min-mapq {params.min_mapq} "
        "--min-segment {params.min_segment} --cluster-tol {params.cluster_tol} "
        "--min-reads {params.min_reads} --max-span {params.max_span} "
        "--out-junctions {output.junctions} --out-insertions {output.insertions} "
        "--out-counts {output.counts} --out-depth {output.depth} --fasta {input.ref} "
        "--min-clip {params.min_clip} --mate-window {params.mate_window} > {log} 2>&1"


rule annotate_cassette_tables:
    """Genes at / between the junctions (e.g. the gene a Keio cassette replaced)."""
    input:
        tsv = OUT_DIR + "{sample}/cassettes/{sample}.cassette_{table}.tsv",
        gb  = OUT_DIR + "{sample}/{sample}.ref.gb",
    output:
        OUT_DIR + "{sample}/annotated/{sample}.cassette_{table}.annotated.tsv",
    wildcard_constraints:
        table = "insertions|junctions",
    log:
        OUT_DIR + "{sample}/logs/annotate_cassette_{table}.log",
    params:
        exe        = script_path("annotate_regions", "annotate_regions.py"),
        skip_types = config.get("annotate_skip_types", "source,gene"),
    conda:
        "../envs/regions.yaml"
    shell:
        "python {params.exe} --genbank {input.gb} --in {input.tsv} --out {output} "
        "--skip-types {params.skip_types} > {log} 2>&1"


# ----------------------------------------------------------------------
# IGV snapshots (run_igv_snapshots; IGV itself is not in conda -- the
# rule uses config igv_exe and xvfb-run when it is on PATH)
# ----------------------------------------------------------------------
if RUN_IGV:

    rule cassette_igv_snapshots:
        """One IGV PNG per insertion site: reference + annotation, all reads
        vs genome, and the genome+cassettes alignment (junction reads)."""
        input:
            regions      = OUT_DIR + "{sample}/annotated/{sample}.cassette_insertions.annotated.tsv",
            all_bam      = OUT_DIR + "{sample}/{sample}.sorted.bam",
            all_bai      = OUT_DIR + "{sample}/{sample}.sorted.bam.bai",
            cassette_bam = OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.sorted.bam",
            cassette_bai = OUT_DIR + "{sample}/cassettes/{sample}.genome_cassettes.sorted.bam.bai",
            ref          = OUT_DIR + "{sample}/{sample}.ref.fasta",
            ref_fai      = OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
            gb           = OUT_DIR + "{sample}/{sample}.ref.gb",
        output:
            snapshots = directory(OUT_DIR + "{sample}/cassettes/igv_snapshots"),
        log:
            OUT_DIR + "{sample}/logs/cassette_igv_snapshots.log",
        params:
            exe     = script_path("run_igv_snapshots", "run_igv_snapshots.py"),
            igv     = config.get("igv_exe", "igv.sh"),
            flank   = config.get("igv_flank", 1000),
            heap_mb = config.get("igv_heap_mb", 4000),
        resources:
            mem_mb   = config.get("igv_mem_mb", 6000),
            igv_jobs = 1,          # run.sh allows 4 IGVs at once (IGV_JOBS=n to change)
        conda:
            "../envs/plasmid_screen.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p "{output.snapshots}"
            python "{params.exe}" \
                --regions "{input.regions}" \
                --flank {params.flank} \
                --all-bam "{input.all_bam}" \
                --cassette-bam "{input.cassette_bam}" \
                --reference "{input.ref}" \
                --genbank "{input.gb}" \
                --outdir "{output.snapshots}" \
                --igv "{params.igv}" \
                --label cassette \
                --memory {params.heap_mb} \
                > "{log}" 2>&1
            """
