# ======================================================================
# rules/plasmid_screen.smk :: plasmid-panel presence/coverage screening
# (ONT and short reads)
#   + cassette read extraction, reference alignment, hotspot detection
#     (run_cassette_hotspots) and automated IGV snapshots (run_igv_snapshots)
#
# Every filtered read set is aligned to the combined plasmid panel
# (minimap2 map-ont for ONT, sr for paired short reads). Reads hitting the
# cassette record of the panel (config cassette_reference) are pulled out
# -- for short reads BOTH mates of every hit name -- and re-aligned to the
# sample's own reference; clusters of such alignments are the "hotspots"
# (e.g. insertion sites), optionally snapshotted in IGV.
#
# Old config keys chlor_cassette_reference / chlor_hotspot_* are still
# accepted as fall-backs for the cassette_* names.
# ======================================================================

CASSETTE_REF = config.get("cassette_reference", config.get("chlor_cassette_reference", "Chlor_cassette"))


rule build_plasmid_panel:
    """Concatenate every fasta in plasmid_panel_dir into one combined reference."""
    output:
        panel = OUT_DIR + "plasmid_panel/combined_panel.fasta",
    log:
        OUT_DIR + "logs/build_plasmid_panel.log",
    params:
        panel_dir = PLASMID_PANEL_DIR,
    run:
        import glob
        import os

        files = sorted(
            glob.glob(os.path.join(params.panel_dir, "*.fasta"))
            + glob.glob(os.path.join(params.panel_dir, "*.fa"))
        )

        if not files:
            raise WorkflowError(
                f"no .fasta/.fa files found in plasmid_panel_dir ({params.panel_dir})"
            )

        os.makedirs(os.path.dirname(output.panel), exist_ok=True)

        with open(output.panel, "w") as out_fh, open(log[0], "w") as log_fh:
            log_fh.write("plasmid panel built from:\n")
            for f in files:
                log_fh.write(f"  {f}\n")
                with open(f) as in_fh:
                    out_fh.write(in_fh.read())


rule plasmid_screen_align:
    """Align the filtered reads (ONT, or R1+R2) to the combined plasmid panel."""
    input:
        reads = filt_reads,
        panel = OUT_DIR + "plasmid_panel/combined_panel.fasta",
    output:
        bam = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam",
        bai = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam.bai",
    log:
        OUT_DIR + "{sample}/logs/plasmid_screen_align.log",
    threads: config.get("plasmid_screen_threads", 4)
    resources:
        mem_mb = config.get("plasmid_screen_mem_mb", 4000),
    params:
        preset = minimap2_preset,
    conda:
        "../envs/plasmid_screen.yaml"
    shell:
        r"""
        set -euo pipefail

        minimap2 -t {threads} -ax {params.preset} \
            "{input.panel}" {input.reads} 2> "{log}" \
        | samtools sort -@ {threads} -o "{output.bam}" - 2>> "{log}"

        samtools index "{output.bam}" 2>> "{log}"
        """


rule plasmid_screen_coverage:
    """Per-plasmid breadth/depth for one sample (one row per panel plasmid)."""
    input:
        bam = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam",
    output:
        tsv = OUT_DIR + "{sample}/{sample}.plasmid_coverage.tsv",
    log:
        OUT_DIR + "{sample}/logs/plasmid_screen_coverage.log",
    conda:
        "../envs/plasmid_screen.yaml"
    shell:
        r"""samtools coverage "{input.bam}" > "{output.tsv}" 2> "{log}" """


rule aggregate_plasmid_coverage:
    """Run-level long-format table: one row per sample x plasmid."""
    input:
        tsvs = expand(
            OUT_DIR + "{sample}/{sample}.plasmid_coverage.tsv",
            sample=SAMPLE_NAMES
        ),
    output:
        summary = OUT_DIR + "plasmid_coverage_summary.tsv",
    log:
        OUT_DIR + "logs/aggregate_plasmid_coverage.log",
    params:
        exe = script_path("coverage_summary", "aggregate_coverage.1.py"),
        samples = SAMPLE_NAMES,
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.exe} -in {input.tsvs} -samples {params.samples} "
        "-name_col replicon -source plasmid_panel -out {output.summary} > {log} 2>&1"


# ----------------------------------------------------------------------
# Cassette hotspots (run_cassette_hotspots)
# ----------------------------------------------------------------------
if RUN_HOTSPOTS:

    def cassette_reads(wildcards):
        """Extracted cassette reads: [fastq] (ONT) or [R1, R2] (short reads)."""
        s = wildcards.sample
        d = f"{OUT_DIR}{s}/plasmidic_reads/{s}"
        if is_ont(s):
            return [d + ".cassette.fastq.gz"]
        return [d + ".cassette_R1.fastq.gz", d + ".cassette_R2.fastq.gz"]


    rule extract_cassette_reads_ont:
        """
        Extract the unique original FASTQ reads having an alignment to the
        cassette. Read names are collected from the panel BAM, then the
        FASTQ is filtered so supplementary/secondary alignments do not
        duplicate a read.
        """
        input:
            bam   = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam",
            bai   = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam.bai",
            reads = filt_reads,
        output:
            reads = OUT_DIR + "{sample}/plasmidic_reads/{sample}.cassette.fastq.gz",
        wildcard_constraints:
            sample = ONT_RE,
        log:
            OUT_DIR + "{sample}/logs/extract_cassette_reads.log",
        params:
            exe = script_path("extract_cassette_reads", "extract_cassette_reads.py"),
            cassette_ref = CASSETTE_REF,
        conda:
            "../envs/plasmid_screen.yaml"
        shell:
            r"""
            set -euo pipefail

            mkdir -p "$(dirname "{output.reads}")"

            # Alignments to the cassette record only; the helper collects
            # the unique query names and filters the FASTQ.
            samtools view -b "{input.bam}" "{params.cassette_ref}" \
                2> "{log}" > "{output.reads}.tmp.bam"

            python "{params.exe}" \
                --bam "{output.reads}.tmp.bam" \
                --reads {input.reads} \
                --out "{output.reads}" \
                >> "{log}" 2>&1

            rm -f "{output.reads}.tmp.bam"
            """


    rule extract_cassette_reads_sr:
        """
        Extract both mates (R1 + R2, kept in sync) of every read name with an
        alignment to the cassette; a pair is kept when either mate hit it.
        """
        input:
            bam   = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam",
            bai   = OUT_DIR + "{sample}/{sample}.plasmid_panel.sorted.bam.bai",
            reads = filt_reads,
        output:
            r1 = OUT_DIR + "{sample}/plasmidic_reads/{sample}.cassette_R1.fastq.gz",
            r2 = OUT_DIR + "{sample}/plasmidic_reads/{sample}.cassette_R2.fastq.gz",
        wildcard_constraints:
            sample = SR_RE,
        log:
            OUT_DIR + "{sample}/logs/extract_cassette_reads.log",
        params:
            exe = script_path("extract_cassette_reads", "extract_cassette_reads.py"),
            cassette_ref = CASSETTE_REF,
        conda:
            "../envs/plasmid_screen.yaml"
        shell:
            r"""
            set -euo pipefail

            mkdir -p "$(dirname "{output.r1}")"

            samtools view -b "{input.bam}" "{params.cassette_ref}" \
                2> "{log}" > "{output.r1}.tmp.bam"

            python "{params.exe}" \
                --bam "{output.r1}.tmp.bam" \
                --reads {input.reads} \
                --out "{output.r1}" "{output.r2}" \
                >> "{log}" 2>&1

            rm -f "{output.r1}.tmp.bam"
            """


    rule align_cassette_reads_to_reference:
        """
        Align the cassette-derived reads against the sample reference genome.
        The resulting BAM is intended for IGV and hotspot detection.
        """
        input:
            reads   = cassette_reads,
            ref     = OUT_DIR + "{sample}/{sample}.ref.fasta",
            ref_fai = OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
        output:
            bam = OUT_DIR + "{sample}/plasmidic_reads/{sample}.cassette.sorted.bam",
            bai = OUT_DIR + "{sample}/plasmidic_reads/{sample}.cassette.sorted.bam.bai",
        log:
            OUT_DIR + "{sample}/logs/align_cassette_to_reference.log",
        threads: config.get("plasmid_screen_threads", 4)
        resources:
            mem_mb = config.get("plasmid_screen_mem_mb", 4000),
        params:
            preset = minimap2_preset,
        conda:
            "../envs/plasmid_screen.yaml"
        shell:
            r"""
            set -euo pipefail

            mkdir -p "$(dirname "{output.bam}")"

            minimap2 -t {threads} -ax {params.preset} \
                "{input.ref}" {input.reads} 2> "{log}" \
            | samtools sort -@ {threads} -o "{output.bam}" - 2>> "{log}"

            samtools index "{output.bam}" 2>> "{log}"
            """


    rule cassette_hotspots:
        """
        Find clusters of cassette-read alignments on the reference genome.
        Each cluster is expanded by cassette_hotspot_flank bp for IGV
        visualization.
        """
        input:
            bam = OUT_DIR + "{sample}/plasmidic_reads/{sample}.cassette.sorted.bam",
        output:
            tsv = OUT_DIR + "{sample}/plasmidic_reads/{sample}.cassette_hotspots.tsv",
            bed = OUT_DIR + "{sample}/plasmidic_reads/{sample}.cassette_hotspots.bed",
        log:
            OUT_DIR + "{sample}/logs/cassette_hotspots.log",
        params:
            exe = script_path("find_cassette_hotspots", "find_cassette_hotspots.py"),
            label = CASSETTE_REF,
            max_gap = config.get("cassette_hotspot_max_gap", config.get("chlor_hotspot_max_gap", 1000)),
            flank = config.get("cassette_hotspot_flank", config.get("chlor_hotspot_flank", 500)),
            min_reads = config.get("cassette_hotspot_min_reads", config.get("chlor_hotspot_min_reads", 2)),
        conda:
            "../envs/plasmid_screen.yaml"
        shell:
            r"""
            set -euo pipefail

            python "{params.exe}" \
                --bam "{input.bam}" \
                --out-tsv "{output.tsv}" \
                --out-bed "{output.bed}" \
                --max-gap {params.max_gap} \
                --flank {params.flank} \
                --min-reads {params.min_reads} \
                --label "{params.label}" \
                > "{log}" 2>&1
            """


# ----------------------------------------------------------------------
# IGV snapshots (run_igv_snapshots; IGV itself is not in conda -- the
# rule uses config igv_exe and xvfb-run when it is on PATH)
# ----------------------------------------------------------------------
if RUN_IGV:

    rule cassette_igv_snapshots:
        """
        Generate one IGV PNG per cassette hotspot.

        The snapshot contains:
          - reference sequence
          - GenBank-derived annotations
          - all-read BAM
          - cassette-read BAM
          - hotspot +/- cassette_hotspot_flank bp
        """
        input:
            hotspots     = OUT_DIR + "{sample}/plasmidic_reads/{sample}.cassette_hotspots.tsv",
            all_bam      = OUT_DIR + "{sample}/{sample}.sorted.bam",
            all_bai      = OUT_DIR + "{sample}/{sample}.sorted.bam.bai",
            cassette_bam = OUT_DIR + "{sample}/plasmidic_reads/{sample}.cassette.sorted.bam",
            cassette_bai = OUT_DIR + "{sample}/plasmidic_reads/{sample}.cassette.sorted.bam.bai",
            ref          = OUT_DIR + "{sample}/{sample}.ref.fasta",
            ref_fai      = OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
            gb           = OUT_DIR + "{sample}/{sample}.ref.gb",
        output:
            snapshots = directory(
                OUT_DIR + "{sample}/plasmidic_reads/igv_cassette_hotspots"
            ),
        log:
            OUT_DIR + "{sample}/logs/cassette_igv_snapshots.log",
        params:
            exe = script_path("run_igv_snapshots", "run_igv_snapshots.py"),
            igv = config.get("igv_exe", "igv.sh"),
            label = CASSETTE_REF,
            heap_mb = config.get("igv_heap_mb", 4000),
        resources:
            mem_mb = config.get("igv_mem_mb", 6000),
        conda:
            "../envs/plasmid_screen.yaml"
        shell:
            r"""
            set -euo pipefail

            mkdir -p "{output.snapshots}"

            python "{params.exe}" \
                --hotspots "{input.hotspots}" \
                --all-bam "{input.all_bam}" \
                --cassette-bam "{input.cassette_bam}" \
                --reference "{input.ref}" \
                --genbank "{input.gb}" \
                --outdir "{output.snapshots}" \
                --igv "{params.igv}" \
                --label "{params.label}" \
                --memory {params.heap_mb} \
                > "{log}" 2>&1
            """
