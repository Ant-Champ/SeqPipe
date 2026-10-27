# ======================================================================
# rules/assembly.smk :: reference-free de novo assembly (Flye)
#
# Runs directly off the same filtered reads used by the alignment branch,
# so it can be enabled on its own or alongside reference-based analysis.
# Samples carrying more than one replicon (e.g. chromosome + plasmid) are
# expected to assemble into multiple contigs -- see the per-sample
# assembly_info.txt for the per-contig breakdown (length, coverage,
# circularity), that is not an error.
#
# Flye hard-crashes ("No disjointigs were assembled") when a sample's
# filtered read yield is too low relative to its target size for reads to
# overlap at all -- this pipeline mixes small plasmids and whole genomes
# in the same run, so "too low" only means something relative to each
# sample's own reference length. assemble_denovo checks that before ever
# calling Flye and skips cleanly (empty fasta/gfa, an explanatory
# assembly_info.txt, exit 0) instead of failing the job.
# ======================================================================

def flye_min_bases(wildcards):
    """Minimum filtered-read bases required to attempt assembly (flye_min_coverage x ref length)."""
    return int(ref_length(wildcards.sample) * config.get("flye_min_coverage", 10))


rule assemble_denovo:
    input:
        reads = OUT_DIR + "{sample}/{sample}.filt.fastq.gz",
    output:
        asm_dir = directory(OUT_DIR + "{sample}/assembly"),
        fasta   = OUT_DIR + "{sample}/{sample}.assembly.fasta",
        info    = OUT_DIR + "{sample}/{sample}.assembly_info.txt",
        gfa     = OUT_DIR + "{sample}/{sample}.assembly_graph.gfa",
    log:
        OUT_DIR + "{sample}/logs/assemble_denovo.log",
    threads: config.get("flye_threads", 4)
    resources:
        mem_mb  = config.get("flye_mem_mb", 8000),
        runtime = config.get("flye_runtime", 180),   # minutes (SLURM)
    params:
        read_mode = config.get("flye_read_mode", "--nano-raw"),
        min_bases = flye_min_bases,
    conda:
        "../envs/assembly.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{output.asm_dir}"
        total_bases=$(zcat "{input.reads}" | awk 'NR % 4 == 2 {{ c += length($0) }} END {{ print c+0 }}')
        min_bases={params.min_bases}
        echo "filtered read bases: $total_bases (minimum required: $min_bases)" > "{log}"
        if [ "$total_bases" -lt "$min_bases" ]; then
            echo "SKIPPED: only $total_bases bp of filtered reads, below the $min_bases bp minimum (flye_min_coverage x reference length) -- not enough read depth to assemble." >> "{log}"
            printf '#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n' > "{output.info}"
            printf '# SKIPPED: insufficient read depth (%s bp < %s bp minimum)\n' "$total_bases" "$min_bases" >> "{output.info}"
            : > "{output.fasta}"
            : > "{output.gfa}"
            exit 0
        fi
        flye {params.read_mode} "{input.reads}" \
             --out-dir "{output.asm_dir}" \
             --threads {threads} >> "{log}" 2>&1
        cp -- "{output.asm_dir}/assembly.fasta"     "{output.fasta}"
        cp -- "{output.asm_dir}/assembly_info.txt"  "{output.info}"
        cp -- "{output.asm_dir}/assembly_graph.gfa" "{output.gfa}"
        """
