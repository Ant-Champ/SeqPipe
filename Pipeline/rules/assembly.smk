# ======================================================================
# rules/assembly.smk :: reference-free de novo assembly
#   ONT          : Flye
#   short reads  : Unicycler (short-read mode, SPAdes-based)
#
# Runs directly off the same filtered reads used by the alignment branch,
# so it can be enabled on its own or alongside reference-based analysis.
# Both rules write the same outputs (assembly.fasta, assembly_info.txt,
# assembly_graph.gfa), so downstream steps do not care how a sample was
# assembled. Samples carrying more than one replicon (e.g. chromosome +
# plasmid) are expected to assemble into multiple contigs -- see the
# per-sample assembly_info.txt for the per-contig breakdown (length,
# coverage, circularity), that is not an error.
#
# Assemblers hard-crash (Flye: "No disjointigs were assembled") when a
# sample's filtered read yield is too low relative to its target size --
# this pipeline mixes small plasmids and whole genomes in the same run,
# so "too low" only means something relative to each sample's own
# reference length. Both rules check that before calling the assembler
# and skip cleanly (empty fasta/gfa, an explanatory assembly_info.txt
# marked '# SKIPPED', exit 0) instead of failing the job.
#
# Unicycler writes no assembly_info.txt: it is built from the contig
# headers of assembly.fasta (">1 length=4641652 depth=1.00x circular=true")
# with the columns  #seq_name  length  cov.  circ.
# ======================================================================

def flye_min_bases(wildcards):
    """Minimum filtered-read bases required to attempt assembly (flye_min_coverage x ref length)."""
    return int(ref_length(wildcards.sample) * config.get("flye_min_coverage", 10))


def unicycler_min_bases(wildcards):
    """Minimum bases (R1 + R2) required to attempt assembly (unicycler_min_coverage x ref length).

    unicycler_min_coverage falls back to flye_min_coverage, then 10.
    """
    cov = config.get("unicycler_min_coverage", config.get("flye_min_coverage", 10))
    return int(ref_length(wildcards.sample) * cov)


rule assemble_ont:
    input:
        reads = OUT_DIR + "{sample}/{sample}.filt.fastq.gz",
    output:
        asm_dir = directory(OUT_DIR + "{sample}/assembly"),
        fasta   = OUT_DIR + "{sample}/{sample}.assembly.fasta",
        info    = OUT_DIR + "{sample}/{sample}.assembly_info.txt",
        gfa     = OUT_DIR + "{sample}/{sample}.assembly_graph.gfa",
    wildcard_constraints:
        sample = ONT_RE,
    log:
        OUT_DIR + "{sample}/logs/assemble.log",
    threads: config.get("flye_threads", 4)
    resources:
        mem_mb = config.get("flye_mem_mb", 8000),
    params:
        read_mode = config.get("flye_read_mode", "--nano-raw"),
        min_bases = flye_min_bases,
    conda:
        "../envs/flye.yaml"
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


rule assemble_sr:
    input:
        r1 = OUT_DIR + "{sample}/{sample}.filt_R1.fastq.gz",
        r2 = OUT_DIR + "{sample}/{sample}.filt_R2.fastq.gz",
    output:
        asm_dir = directory(OUT_DIR + "{sample}/assembly"),
        fasta   = OUT_DIR + "{sample}/{sample}.assembly.fasta",
        info    = OUT_DIR + "{sample}/{sample}.assembly_info.txt",
        gfa     = OUT_DIR + "{sample}/{sample}.assembly_graph.gfa",
    wildcard_constraints:
        sample = SR_RE,
    log:
        OUT_DIR + "{sample}/logs/assemble.log",
    threads: config.get("unicycler_threads", 8)
    resources:
        mem_mb = config.get("unicycler_mem_mb", 16000),
    params:
        mode      = config.get("unicycler_mode", "normal"),
        min_bases = unicycler_min_bases,
    conda:
        "../envs/unicycler.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{output.asm_dir}"
        # bases in both mates (FASTQ sequence lines)
        total_bases=$(zcat "{input.r1}" "{input.r2}" | awk 'NR % 4 == 2 {{ c += length($0) }} END {{ print c+0 }}')
        min_bases={params.min_bases}
        echo "filtered read bases (R1+R2): $total_bases (minimum required: $min_bases)" > "{log}"
        if [ "$total_bases" -lt "$min_bases" ]; then
            echo "SKIPPED: only $total_bases bp of filtered reads, below the $min_bases bp minimum (unicycler_min_coverage x reference length) -- not enough read depth to assemble." >> "{log}"
            printf '#seq_name\tlength\tcov.\tcirc.\n' > "{output.info}"
            printf '# SKIPPED: insufficient read depth (%s bp < %s bp minimum)\n' "$total_bases" "$min_bases" >> "{output.info}"
            : > "{output.fasta}"
            : > "{output.gfa}"
            exit 0
        fi
        unicycler -1 "{input.r1}" -2 "{input.r2}" \
                  -o "{output.asm_dir}" \
                  -t {threads} \
                  --mode {params.mode} >> "{log}" 2>&1
        cp -- "{output.asm_dir}/assembly.fasta" "{output.fasta}"
        cp -- "{output.asm_dir}/assembly.gfa"   "{output.gfa}"

        # Flye-like per-contig table from the Unicycler FASTA headers:
        #   >1 length=4641652 depth=1.00x circular=true
        awk 'BEGIN {{ print "#seq_name\tlength\tcov.\tcirc." }}
             /^>/ {{
                 name = substr($1, 2); len = "NA"; cov = "NA"; circ = "N"
                 for (i = 2; i <= NF; i++) {{
                     if ($i ~ /^length=/)        {{ len = substr($i, 8) }}
                     else if ($i ~ /^depth=/)    {{ cov = substr($i, 7); sub(/x$/, "", cov) }}
                     else if ($i ~ /^circular=/) {{ circ = ($i == "circular=true") ? "Y" : "N" }}
                 }}
                 printf "%s\t%s\t%s\t%s\n", name, len, cov, circ
             }}' "{output.fasta}" > "{output.info}"
        """
