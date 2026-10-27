# ======================================================================
# rules/common.smk :: shared helpers + input staging
# ======================================================================

import functools
from Bio import SeqIO


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def is_ont(sample):
    return samples.at[sample, "platform"] == "ont"


@functools.lru_cache(maxsize=None)
def _genbank_length(path):
    return sum(len(rec.seq) for rec in SeqIO.parse(path, "genbank"))


def ref_length(sample):
    """Total reference length (bp, all records) for a sample. Read from the
    source GenBank file (not the staged FASTA) because params functions are
    evaluated while the DAG is built, before staging has run."""
    return _genbank_length(ref_gb_path(sample))


def _resolve_fastq(value):
    """A sheet entry is either a path (absolute, or relative to fastq_dir)
    or, for ONT, a barcode -> {fastq_dir}{fastq_base_name}{barcode}.fastq.gz"""
    if re.search(r"\.f(ast)?q(\.gz)?$", value):
        return value if os.path.isabs(value) else FASTQ_DIR + value
    return f"{FASTQ_DIR}{FASTQ_BASE_NAME}{value}.fastq.gz"


def raw_reads(sample):
    """Raw input FASTQ(s) of a sample: [reads] for ONT, [R1, R2] for short reads."""
    r = [_resolve_fastq(samples.at[sample, "reads1"])]
    if not is_ont(sample):
        r.append(_resolve_fastq(samples.at[sample, "reads2"]))
    return r


def staged_reads(sample):
    if is_ont(sample):
        return [OUT_DIR + f"{sample}/{sample}.fastq.gz"]
    return [OUT_DIR + f"{sample}/{sample}_R1.fastq.gz", OUT_DIR + f"{sample}/{sample}_R2.fastq.gz"]


def filt_reads(wildcards):
    """Filtered reads used by every analysis: [reads] (ONT) or [R1, R2] (short reads)."""
    s = wildcards.sample
    if is_ont(s):
        return [OUT_DIR + f"{s}/{s}.filt.fastq.gz"]
    return [OUT_DIR + f"{s}/{s}.filt_R1.fastq.gz", OUT_DIR + f"{s}/{s}.filt_R2.fastq.gz"]


def minimap2_preset(wildcards):
    """minimap2 -x preset for this sample's read type."""
    return "map-ont" if is_ont(wildcards.sample) else "sr"


def ref_gb_path(sample):
    """Resolve the GenBank reference for a sample (name with or without .gb)."""
    ref = samples.at[sample, "reference"]
    if not re.search(r"\.(gb|gbk|genbank)$", ref):
        ref += ".gb"
    return ref if os.path.isabs(ref) else REF_DIR + ref


def sample_ref_gb(wildcards):
    return ref_gb_path(wildcards.sample)


# ----------------------------------------------------------------------
# Stage per-sample inputs
# ----------------------------------------------------------------------
rule stage_reads_ont:
    input:
        lambda wc: raw_reads(wc.sample),
    output:
        OUT_DIR + "{sample}/{sample}.fastq.gz",
    wildcard_constraints:
        sample = ONT_RE,
    params:
        mode = config.get("stage_mode", "symlink"),
    shell:
        r"""
        if [ "{params.mode}" = "copy" ]; then
            cp -- "{input}" "{output}"
        else
            ln -sf -- "$(readlink -f "{input}")" "{output}"
        fi
        """


rule stage_reads_sr:
    input:
        r1 = lambda wc: raw_reads(wc.sample)[0],
        r2 = lambda wc: raw_reads(wc.sample)[1],
    output:
        r1 = OUT_DIR + "{sample}/{sample}_R1.fastq.gz",
        r2 = OUT_DIR + "{sample}/{sample}_R2.fastq.gz",
    wildcard_constraints:
        sample = SR_RE,
    params:
        mode = config.get("stage_mode", "symlink"),
    shell:
        r"""
        for pair in "{input.r1}|{output.r1}" "{input.r2}|{output.r2}"; do
            src="${{pair%%|*}}"; dst="${{pair##*|}}"
            if [ "{params.mode}" = "copy" ]; then
                cp -- "$src" "$dst"
            else
                ln -sf -- "$(readlink -f "$src")" "$dst"
            fi
        done
        """


rule stage_ref:
    input:
        sample_ref_gb,
    output:
        OUT_DIR + "{sample}/{sample}.ref.gb",
    shell:
        r"""cp -- "{input}" "{output}" """


rule genbank_to_fasta:
    input:
        gb = OUT_DIR + "{sample}/{sample}.ref.gb",
    output:
        fasta = OUT_DIR + "{sample}/{sample}.ref.fasta",
    run:
        # Contig names = GenBank record.id (e.g. NC_000913.3): every table
        # is annotated by matching on it.
        from Bio import SeqIO
        SeqIO.convert(input.gb, "genbank", output.fasta, "fasta")


rule faidx_reference:
    input:
        OUT_DIR + "{sample}/{sample}.ref.fasta",
    output:
        OUT_DIR + "{sample}/{sample}.ref.fasta.fai",
    conda:
        "../envs/align.yaml"
    shell:
        r"""samtools faidx "{input}" """
