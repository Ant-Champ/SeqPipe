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


def ok_samples(names=None):
    """The samples (default: all) that passed the depth check (rule depth_check),
    in sheet order. Only call inside input / params functions: it waits for the
    depth_checks checkpoint (every sample's QC) before answering."""
    names = SAMPLE_NAMES if names is None else names
    if not MIN_THEORETICAL_DEPTH:
        return list(names)
    with open(checkpoints.depth_checks.get().output[0]) as fh:
        status = dict(l.split("\t")[0::5] for l in fh.read().splitlines()[1:] if l.strip())
    return [s for s in names if status.get(s) == "ok"]


def depth_checks():
    """depth_check tables ([] when the check is off)."""
    if not MIN_THEORETICAL_DEPTH:
        return []
    return [OUT_DIR + "depth_check.tsv"] + \
        expand(OUT_DIR + "{sample}/qc/{sample}.depth_check.tsv", sample=SAMPLE_NAMES)


def raw_reads(sample):
    """Raw input FASTQ(s) of a sample, as written in the sheet: [reads] for ONT, [R1, R2] for short reads."""
    r = [samples.at[sample, "reads1"]]
    if not is_ont(sample):
        r.append(samples.at[sample, "reads2"])
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
    return _ref_path(samples.at[sample, "reference"])


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
