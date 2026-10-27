# ======================================================================
# rules/common.smk :: shared helpers + input staging
# ======================================================================

from Bio import SeqIO


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def ref_length(sample):
    """Total reference length (bp) for a sample, from its staged FASTA."""
    fasta = f"{OUT_DIR}{sample}/{sample}.ref.fasta"
    return sum(len(rec.seq) for rec in SeqIO.parse(fasta, "fasta"))

def barcode_fastq(wildcards):
    """Resolve the demultiplexed FASTQ for a sample's barcode."""
    bc = samples.at[wildcards.sample, "barcode"]
    return f"{FASTQ_DIR}{FASTQ_BASE_NAME}{bc}.fastq.gz"

def sample_ref_gb(wildcards):
    """Resolve the GenBank reference for a sample."""
    return f"{REF_DIR}{samples.at[wildcards.sample, 'ref_fasta']}.gb"


# ----------------------------------------------------------------------
# Stage per-sample inputs
# ----------------------------------------------------------------------
rule stage_reads:
    input:
        barcode_fastq,
    output:
        OUT_DIR + "{sample}/{sample}.fastq.gz",
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
        from Bio import SeqIO
        SeqIO.convert(input.gb, "genbank", output.fasta, "fasta")
