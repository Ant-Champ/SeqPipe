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


def safe_record_names(records):
    """Give every GenBank record a usable contig name, in place. GenBank files
    exported by SnapGene / Benchling often have 'ACCESSION .' / 'VERSION .', so
    Biopython names the record '.' -- and htslib reads a contig called '.' as
    "the whole file": every BAM lookup on it returns all reads (cassette reads
    counted as genome reads, pictures that never finish). A bad id is replaced
    by the LOCUS name, else contig<N>; names are made unique. Returns
    [(old, new)] for the names that changed."""
    import re
    ok = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*$")
    used, changed = set(), []
    for i, rec in enumerate(records, 1):
        name = rec.id
        if not ok.match(name or ""):
            name = rec.name if ok.match(rec.name or "") else re.sub(r"[^A-Za-z0-9._-]", "_", rec.name or "")
            if not ok.match(name):
                name = f"contig{i}"
        base, k = name, 2
        while name in used:
            name, k = f"{base}_{k}", k + 1
        used.add(name)
        if name != rec.id:
            changed.append((rec.id, name))
        rec.id = rec.name = name
        rec.annotations["accessions"] = [name]
    return changed


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
    """One FASTQ per sample: the sheet's file (symlink, or copy with
    stage_mode: copy), or every part of a split FASTQ merged in numeric order
    (see ONT_PARTS in the Snakefile). Plain .fastq input is gzipped."""
    input:
        lambda wc: ONT_PARTS[wc.sample],
    output:
        OUT_DIR + "{sample}/{sample}.fastq.gz",
    wildcard_constraints:
        sample = ONT_RE,
    log:
        OUT_DIR + "{sample}/logs/stage_reads.log",
    params:
        mode   = config.get("stage_mode", "symlink"),
        nparts = lambda wc: len(ONT_PARTS[wc.sample]),
        gz     = lambda wc: "yes" if ONT_PARTS[wc.sample][0].endswith(".gz") else "no",
    shell:
        r"""
        set -euo pipefail
        if [ {params.nparts} -gt 1 ]; then
            # gzip members can simply be concatenated; plain parts are compressed
            if [ "{params.gz}" = yes ]; then cat -- {input:q} > "{output}"
            else cat -- {input:q} | gzip > "{output}"; fi
            echo "merged {params.nparts} parts:" > "{log}"
            printf '%s
' {input:q} >> "{log}"
        elif [ "{params.gz}" = no ]; then
            gzip -c -- "{input}" > "{output}"
            echo "gzipped {input}" > "{log}"
        elif [ "{params.mode}" = "copy" ]; then
            cp -- "{input}" "{output}"
            echo "copied {input}" > "{log}"
        else
            ln -sf -- "$(readlink -f "{input}")" "{output}"
            echo "linked {input}" > "{log}"
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
    """The sample's GenBank reference, with usable contig names (see
    safe_record_names): every later step takes its contig names from here."""
    input:
        sample_ref_gb,
    output:
        OUT_DIR + "{sample}/{sample}.ref.gb",
    log:
        OUT_DIR + "{sample}/logs/stage_ref.log",
    run:
        import warnings
        from Bio import SeqIO, BiopythonWarning
        records = list(SeqIO.parse(input[0], "genbank"))
        changed = safe_record_names(records)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)     # long LOCUS names
            SeqIO.write(records, output[0], "genbank")
        with open(log[0], "w") as fh:
            fh.write(f"{len(records)} record(s): {', '.join(r.id for r in records)}\n")
            for old, new in changed:
                fh.write(f"renamed contig '{old}' -> '{new}' (unusable name in {input[0]})\n")


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
