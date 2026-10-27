import argparse
import pandas as pd
from Bio import SeqIO

def is_near_gatc(seq, pos, window=10):
    """
    Check if a position is near a GATC site within a given window.
    """
    start = max(0, pos - window)
    end = pos + window
    segment = seq[start:end]

    return "GATC" in segment

def main():
    parser = argparse.ArgumentParser(description="Annotate mutations near GATC sites.")
    parser.add_argument("-ref_fasta", required=True, help="Path to reference FASTA file.")
    parser.add_argument("-variants", required=True, help="Path to variants file.")
    parser.add_argument("-out", required=True, help="Path to output file.")

    args = parser.parse_args()

    # Read all reference records (id -> upper-case sequence)
    seqs = {rec.id: str(rec.seq.upper()) for rec in SeqIO.parse(args.ref_fasta, "fasta")}
    if not seqs:
        raise SystemExit(f"No sequences found in {args.ref_fasta}")
    single = next(iter(seqs.values())) if len(seqs) == 1 else None

    # Read variants file into a DataFrame (CHROM kept as text; may be header-only)
    df = pd.read_csv(args.variants, sep="\t", dtype={"CHROM": str})

    def near(chrom, pos):
        # the variant's own record; single-record FASTA is the only fallback
        seq = seqs.get(chrom, single) if isinstance(chrom, str) else single
        if seq is None:
            return "NA"
        return "Yes" if is_near_gatc(seq, int(pos) - 1) else "No"

    chroms = df["CHROM"] if "CHROM" in df.columns else [None] * len(df)
    df["NearGATC"] = [near(c, p) for c, p in zip(chroms, df["POS"])]

    df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
