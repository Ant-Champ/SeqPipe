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

    # Read reference fasta file
    record = next(SeqIO.parse(args.ref_fasta, "fasta"))
    ref_sequence = str(record.seq.upper())

    # Read variants file into a DataFrame
    df = pd.read_csv(args.variants, sep="\t")

    # Apply the function to each position in the DataFrame
    df["NearGATC"] = df["POS"].apply(lambda pos: "Yes" if is_near_gatc(ref_sequence, pos-1) else "No")

    # Write the modified DataFrame to the output file
    df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
