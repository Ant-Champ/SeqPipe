#!/usr/bin/env python3
"""
Convert a GenBank file to a simple GFF3 annotation track for IGV.

The reference sequence itself remains the FASTA. This script exports
GenBank features (CDS, genes, tRNA, rRNA, misc_feature, etc.) as GFF3.
"""
import argparse
import re
from Bio import SeqIO


def esc(value):
    value = str(value)
    return value.replace("%", "%25").replace(";", "%3B").replace("=", "%3D").replace(",", "%2C")


def feature_name(feature):
    q = feature.qualifiers
    for key in ("gene", "locus_tag", "label", "name", "product"):
        if key in q and q[key]:
            return str(q[key][0])
    return feature.type


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--genbank", required=True)
    p.add_argument("--gff", required=True)
    args = p.parse_args()

    with open(args.gff, "w") as out:
        out.write("##gff-version 3\n")

        for record in SeqIO.parse(args.genbank, "genbank"):
            seqid = record.id

            for i, feature in enumerate(record.features, start=1):
                if feature.type == "source":
                    continue

                if feature.location is None:
                    continue

                # Handle simple and compound locations.
                parts = list(feature.location.parts)
                for part_index, part in enumerate(parts, start=1):
                    start = int(part.start) + 1
                    end = int(part.end)

                    strand = "."
                    if feature.strand == 1:
                        strand = "+"
                    elif feature.strand == -1:
                        strand = "-"

                    name = feature_name(feature)
                    fid = f"{feature.type}_{i}_{part_index}"

                    attrs = [
                        f"ID={esc(fid)}",
                        f"Name={esc(name)}",
                        f"gene_biotype={esc(feature.type)}",
                    ]

                    q = feature.qualifiers
                    for key in ("locus_tag", "gene", "product", "note"):
                        if key in q and q[key]:
                            attrs.append(
                                f"{key}={esc('|'.join(str(x) for x in q[key]))}"
                            )

                    out.write(
                        f"{seqid}\tGenBank\t{feature.type}\t"
                        f"{start}\t{end}\t.\t{strand}\t.\t"
                        f"{';'.join(attrs)}\n"
                    )


if __name__ == "__main__":
    main()
