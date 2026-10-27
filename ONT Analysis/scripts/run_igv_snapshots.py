#!/usr/bin/env python3
"""
Create and execute IGV batch scripts for Chlor_cassette hotspots.

This script:
  1. Converts the sample GenBank annotation to GFF3.
  2. Creates an IGV genome JSON using the sample FASTA.
  3. Creates an IGV batch file loading:
       - all-read BAM
       - Chlor_cassette-read BAM
       - GenBank-derived GFF3 annotations
  4. Visits each hotspot +/- 500 bp (or the configured hotspot interval)
     and saves a PNG snapshot.
"""
import argparse
import json
import os
import shlex
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--hotspots", required=True)
    p.add_argument("--all-bam", required=True)
    p.add_argument("--chlor-bam", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--genbank", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--igv", required=True)
    p.add_argument("--memory", default="4g")
    args = p.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    reference = Path(args.reference).resolve()
    genbank = Path(args.genbank).resolve()
    all_bam = Path(args.all_bam).resolve()
    chlor_bam = Path(args.chlor_bam).resolve()

    # Verify the FASTA has an index. Snakemake's normal reference branch
    # should already create one, but fail clearly if it is missing.
    fai = Path(str(reference) + ".fai")
    if not fai.exists():
        raise SystemExit(
            f"Missing FASTA index: {fai}\n"
            "Create it with: samtools faidx <reference.fasta>"
        )

    # Import the conversion function locally so this script is standalone
    # and does not depend on Python path configuration.
    from Bio import SeqIO

    gff = outdir / f"{genbank.stem}.gff3"

    def esc(value):
        value = str(value)
        return (
            value.replace("%", "%25")
            .replace(";", "%3B")
            .replace("=", "%3D")
            .replace(",", "%2C")
        )

    def feature_name(feature):
        q = feature.qualifiers
        for key in ("gene", "locus_tag", "label", "name", "product"):
            if key in q and q[key]:
                return str(q[key][0])
        return feature.type

    with open(gff, "w") as out:
        out.write("##gff-version 3\n")
        for record in SeqIO.parse(str(genbank), "genbank"):
            seqid = record.id
            for i, feature in enumerate(record.features, start=1):
                if feature.type == "source" or feature.location is None:
                    continue

                for part_index, part in enumerate(feature.location.parts, start=1):
                    start = int(part.start) + 1
                    end = int(part.end)

                    strand = "."
                    if feature.location.strand == 1:
                        strand = "+"
                    elif feature.location.strand == -1:
                        strand = "-"

                    name = feature_name(feature)
                    fid = f"{feature.type}_{i}_{part_index}"

                    attrs = [
                        f"ID={esc(fid)}",
                        f"Name={esc(name)}",
                        f"gene_biotype={esc(feature.type)}",
                    ]

                    for key in ("locus_tag", "gene", "product", "note"):
                        q = feature.qualifiers
                        if key in q and q[key]:
                            attrs.append(
                                f"{key}={esc('|'.join(str(x) for x in q[key]))}"
                            )

                    out.write(
                        f"{seqid}\tGenBank\t{feature.type}\t"
                        f"{start}\t{end}\t.\t{strand}\t.\t"
                        f"{';'.join(attrs)}\n"
                    )

    # Read hotspots.
    hotspots = []
    with open(args.hotspots) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        required = {
            "hotspot", "chrom", "start", "end",
            "expanded_start", "expanded_end"
        }
        missing = required - set(header)
        if missing:
            raise SystemExit(
                f"Hotspot table is missing columns: {', '.join(sorted(missing))}"
            )

        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            hotspots.append({
                "id": fields[col["hotspot"]],
                "chrom": fields[col["chrom"]],
                "start": int(fields[col["start"]]),
                "end": int(fields[col["end"]]),
                "expanded_start": int(fields[col["expanded_start"]]),
                "expanded_end": int(fields[col["expanded_end"]]),
            })

    # If there are no hotspots, leave a clear marker and finish successfully.
    if not hotspots:
        (outdir / "NO_CHLOR_HOTSPOTS.txt").write_text(
            "No Chlor_cassette hotspots passed the configured minimum-read threshold.\n"
        )
        return

    # Build IGV genome JSON.
    genome_json = outdir / "igv_genome.json"

    records = list(SeqIO.parse(str(reference), "fasta"))
    chromosome_order = [r.id for r in records]

    genome = {
        "id": reference.stem,
        "name": reference.stem,
        "fastaURL": str(reference),
        "indexURL": str(fai),
        "chromosomeOrder": chromosome_order,
        "tracks": [
            {
                "name": "GenBank annotations",
                "url": str(gff),
                "format": "gff3",
            }
        ],
    }

    genome_json.write_text(json.dumps(genome, indent=2))

    batch = outdir / "igv_batch.txt"

    lines = [
        "new",
        f"genome {genome_json}",
        "setSleepInterval 1000",
        f"snapshotDirectory {outdir}",
        "maxPanelHeight 5000",
        f"load {all_bam}",
        "setSleepInterval 1000",
        f"load {chlor_bam}",
        "setSleepInterval 1000",
        f"load {gff}",
        "setSleepInterval 1500",
    ]

    for h in hotspots:
        locus = (
            f'{h["chrom"]}:{h["expanded_start"] + 1}-'
            f'{h["expanded_end"]}'
        )
        snapshot_name = (
            f'{h["id"]}_{h["chrom"]}_'
            f'{h["start"] + 1}-{h["end"]}.png'
        )

        lines.extend([
            f"goto {locus}",
            "setSleepInterval 1000",
            "sort position",
            f"snapshot {snapshot_name}",
        ])

    # Exit IGV after all snapshots have been generated so the subprocess
    # returns and the Snakemake rule can finish.
    lines.append("exit")
    batch.write_text("\n".join(lines) + "\n")

    igv_cmd = shlex.split(args.igv)
    if not igv_cmd:
        raise SystemExit("Empty --igv command.")

    # IGV's official launcher supports: igv.sh -b <batch file>.
    # Use xvfb-run when available so this also works on headless clusters.
    xvfb = subprocess.run(
        ["bash", "-lc", "command -v xvfb-run >/dev/null 2>&1"],
        check=False,
    ).returncode == 0

    if xvfb:
        command = [
            "xvfb-run", "-a",
            *igv_cmd,
            "-b", str(batch),
        ]
    else:
        command = [
            *igv_cmd,
            "-b", str(batch),
        ]

    print("Running IGV:")
    print(" ".join(shlex.quote(x) for x in command))
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
