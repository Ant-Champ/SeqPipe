#!/usr/bin/env python3
"""
Create and execute IGV batch scripts for cassette insertion sites.

This script:
  1. Converts the sample GenBank annotation to GFF3.
  2. Creates an IGV genome JSON using the sample FASTA.
  3. Creates an IGV batch file loading:
       - all-read BAM (sample vs genome)
       - cassette BAM (sample vs genome + cassette contigs)
       - GenBank-derived GFF3 annotations
  4. Visits each insertion region (--regions = the insertions table of
     find_cassette_insertions.py, '#' header skipped, first three columns
     contig/start/end 0-based half-open, status = column 5, cassette = column
     6; an extra trailing 'features' column is ignored) +/- --flank bp and
     saves a PNG snapshot per region.

If there are no regions, a NO_CASSETTE_INSERTIONS.txt marker is written
instead and IGV is not run.
"""
import argparse
import json
import os
import random
import re
import shlex
import subprocess
import time
from pathlib import Path

from Bio import SeqIO


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--regions", required=True,
                   help="insertions TSV from find_cassette_insertions.py")
    p.add_argument("--flank", type=int, default=1000,
                   help="bp added on each side of every region (default: 1000)")
    p.add_argument("--all-bam", required=True,
                   help="BAM of the sample mapped to its genome")
    p.add_argument("--cassette-bam", required=True,
                   help="BAM of the sample mapped to genome + cassette contigs")
    p.add_argument("--reference", required=True)
    p.add_argument("--genbank", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--igv", required=True)
    p.add_argument("--label", default="Cassette",
                   help="name used in messages and the "
                        "NO_CASSETTE_INSERTIONS.txt marker (default: Cassette)")
    p.add_argument("--memory", type=int, default=4000,
                   help="IGV Java heap in MB, applied via _JAVA_OPTIONS=-Xmx<memory>m "
                        "(default: 4000)")
    p.add_argument("--attempts", type=int, default=3,
                   help="tries before giving up (IGV_FAILED.txt, exit 0) (default: 3)")
    args = p.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    reference = Path(args.reference).resolve()
    genbank = Path(args.genbank).resolve()
    all_bam = Path(args.all_bam).resolve()
    cassette_bam = Path(args.cassette_bam).resolve()

    # Verify the FASTA has an index. Snakemake's normal reference branch
    # should already create one, but fail clearly if it is missing.
    fai = Path(str(reference) + ".fai")
    if not fai.exists():
        raise SystemExit(
            f"Missing FASTA index: {fai}\n"
            "Create it with: samtools faidx <reference.fasta>"
        )

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

    # Read regions (insertions table: contig start end length status cassette ...).
    chrom_len = {r.id: len(r.seq) for r in SeqIO.parse(str(reference), "fasta")}
    regions = []
    with open(args.regions) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise SystemExit(
                    f"Regions table needs >= 6 columns (contig start end length status "
                    f"cassette), got {len(fields)}: {line.strip()}"
                )
            start, end = int(fields[1]), int(fields[2])
            chrom = fields[0]
            top = chrom_len.get(chrom, end + args.flank)
            regions.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "label": re.sub(r"[^A-Za-z0-9._-]+", "_", f"{fields[4]}_{fields[5]}").strip("_"),
                "expanded_start": max(0, start - args.flank),
                "expanded_end": min(top, end + args.flank),
            })

    # If there are no regions, leave a clear marker and finish successfully.
    if not regions:
        (outdir / "NO_CASSETTE_INSERTIONS.txt").write_text(
            f"No {args.label} insertion sites were found.\n"
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
        # The GFF3 is loaded once, via the batch file below (not as a
        # genome track as well, which would show it twice).
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
        f"load {cassette_bam}",
        "setSleepInterval 1000",
        f"load {gff}",
        "setSleepInterval 1500",
    ]

    for i, h in enumerate(regions, start=1):
        locus = (
            f'{h["chrom"]}:{h["expanded_start"] + 1}-'
            f'{h["expanded_end"]}'
        )
        chrom_tag = re.sub(r"[^A-Za-z0-9._-]+", "_", h["chrom"])
        snapshot_name = (
            f'{i:02d}_{h["label"]}_{chrom_tag}_'
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

    # Cap the IGV JVM heap (the igv.sh launcher honours _JAVA_OPTIONS).
    env = dict(os.environ)
    env["_JAVA_OPTIONS"] = f"-Xmx{args.memory}m"

    # Many IGV jobs run at once: `xvfb-run -a` alone lets two of them pick the
    # same free display and one dies. Start each attempt at a random display
    # number, and retry a few times before giving up.
    log_tail = ""
    for attempt in range(1, args.attempts + 1):
        for png in outdir.glob("*.png"):
            png.unlink()
        if xvfb:
            command = ["xvfb-run", "-a", "-n", str(random.randint(100, 30000)),
                       "-s", "-screen 0 1600x1200x24", *igv_cmd, "-b", str(batch)]
        else:
            command = [*igv_cmd, "-b", str(batch)]
        print(f"Running IGV (attempt {attempt}/{args.attempts}):", flush=True)
        print(" ".join(shlex.quote(x) for x in command), flush=True)
        res = subprocess.run(command, env=env, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True)
        print(res.stdout, flush=True)
        n_png = len(list(outdir.glob("*.png")))
        if res.returncode == 0 and n_png >= len(regions):
            print(f"{n_png} snapshot(s) written")
            (outdir / "IGV_FAILED.txt").unlink(missing_ok=True)
            return
        log_tail = "\n".join(res.stdout.splitlines()[-30:])
        print(f"IGV exit status {res.returncode}, {n_png}/{len(regions)} snapshot(s)", flush=True)
        time.sleep(random.uniform(5, 20))

    # Pictures are a nice-to-have: do not fail the whole run for them.
    (outdir / "IGV_FAILED.txt").write_text(
        f"IGV failed {args.attempts} times; the insertion tables are still valid. "
        f"See logs/cassette_igv_snapshots.log. Last output:\n{log_tail}\n")
    print("WARNING: IGV failed, IGV_FAILED.txt written")


if __name__ == "__main__":
    main()
