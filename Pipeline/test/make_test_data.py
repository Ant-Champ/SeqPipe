#!/usr/bin/env python3
"""Generate the small synthetic dataset used by `./run.sh test` and `./run.sh setup`.

Reference (test/data/ref/testref.gb): a 20 kb circular "chromosome" with three
CDS + a 4 kb plasmid. The sequenced strain is the reference with 2 SNPs in
geneA, geneB (bp 9001-10200) replaced by a Kan cassette and geneC
(bp 15001-16800) by a Chlor cassette in reverse orientation, Keio-style. The
two cassettes share identical 150 bp ends, like real Keio cassettes. The run
should report: 2 SNPs in geneA, low coverage + breakpoints over geneB and
geneC, and two cassette insertions (Kan in geneB, Chlor in geneC).

Writes (all small, committed to the repository):
  data/ref/testref.gb           GenBank reference (2 records)
  data/cassettes/*.gb           TEST_Kan_Cassette, TEST_Chlor_Cassette
  data/panel/*.gb, *.fasta      plasmid panel: the strain's plasmid + an absent plasmid
  data/fastq/barcode01.fastq.gz ONT reads (~40x, 1 % errors)
  data/fastq/sr_R1/R2.fastq.gz  paired 75 bp short reads (~60x)
Only the python standard library is used; output is deterministic (seed 1).
"""
import gzip
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent / "data"
rng = random.Random(1)


def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def genbank(name, seq, feats, circular=True):
    """Minimal GenBank record that Biopython parses; feats = [(start1, end1, gene, product)]."""
    topo = "circular" if circular else "linear"
    lines = [f"LOCUS       {name:<16} {len(seq):>11} bp    DNA     {topo:<8} BCT 01-JAN-2026",
             f"DEFINITION  synthetic test sequence {name}.",
             f"ACCESSION   {name}", f"VERSION     {name}.1", "KEYWORDS    .",
             "SOURCE      synthetic", "  ORGANISM  synthetic", "            .",
             "FEATURES             Location/Qualifiers",
             f"     source          1..{len(seq)}", '                     /organism="synthetic"']
    for s, e, gene, prod in feats:
        for ftype in ("gene", "CDS"):
            lines += [f"     {ftype:<16}{s}..{e}", f'                     /gene="{gene}"']
            if ftype == "CDS":
                lines.append(f'                     /product="{prod}"')
    lines.append("ORIGIN")
    low = seq.lower()
    for i in range(0, len(low), 60):
        blocks = " ".join(low[j:j + 10] for j in range(i, min(i + 60, len(low)), 10))
        lines.append(f"{i + 1:>9} {blocks}")
    lines.append("//")
    return "\n".join(lines) + "\n"


def mutate(seq, rate):
    """Random substitutions / 1-bp indels at `rate` per base (ONT-like noise)."""
    out = []
    for b in seq:
        r = rng.random()
        if r < rate / 3:
            out.append(rng.choice("ACGT"))
        elif r < 2 * rate / 3:
            continue
        elif r < rate:
            out.append(b + rng.choice("ACGT"))
        else:
            out.append(b)
    return "".join(out)


def sample_circular(genome, start, length):
    g = genome * 2
    return g[start:start + length]


def main():
    chrom = rand_seq(20000)
    plasmid = rand_seq(4000)
    absent = rand_seq(3000)
    # Keio-like cassettes: identical 150 bp ends (FRT / priming sites), unique middles
    end5, end3 = rand_seq(150), rand_seq(150)
    kan = end5 + rand_seq(1200) + end3            # 1500 bp
    chlor = end5 + rand_seq(1000) + end3          # 1300 bp

    feats_chr = [(2001, 3500, "geneA", "test protein A"),
                 (9001, 10200, "geneB", "test protein B (replaced by Kan cassette in strain)"),
                 (15001, 16800, "geneC", "test protein C (replaced by Chlor cassette in strain)")]
    feats_pl = [(501, 1300, "repA", "replication protein")]

    # strain = reference with 2 SNPs in geneA, geneB replaced by the Kan cassette
    # (forward) and geneC by the Chlor cassette (reverse orientation)
    strain = list(chrom)
    for pos in (2500, 3000):                       # 1-based positions
        strain[pos - 1] = "A" if strain[pos - 1] != "A" else "C"
    strain = "".join(strain)
    strain = strain[:9000] + kan + strain[10200:15000] + revcomp(chlor) + strain[16800:]
    replicons = [(strain, 1.0), (plasmid, 2.0)]    # plasmid at 2 copies / cell

    for d in ("ref", "panel", "cassettes", "fastq"):
        (HERE / d).mkdir(parents=True, exist_ok=True)

    with open(HERE / "ref" / "testref.gb", "w", newline="\n") as fh:
        fh.write(genbank("testchr", chrom, feats_chr))
        fh.write(genbank("testplasmid", plasmid, feats_pl))
    with open(HERE / "cassettes" / "TEST_Kan_Cassette.gb", "w", newline="\n") as fh:
        fh.write(genbank("TEST_Kan_Cassette", kan, [(151, 1350, "kan", "kanamycin resistance")], circular=False))
    with open(HERE / "cassettes" / "TEST_Chlor_Cassette.gb", "w", newline="\n") as fh:
        fh.write(genbank("TEST_Chlor_Cassette", chlor, [(151, 1150, "cat", "chloramphenicol resistance")], circular=False))
    # plasmid panel: the sample's plasmid (present, GenBank) + one absent plasmid (FASTA)
    with open(HERE / "panel" / "present_plasmid.gb", "w", newline="\n") as fh:
        fh.write(genbank("present_plasmid", plasmid, feats_pl))
    with open(HERE / "panel" / "absent_plasmid.fasta", "w", newline="\n") as fh:
        fh.write(f">absent_plasmid\n{absent}\n")

    # ONT reads: lengths 1-6 kb, ~40x per copy, 1 % errors, both strands
    with gzip.open(HERE / "fastq" / "barcode01.fastq.gz", "wt", newline="\n") as fh:
        n = 0
        for seq, copies in replicons:
            target = 40 * len(seq) * copies
            done = 0
            while done < target:
                ln = min(rng.randint(1000, 6000), len(seq))
                read = mutate(sample_circular(seq, rng.randrange(len(seq)), ln), 0.01)
                if rng.random() < 0.5:
                    read = revcomp(read)
                n += 1
                done += len(read)
                fh.write(f"@ont_{n}\n{read}\n+\n{'5' * len(read)}\n")

    # Short reads: PE75, insert 250-400 bp, ~60x, 0.1 % substitutions
    with gzip.open(HERE / "fastq" / "sr_R1.fastq.gz", "wt", newline="\n") as f1, \
         gzip.open(HERE / "fastq" / "sr_R2.fastq.gz", "wt", newline="\n") as f2:
        n = 0
        for seq, copies in replicons:
            pairs = int(60 * len(seq) * copies / 150)
            for _ in range(pairs):
                ins = rng.randint(250, 400)
                frag = sample_circular(seq, rng.randrange(len(seq)), ins)
                if rng.random() < 0.5:
                    frag = revcomp(frag)
                r1 = mutate(frag[:75], 0.001)[:75]
                r2 = mutate(revcomp(frag)[:75], 0.001)[:75]
                n += 1
                f1.write(f"@sr_{n}/1\n{r1}\n+\n{'I' * len(r1)}\n")
                f2.write(f"@sr_{n}/2\n{r2}\n+\n{'I' * len(r2)}\n")
    print(f"test data written to {HERE}")


if __name__ == "__main__":
    main()
