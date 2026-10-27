#!/usr/bin/env python3
"""
simplify_bcftools.py -- bcftools VCF -> simple TSV (short-read counterpart
of simplify.4.py, same CLI and same output).

Input : VCF from `bcftools mpileup -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP
        | bcftools call -m -v --ploidy 1` (plain or .gz, one sample).
Output: TSV with columns CHROM POS REF ALT %Reads MutationType
        (identical header / formatting to simplify.4.py).

  %Reads       = round(100 * AD[called allele] / sum(AD)), an integer, like
                 simplify.4.py (which uses ALT reads / all reads).
  MutationType = INDEL if len(REF) != len(ALT) else SNP.
  -min         = minimum QUAL (records with QUAL < min are dropped; a missing
                 QUAL '.' is treated as passing). Default 0.

Records whose GT is reference (0 / 0/0 / ./.) are skipped. Multi-allelic
sites: one row per distinct non-reference allele named in GT (with ploidy 1
that is exactly one row: the allele bcftools called), each with its own
AD[allele]/sum(AD). Symbolic alleles (<*>, <NON_REF>, *) are skipped.
A VCF with no records gives a header-only output.
"""

import argparse
import csv
import gzip

COLUMNS = ['CHROM', 'POS', 'REF', 'ALT', '%Reads', 'MutationType']


def mutation_type(ref, alt):
    return 'INDEL' if len(ref) != len(alt) else 'SNP'


def open_text(path):
    with open(path, 'rb') as fh:
        magic = fh.read(2)
    if magic == b'\x1f\x8b':
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_info_free_int_list(text):
    out = []
    for x in text.split(','):
        try:
            out.append(int(x))
        except ValueError:
            out.append(0)
    return out


def called_alleles(gt, n_alt):
    """Distinct non-reference allele indices (1-based) in a GT string."""
    if gt is None:
        return list(range(1, n_alt + 1))      # no GT column: report every ALT
    alleles = []
    for a in gt.replace('|', '/').split('/'):
        if a.isdigit() and int(a) > 0 and int(a) not in alleles:
            alleles.append(int(a))
    return alleles


def percent_reads(ad, allele):
    total = sum(ad)
    if total == 0 or allele >= len(ad):
        return 0
    return round(ad[allele] / total * 100)


def parse_vcf(vcf_file, min_qual):
    rows = []
    with open_text(vcf_file) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\r\n').split('\t')
            if len(f) < 5:
                continue
            chrom, pos, ref, alts = f[0], f[1], f[3], f[4].split(',')
            qual = f[5] if len(f) > 5 else '.'
            if qual not in ('.', ''):
                try:
                    if float(qual) < min_qual:
                        continue
                except ValueError:
                    pass
            if alts == ['.']:
                continue

            gt, ad = None, None
            if len(f) >= 10:
                keys = f[8].split(':')
                vals = f[9].split(':')
                d = dict(zip(keys, vals))
                gt = d.get('GT')
                if 'AD' in d:
                    ad = parse_info_free_int_list(d['AD'])

            for allele in called_alleles(gt, len(alts)):
                if allele > len(alts):
                    continue
                alt = alts[allele - 1]
                if alt.startswith('<') or alt == '*':
                    continue
                pct = percent_reads(ad, allele) if ad else 0
                rows.append([chrom, int(pos), ref, alt, pct, mutation_type(ref, alt)])
    return rows


def main():
    parser = argparse.ArgumentParser(description='Simplify bcftools VCF file.')
    parser.add_argument('-in', dest='input_file', required=True, help='Input VCF file (plain or gz)')
    parser.add_argument('-out', dest='output_file', required=True, help='Output simplified tsv file')
    parser.add_argument('-min', dest='min_qual', default=0, help='minimum QUAL threshold, default:0')
    args = parser.parse_args()

    rows = parse_vcf(args.input_file, float(args.min_qual))
    with open(args.output_file, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t', lineterminator='\n')
        w.writerow(COLUMNS)
        w.writerows(rows)


if __name__ == '__main__':
    main()
