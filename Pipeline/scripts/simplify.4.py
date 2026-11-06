import pandas as pd
import argparse

def mutation_type(row):
    if any(len(alt) != len(row['REF']) for alt in row['ALT'].split(',')):
        return 'INDEL'
    else:
        return 'SNP'

def extract_imf(info):
    items = dict(item.split("=") for item in info.split(";") if '=' in item)
    return float(items.get('IMF', 0))

def calculate_reads_percentage(row):
    #if row['MutationType'] == 'INDEL':
    #    return round(extract_imf(row['INFO']) * 100)
    denominator = row['REF_F'] + row['REF_R'] + row['ALT_F'] + row['ALT_R']
    if denominator == 0:
        return 0
    return round((row['ALT_F'] + row['ALT_R']) / denominator * 100)

def parse_vcf(vcf_file, min_allele):
    try:
        # CHROM as text: a contig named 01 or NA must keep its name
        df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None, dtype={0: str}, keep_default_na=False)
    except pd.errors.EmptyDataError:
        columns = ['CHROM', 'POS', 'REF', 'ALT', '%Reads', 'MutationType']
        return pd.DataFrame(columns=columns)

    df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
    #df['DP4'] = df['INFO'].apply(lambda x: dict(item.split("=") for item in x.split(";") if '=' in item).get('DP4', '0,0,0,0'))
    df['SR'] = df['INFO'].apply(lambda x: dict(item.split("=") for item in x.split(";") if '=' in item).get('SR', '0,0,0,0'))
    # SR = ref fwd, ref rev, then fwd, rev for each ALT allele (a multi-allelic
    # record has more than 4 values): the ALT counts are summed over the alleles
    sr = df['SR'].apply(lambda v: [int(x) for x in v.split(',')])
    df['REF_F'] = sr.str[0]
    df['REF_R'] = sr.str[1]
    df['ALT_F'] = sr.apply(lambda v: sum(v[2::2]))
    df['ALT_R'] = sr.apply(lambda v: sum(v[3::2]))

    df['MutationType'] = df.apply(mutation_type, axis=1)
    df['%Reads'] = df.apply(calculate_reads_percentage, axis=1)

    df = df[df['%Reads'] > float(min_allele)]

    if df.empty:
        columns = ['CHROM', 'POS', 'REF', 'ALT', '%Reads', 'MutationType']
        return pd.DataFrame(columns=columns)

    df = df[['CHROM', 'POS', 'REF', 'ALT', '%Reads', 'MutationType']]
    return df

def main():
    parser = argparse.ArgumentParser(description='Simplify VCF file.')
    parser.add_argument('-in', dest='input_file', required=True, help='Input VCF file')
    parser.add_argument('-out', dest='output_file', required=True, help='Output simplified tsv file')
    parser.add_argument('-min', dest='min_allele', help='minimum allelic frequency threshold %%, default:0', default=0)
    args = parser.parse_args()

    df = parse_vcf(args.input_file, args.min_allele)
    df.to_csv(args.output_file, sep='\t', index=False)

if __name__ == '__main__':
    main()
