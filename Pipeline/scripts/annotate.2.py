import argparse
from Bio import SeqIO
import csv

def get_annotation_description(feature):
    # This function returns a string with selected qualifiers for the feature
    descriptions = []
    if 'gene' in feature.qualifiers:
        descriptions.append("gene:" + feature.qualifiers['gene'][0])
    if 'product' in feature.qualifiers:
        descriptions.append("product:" + feature.qualifiers['product'][0])
    if 'label' in feature.qualifiers:
        descriptions.append("label:" + feature.qualifiers['label'][0])
    return '; '.join(descriptions) if descriptions else 'No specific annotation'

def mutation_in_annotation(mutations_file, gb_file, output_file):
    # Load annotations from the GenBank file, per record id (= CHROM)
    annotations = {}
    for record in SeqIO.parse(gb_file, "genbank"):
        feats = annotations.setdefault(record.id, [])
        for feature in record.features:
            if feature.type != "source":  # Exclude the source feature
                description = get_annotation_description(feature)
                # one interval per part: a join() or a feature across the origin
                # is not one block from its lowest to its highest base
                for part in feature.location.parts:
                    feats.append({
                        'start': int(part.start),  # 0-based
                        'end': int(part.end),      # exclusive
                        'description': description
                    })
    single = next(iter(annotations.values())) if len(annotations) == 1 else None

    # Process mutations TSV
    with open(mutations_file) as mutations, open(output_file, 'w', newline='') as out_file:
        reader = csv.DictReader(mutations, delimiter='\t')
        if reader.fieldnames is None:  # completely empty input
            return
        fieldnames = reader.fieldnames + ['Annotation Descriptions']
        writer = csv.DictWriter(out_file, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()
        for row in reader:
            pos = int(row['POS'])
            feats = annotations.get(row.get('CHROM'), single)
            if feats is None:
                feats = []
            # Find annotations that encompass the mutation (1-based POS vs 0-based start)
            overlapping_annotations = list(dict.fromkeys(
                annotation['description'] for annotation in feats if annotation['start'] < pos <= annotation['end']))

            # Join overlapping annotations into a single string
            row['Annotation Descriptions'] = '; '.join(overlapping_annotations) if overlapping_annotations else 'None'
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Check if mutations are in annotated regions of a GenBank file and list descriptive annotations.")
    parser.add_argument('-t', '--tsv', required=True, help="Path to the input TSV file containing mutations.")
    parser.add_argument('-a', '--annotation', required=True, help="Path to the GenBank file with annotations.")
    parser.add_argument('-o', '--output', required=True, help="Path to the output TSV file.")

    args = parser.parse_args()

    mutation_in_annotation(args.tsv, args.annotation, args.output)

if __name__ == "__main__":
    main()
