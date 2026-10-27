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
    # Load annotations from the GenBank file
    annotations = []
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type != "source":  # Exclude the source feature
                description = get_annotation_description(feature)
                annotations.append({
                    'start': feature.location.start,
                    'end': feature.location.end,
                    'description': description
                })

    # Process mutations TSV
    with open(mutations_file) as mutations, open(output_file, 'w', newline='') as out_file:
        reader = csv.DictReader(mutations, delimiter='\t')
        fieldnames = reader.fieldnames + ['Annotation Descriptions']
        writer = csv.DictWriter(out_file, fieldnames=fieldnames, delimiter='\t')
        
        writer.writeheader()
        for row in reader:
            pos = int(row['POS'])
            # Find annotations that encompass the mutation
            overlapping_annotations = [annotation['description'] for annotation in annotations if annotation['start'] <= pos <= annotation['end']]
            
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
