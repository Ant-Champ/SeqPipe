import numpy as np
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Look for coverage anomalies (outside 10 std_dev or depth of 0)')
parser.add_argument('-cov', dest='coverage_file', required=True, help='Input 1bp resolution coverage file REF POS DEPTH')
parser.add_argument('-ref', dest='ref_file', required=True, help='Reference fasta file')
parser.add_argument('-l', '--min_length', type=int, default=0, help='Minimum length of anomalous region to report')
parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
args = parser.parse_args()

ref_sequence = list(SeqIO.parse(args.ref_file, "fasta"))[0]

# Load coverage data
with open(args.coverage_file, 'r') as f:
    lines = f.readlines()

if not lines:
    print("Error: The coverage file is empty. This means that no reads could be aligned on the reference.")
    exit()

# Split into position and depth
references, positions, depths = zip(*(line.split('\t') for line in lines))

# Convert strings to integers
positions = [int(pos) for pos in positions]
depths = [int(depth) for depth in depths]

# Calculate mean coverage and standard deviation
mean_coverage = np.mean(depths)
std_dev = np.std(depths)

if args.verbose:
    print(f"Mean coverage: {round(mean_coverage)}")
    print(f"Standard deviation: {round(std_dev)}")

# Define ranges
lower_limit = max(0,mean_coverage - 1 * std_dev)
upper_limit = mean_coverage + 1 * std_dev

anomaly_detected = False
anomaly_start = None
for ref, pos, depth in zip(references, positions, depths):
    if depth < lower_limit or depth > upper_limit or depth == 0:   # Added check for depth of 0 here
        if anomaly_start is None:
            anomaly_start = pos
    else:
        if anomaly_start is not None and pos - anomaly_start >= args.min_length:
            print(f"Anomalous region from position {anomaly_start} to {pos-1} on reference {ref} with length {pos - anomaly_start}")
            anomaly_detected=True
        anomaly_start = None

# Check if the last position was part of an anomaly
if anomaly_start is not None and positions[-1] - anomaly_start >= args.min_length:
    print(f"Anomalous region from position {anomaly_start} to {positions[-1]} on reference {references[-1]} with length {positions[-1] - anomaly_start}")
    anomaly_detected=True
