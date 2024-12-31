import csv
import sys

def paf_to_bed(paf_file, bed_file):
    with open(paf_file, 'r') as paf, open(bed_file, 'w') as bed:
        reader = csv.reader(paf, delimiter='\t')
        writer = csv.writer(bed, delimiter='\t', lineterminator='\n')
        
        for row in reader:
            target_name = row[0]
            target_start = row[2]
            target_end = row[3]
            writer.writerow([target_name, target_start, target_end])

def filter_paf(paf_file, overlapping_file, output_file):
    overlapping_intervals = set()

    with open(overlapping_file, 'r') as overlaps:
        reader = csv.reader(overlaps, delimiter='\t')
        for row in reader:
            overlapping_intervals.add((row[0], row[1], row[2]))

    with open(paf_file, 'r') as paf, open(output_file, 'w') as output:
        reader = csv.reader(paf, delimiter='\t')
        writer = csv.writer(output, delimiter='\t', lineterminator='\n')

        for row in reader:
            target_name = row[0]
            target_start = row[2]
            target_end = row[3]
            if (target_name, target_start, target_end) in overlapping_intervals:
                writer.writerow(row)

# Usage
paf_file = sys.argv[1]
bed_file = sys.argv[2]
features_file = sys.argv[3]
overlapping_file = sys.argv[4]
output_file = sys.argv[5]

paf_to_bed(paf_file, bed_file)

# Use bedtools to find overlapping regions
import os
os.system(f"bedtools intersect -a {bed_file} -b {features_file} -v -f 0.8 > {overlapping_file}")

# Filter the original PAF file
filter_paf(paf_file, overlapping_file, output_file)
